#!/usr/bin/env Rscript
# Native PacBio DADA2. Input FASTQs have already been trimmed/oriented by cutadapt.
args <- commandArgs(trailingOnly = TRUE)
# TSV counts must remain plain integers even for deeply sequenced samples.
options(scipen = 999)
if (!requireNamespace("dada2", quietly = TRUE)) stop("Install the R dada2 package with Bioconductor")
ns <- asNamespace("dada2")
if (!exists("PacBioErrfun", envir = ns, inherits = FALSE)) {
    stop("This dada2 version does not provide PacBioErrfun; install a compatible version")
}
cat(R.version.string, "\nDADA2:", as.character(utils::packageVersion("dada2")), "\n")
if (identical(args, "--check")) quit(status = 0)
if (length(args) != 12L) stop("Expected manifest, temp, fasta, table, stats, filtered, transitions, minLen, maxLen, maxEE, abundance, threads")
manifest <- read.delim(args[1], check.names = FALSE, colClasses = "character", na.strings = NULL)
writeLines(c(R.version.string, paste("DADA2:", as.character(utils::packageVersion("dada2")))),
           file.path(dirname(args[4]), "post_dada2_versions.txt"))
temp <- args[2]
minLen <- as.integer(args[8]); maxLen <- as.integer(args[9])
maxEE <- as.numeric(args[10]); abundance <- as.numeric(args[11])
threads <- as.integer(args[12])
multi <- if (threads == 1L) FALSE else threads
filter.multi <- if (.Platform$OS.type == "windows") FALSE else multi
dir.create(file.path(temp, "dada2_filtered"), showWarnings = FALSE)
filtered.paths <- file.path(temp, "dada2_filtered", sprintf("sample_%06d.fastq", seq_len(nrow(manifest))))
stats <- data.frame(sample = manifest$sample, input = 0L, filtered = 0L,
                    denoised = 0L, non_chimeric = 0L, check.names = FALSE)
for (i in seq_len(nrow(manifest))) {
    if (!nzchar(manifest$fastq[i])) next
    out <- dada2::filterAndTrim(manifest$fastq[i], filtered.paths[i],
        minLen = minLen, maxLen = maxLen, truncLen = 0, trimLeft = 0,
        maxEE = maxEE, truncQ = 2, minQ = 3, maxN = 0, rm.phix = FALSE,
        compress = FALSE, multithread = filter.multi)
    stats$input[i] <- out[1, "reads.in"]
    stats$filtered[i] <- out[1, "reads.out"]
}
active <- which(stats$filtered > 0)
if (!length(active)) stop("No reads passed DADA2 filtering")
err <- dada2::learnErrors(filtered.paths[active], nreads = 1000000,
    errorEstimationFunction = get("PacBioErrfun", envir = ns),
    multithread = multi, BAND_SIZE = 32)
saveRDS(err, file.path(temp, "dada2_error_model.rds"))
dds <- vector("list", length(active))
names(dds) <- manifest$sample[active]
for (j in seq_along(active)) {
    drp <- dada2::derepFastq(filtered.paths[active[j]])
    dds[[j]] <- dada2::dada(drp, err = err, pool = FALSE,
                          multithread = multi, BAND_SIZE = 32)
    stats$denoised[active[j]] <- sum(dada2::getUniques(dds[[j]]))
}
seqtab <- dada2::makeSequenceTable(dds)
nochim <- dada2::removeBimeraDenovo(seqtab, method = "consensus",
    minFoldParentOverAbundance = abundance, allowOneOff = FALSE,
    multithread = multi)
if (!ncol(nochim) || !sum(nochim)) stop("No ASVs remain after DADA2 chimera removal")
# Lexical sequence order keeps IDs stable for the same set of ASVs.
seqs <- sort(colnames(nochim))
ids <- paste0("ASV_", seq_along(seqs))
full <- matrix(0L, nrow = nrow(manifest), ncol = length(seqs),
               dimnames = list(manifest$sample, seqs))
full[rownames(nochim), ] <- nochim[, seqs, drop = FALSE]
stats$non_chimeric <- rowSums(full)
writeLines(as.vector(rbind(paste0(">", ids), seqs)), args[3])
table <- data.frame(feature = ids, t(full), check.names = FALSE)
names(table)[1] <- "#ASV ID"
write.table(table, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(stats, args[5], sep = "\t", quote = FALSE, row.names = FALSE)
# Preserve filtered reads with their original headers and quality strings for splitting.
dest <- file(args[6], "wb")
for (path in filtered.paths[active]) {
    src <- file(path, "rb")
    repeat {
        chunk <- readBin(src, "raw", n = 1048576L)
        if (!length(chunk)) break
        writeBin(chunk, dest)
    }
    close(src)
}
close(dest)
errors <- dada2::getErrors(err, detailed = TRUE)
trans <- errors$trans
estimated <- errors$err_out
if (is.null(trans) || is.null(estimated)) stop("DADA2 error model lacks transition statistics")
grid <- expand.grid(transition = rownames(trans), quality = colnames(trans),
                    stringsAsFactors = FALSE)
grid$count <- as.vector(trans)
grid$estimated <- as.vector(estimated[rownames(trans), colnames(trans), drop = FALSE])
denom <- tapply(grid$count, list(substr(grid$transition, 1, 1), grid$quality), sum)
totals <- mapply(function(tr, q) denom[substr(tr, 1, 1), q], grid$transition, grid$quality)
grid$observed <- ifelse(totals > 0, grid$count / totals, 0)
grid$nominal <- ifelse(substr(grid$transition, 1, 1) == substr(grid$transition, 3, 3),
    1 - 10^(-as.numeric(grid$quality) / 10), 10^(-as.numeric(grid$quality) / 10) / 3)
write.table(grid, args[7], sep = "\t", quote = FALSE, row.names = FALSE)
cat("Completed native PacBio DADA2\n")
