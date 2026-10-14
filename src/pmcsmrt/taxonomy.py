"""Parse each SINTAX rank independently, including partial assignments."""
import re

RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
TAX_COLUMNS = [field for rank in RANKS for field in (rank, rank + '_conf')]


def taxonomy_fields(assignment: str) -> list[str]:
    values = {rank: (name, confidence) for rank, name, confidence in
              re.findall(r'(?:^|,)([dkpcofgs]):([^,()]+)\(([^)]+)\)', assignment)}
    # Some SINTAX databases use k instead of d for the top rank.
    return [value for rank in 'dpcofgs'
            for value in values.get(rank, values.get('k', ('NA', 'NA'))
                                    if rank == 'd' else ('NA', 'NA'))]
