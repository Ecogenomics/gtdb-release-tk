from enum import Enum
from typing import Optional

from gtdb_release_tk.exceptions import InvalidTaxonomyRank


class RankEnum(Enum):
    DOMAIN = 'd__'
    PHYLUM = 'p__'
    CLASS = 'c__'
    ORDER = 'o__'
    FAMILY = 'f__'
    GENUS = 'g__'
    SPECIES = 's__'


VALID_PREFIX = {x.value for x in RankEnum}


class TaxonomyRank(object):
    """Represents a rank (e.g. d__Bacteria)"""
    __slots__ = ('rank', 'name', 'type')

    def __init__(self, rank: str):
        """Create the rank (rank=d__Archaea, name=Archaea, type=RankEnum)"""
        if len(rank) < 3 or rank[0:3] not in VALID_PREFIX:
            raise InvalidTaxonomyRank(f'Not a valid rank: {rank}')
        self.rank: str = rank
        self.name: Optional[str] = rank[3:] if len(rank) > 3 else None
        self.type: RankEnum = RankEnum(rank[0:3])

    def __eq__(self, other):
        return isinstance(other, TaxonomyRank) and self.rank == other.rank and \
               self.name == other.name and self.type == other.type

    def __lt__(self, other):
        return isinstance(other, TaxonomyRank) and self.rank < other.rank

    def __hash__(self):
        return hash(self.rank)

    def __str__(self):
        return self.rank

    def __repr__(self):
        return self.rank

    def is_empty(self) -> bool:
        """True if the rank is empty (e.g. `d__`), False otherwise"""
        return self.name is None
