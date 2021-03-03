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
    __slots__ = ('full', 'name', 'type')

    def __init__(self, rank: str):
        """Create the rank (full=d__Archaea, name=Archaea, type=RankEnum)"""
        if len(rank) < 3 or rank[0:3] not in VALID_PREFIX:
            raise InvalidTaxonomyRank(f'Not a valid rank: {rank}')
        self.full: str = rank
        self.name: Optional[str] = rank[3:] if len(rank) > 3 else None
        self.type: RankEnum = RankEnum(rank[0:3])

    def __eq__(self, other):
        return isinstance(other, TaxonomyRank) and self.full == other.full and \
               self.name == other.name and self.type == other.type

    def __hash__(self):
        return hash(self.full)

    def __str__(self):
        return self.full

    def __repr__(self):
        return f'{self.full} ({self.type.name})'

    def is_empty(self) -> bool:
        """True if the rank is empty (e.g. `d__`), False otherwise"""
        return self.name is None
