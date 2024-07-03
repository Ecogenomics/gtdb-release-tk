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


VALID_PREFIX = frozenset({x.value for x in RankEnum})

RANK_ORDERING = (RankEnum.DOMAIN, RankEnum.PHYLUM, RankEnum.CLASS, RankEnum.ORDER,
                 RankEnum.FAMILY, RankEnum.GENUS, RankEnum.SPECIES)


class TaxonomyRank(object):
    """Represents a rank (e.g. d__Bacteria)"""
    __slots__ = ('rank', 'type')

    def __init__(self, rank: str):
        """Create the rank (rank=d__Archaea, name=Archaea, type=RankEnum)"""
        if len(rank) < 3 or rank[0:3] not in VALID_PREFIX:
            raise InvalidTaxonomyRank(f'Not a valid rank: {rank}')
        self.rank: str = rank
        self.type: RankEnum = RankEnum(rank[0:3])

    @property
    def name(self) -> Optional[str]:
        """e.g. d__Archaea -> Archaea"""
        return self.rank[3:] if len(self.rank) > 3 else None

    @property
    def canonical(self) -> str:
        """e.g. p__Firmicutes_B -> p__Firmicutes"""
        return self.rank_to_canonical(self.rank)

    @staticmethod
    def rank_to_canonical(rank: str) -> str:
        return rank[0:rank.rfind('_')] if len(rank) > 3 and '_' in rank[3:] else rank

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
