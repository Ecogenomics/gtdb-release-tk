import re
from typing import Tuple

from gtdb_release_tk.exceptions import InvalidTaxonomy

re_tax_string = re.compile(r'(d__.*); ?(p__.*); ?(c__.*); ?(o__.*); ?(f__.*); ?(g__.*); ?(s__.*)')


class Taxonomy(object):
    __slots__ = ('d', 'p', 'c', 'o', 'f', 'g', 's')

    def __init__(self, tax_string: str):
        ranks = self.str_to_ranks(tax_string)
        self.d: str = ranks[0]
        self.p: str = ranks[1]
        self.c: str = ranks[2]
        self.o: str = ranks[3]
        self.f: str = ranks[4]
        self.g: str = ranks[5]
        self.s: str = ranks[6]

    @property
    def ranks(self) -> Tuple[str, str, str, str, str, str, str]:
        return self.d, self.p, self.c, self.o, self.f, self.g, self.s

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f'{self.d};{self.p};{self.c};{self.o};{self.f};{self.g};{self.s}'

    def __eq__(self, other):
        return isinstance(other, Taxonomy) and str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    @staticmethod
    def empty(tax: str):
        return len(tax) <= 3

    @staticmethod
    def str_to_ranks(tax_string) -> Tuple[str, str, str, str, str, str, str]:
        """Converts a 7 rank tax string to a tuple and verifies completeness."""
        ranks = re_tax_string.search(tax_string)
        if ranks is None or len(ranks.groups()) != 7:
            raise InvalidTaxonomy(f'Invalid taxonomy string: {tax_string}')
        return (ranks.group(1), ranks.group(2), ranks.group(3), ranks.group(4),
                ranks.group(5), ranks.group(6), ranks.group(7))
