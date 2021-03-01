from gtdb_release_tk.common import assert_file_exists

from typing import Tuple, Collection, Dict


class TaxonomyString(object):
    __slots__ = ('d', 'p', 'c', 'o', 'f', 'g', 's')

    def __init__(self, tax_string: str):
        ranks = tax_string.replace('; ', ';').split(';')
        assert(len(ranks) == 7)
        self.d: str = ranks[0]
        self.p: str = ranks[1]
        self.c: str = ranks[2]
        self.o: str = ranks[3]
        self.f: str = ranks[4]
        self.g: str = ranks[5]
        self.s: str = ranks[6]

    def __repr__(self):
        return f'{self.d};{self.p};{self.c};{self.o};{self.f};{self.g};{self.s}'

    @staticmethod
    def empty(tax: str):
        return len(tax) <= 3



class TaxonomyFile(object):

    def __init__(self, path: str, data: Dict[str, TaxonomyString]):
        self.path: str = path
        self.data: Dict[str, TaxonomyString] = data

    def __repr__(self):
        return f'{self.path[-50:]} ({len(self.data):,} genomes)'

    @classmethod
    def read(cls, path: str):
        assert_file_exists(path)
        data = dict()
        with open(path, 'r') as f:
            for line in f.readlines():
                gid, tax = line.rstrip().split('\t')
                data[gid] = TaxonomyString(tax)
        return cls(path, data)


