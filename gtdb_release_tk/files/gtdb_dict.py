import os
from typing import List

from gtdb_release_tk.files.taxonomy import TaxonomyFile


class GTDBDictFile(object):
    NAME = 'gtdb_r{release}.dic'

    def __init__(self, data: frozenset):
        self.data: frozenset = data

    @classmethod
    def create(cls, taxonomy_file: TaxonomyFile):
        taxa = list()
        for gid, tax in taxonomy_file.data.items():
            for rank in tax.ranks:
                taxa.extend(GTDBDictFile.parse_rank(rank))
        return cls(frozenset(taxa))

    def write(self, root_dir: str, release: str):
        path = os.path.join(root_dir, self.NAME.format(release=release))
        with open(path, 'w', encoding='utf-8') as f:
            for name in sorted(self.data):
                f.write(f'{name}\n')

    @staticmethod
    def parse_rank(rank: str) -> List[str]:
        taxa = list()
        taxa.append(rank[3:])  # without prefix
        taxa.append(rank)      # full rank

        if rank.startswith('s__'):
            taxa.extend(GTDBDictFile.parse_rank_species(rank))
        return taxa

    @staticmethod
    def parse_rank_species(rank: str) -> List[str]:
        taxa = list()
        generic, specific = rank[3:].split()
        taxa.append(generic.strip())
        taxa.append(specific.strip())
        taxa.append(f's__{generic.strip()}')
        taxa.append(f'{generic[0]}. {specific}')
        taxa.append(f's__{generic[0]}. {specific}')
        return taxa
