import logging
import os
from typing import List, Tuple

from gtdb_release_tk.common import summarise_file
from gtdb_release_tk.files.taxonomy import TaxonomyFile
from gtdb_release_tk.models.taxonomy_rank import TaxonomyRank, RankEnum

logger = logging.getLogger('timestamp')


class GTDBDictFile(object):
    NAME = 'gtdb_r{release}.dic'

    def __init__(self, data: frozenset):
        self.data: frozenset = data

    @classmethod
    def create(cls, taxonomy_file: TaxonomyFile):
        logger.info(f'Creating dictionary from {len(taxonomy_file.data):,} genomes')
        taxa = list()
        for gid, tax in taxonomy_file.data.items():
            for rank in tax.ranks:
                taxa.extend(GTDBDictFile.parse_rank(rank))
                if rank.type is RankEnum.SPECIES:
                    taxa.extend(GTDBDictFile.parse_rank_species(rank))
        taxa_set = frozenset(taxa)
        logger.info(f'Created {len(taxa_set):,} unique definitions (from {len(taxa):,})')
        return cls(taxa_set)

    def write(self, root_dir: str, release: str):
        path = os.path.join(root_dir, self.NAME.format(release=release))
        logger.info(f'Writing {len(self.data):,} definitions to disk: {path}')
        with open(path, 'w', encoding='utf-8') as f:
            for name in sorted(self.data):
                f.write(f'{name}\n')
        logger.info(summarise_file(path))

    @staticmethod
    def parse_rank(rank: TaxonomyRank) -> Tuple[str, str]:
        taxa = list()
        taxa.append(rank.name)  # e.g. Archaea
        taxa.append(rank.full)  # e.g. d__Archaea
        return rank.name, rank.full

    @staticmethod
    def parse_rank_species(rank: TaxonomyRank) -> List[str]:
        taxa = list()
        generic, specific = rank.name.split()
        generic, specific = generic.strip(), specific.strip()
        taxa.append(generic)
        taxa.append(specific)
        taxa.append(f'{RankEnum.SPECIES.value}{generic}')
        taxa.append(f'{generic[0]}. {specific}')
        taxa.append(f'{RankEnum.SPECIES.value}{generic[0]}. {specific}')
        return taxa
