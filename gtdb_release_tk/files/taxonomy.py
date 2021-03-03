import logging
from typing import Dict

from gtdb_release_tk.common import assert_file_exists, assert_dir_exists
from gtdb_release_tk.files.metadata import MetadataFile
from gtdb_release_tk.files.sp_clusters import SpClustersFile
from gtdb_release_tk.models.taxonomy_string import TaxonomyString
import os
logger = logging.getLogger('timestamp')


class TaxonomyFile(object):

    def __init__(self, data: Dict[str, TaxonomyString]):
        self.data: Dict[str, TaxonomyString] = data

    def __repr__(self):
        return f'{len(self.data):,} genomes'

    def __eq__(self, other):
        return isinstance(other, TaxonomyFile) and self.as_string() == other.as_string()

    def __hash__(self):
        return hash(self.as_string())

    def as_string(self):
        out = list()
        for gid, tax in sorted(self.data.items()):
            out.append(f'{gid}\t{tax}')
        return '\n'.join(out)


    @classmethod
    def read(cls, path: str):
        assert_file_exists(path)
        data = dict()
        with open(path, 'r') as f:
            for line in f.readlines():
                gid, tax = line.rstrip().split('\t')
                data[gid] = TaxonomyString(tax)
        return cls(data)

    def write(self, path):
        with open(path, 'w') as f:
            for gid, tax in self.data.items():
                f.write(f'{gid}\t{tax}\n')

    @classmethod
    def create(cls, domain: str, metadata_file: MetadataFile, sp_clusters_file: SpClustersFile):
        """Generate taxonomy files for GTDB website."""
        logger.info('Parsing all species clusters from GTDB metadata.')
        genome_rep_id = sp_clusters_file.genome_rep_id()

        logger.info(' ...identified {:,} representative genomes spanning {:,} genomes.'.format(
            len(sp_clusters_file.rows),
            len(genome_rep_id)))

        # get representative genome for each GTDB species cluster
        logger.info('Parsing representative genomes from GTDB metadata.')
        reps = metadata_file.get_rep_tax()  # This just gets the species reps and their tax
        logger.info(' ...identified {:,} representative genomes ({:,} bacterial, {:,} archaeal).'.format(
            sum([1 for t in reps.values()]),
            sum([1 for t in reps.values() if t.d == 'd__Bacteria']),
            sum([1 for t in reps.values() if t.d == 'd__Archaea'])))

        # Create the data
        out = dict()
        for gid in genome_rep_id:
            assert (not gid.startswith('U_'))

            rep_id = genome_rep_id[gid]
            taxa = reps[rep_id]
            assert (taxa.d in ('d__Bacteria', 'd__Archaea'))

            if taxa.d == domain:
                out[gid] = taxa

        return cls(out)


class ArcTaxonomyFile(TaxonomyFile):

    @classmethod
    def create(cls, metadata_file: MetadataFile, sp_clusters_file: SpClustersFile):
        return super().create('d__Archaea', metadata_file, sp_clusters_file)


class BacTaxonomyFile(TaxonomyFile):

    @classmethod
    def create(cls, metadata_file: MetadataFile, sp_clusters_file: SpClustersFile):
        return super().create('d__Bacteria', metadata_file, sp_clusters_file)
