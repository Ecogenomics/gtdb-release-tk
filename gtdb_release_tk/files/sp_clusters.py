
from typing import Dict

from gtdb_release_tk.common import optional_float, optional_int


class SpClustersFile(object):

    def __init__(self, header, rows):
        self.header = header
        self.rows = rows

    @classmethod
    def read(cls, path):
        rows = dict()
        with open(path) as f:
            header = f.readline().strip().split('\t')
            header = {x: i for i, x in enumerate(header)}
            for line in f.readlines():
                row = SpClustersRow(header, line.strip().split('\t'))
                rows[row.rep_genome] = row
        return cls(header, rows)

    def __repr__(self):
        return f'{len(self.rows):,} rows'

    def genome_rep_id(self):
        out = dict()
        for rep_id, sp_cluster_row in self.rows.items():
            for gid in sp_cluster_row.clustered_genomes:
                out[gid] = rep_id
        return out

class SpClustersRow(object):

    __slots__ = ('rep_genome', 'gtdb_species', 'gtdb_taxonomy', 'ani_radius',
                 'mean_intra_species_ani', 'min_intra_species_ani', 'mean_intra_species_af',
                 'min_intra_species_af', 'no_clustered_genomes', 'clustered_genomes')

    def __init__(self, header: Dict[str, int], row: list):
        self.rep_genome = row[header['Representative genome']]
        self.gtdb_species = row[header['GTDB species']]
        self.gtdb_taxonomy = row[header['GTDB taxonomy']]
        self.ani_radius = optional_float(row[header['ANI circumscription radius']])
        self.mean_intra_species_ani = optional_float(row[header['Mean intra-species ANI']])
        self.min_intra_species_ani = optional_float(row[header['Min intra-species ANI']])
        self.mean_intra_species_af = optional_float(row[header['Mean intra-species AF']])
        self.min_intra_species_af = optional_float(row[header['Min intra-species AF']])
        self.no_clustered_genomes = optional_int(row[header['No. clustered genomes']])
        self.clustered_genomes = set(row[header['Clustered genomes']].split(','))

        assert(len(self.clustered_genomes) == self.no_clustered_genomes)

    def __repr__(self):
        return f'{self.rep_genome} ({self.no_clustered_genomes:,} clustered genomes)'
