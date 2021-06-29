###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2021'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import logging
import itertools
from pathlib import PurePath
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy
from biolib.common import canonical_gid

from gtdb_release_tk.lpsn import LPSN
from gtdb_release_tk.taxon_utils import is_latin_taxon, generic_name


class Synonyms(object):
    """Generate files indicating NCBI names considered synonyms under the GTDB taxonomy."""

    def __init__(self, output_dir):

        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def parse_untrustworthy_type_genomes(self, untrustworthy_type_material):
        """Parse file indicating genomes considered to be untrustworthy as type material in GTDB.

        Ideally this information would just be a field in the GTDB and could be read from the
        GTDB metadata file. This is a work in progress.
        """

        untrustworthy_type_gids = set()
        with open(untrustworthy_type_material) as f:
            f.readline()
            for line in f:
                tokens = line.strip().split('\t')
                untrustworthy_type_gids.add(canonical_gid(tokens[0]))

        return untrustworthy_type_gids

    def read_gtdb_metadata(self, gtdb_metadata_file, untrustworthy_type_gids):
        """Read metadata for genomes."""

        Metadata = namedtuple('Metadata', ['gtdb_is_rep',
                                           'gtdb_genome_rep',
                                           'gtdb_type_strain',
                                           'gtdb_type_species',
                                           'gtdb_taxonomy',
                                           'ncbi_taxonomy'])

        metadata = {}
        with open(gtdb_metadata_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')

            genome_idx = header.index('accession')

            gtdb_type_designation_idx = header.index('gtdb_type_designation')
            gtdb_type_species_idx = header.index('gtdb_type_species_of_genus')

            gtdb_taxonomy_idx = header.index('gtdb_taxonomy')
            ncbi_taxonomy_idx = header.index('ncbi_taxonomy')

            gtdb_rep_idx = header.index('gtdb_representative')
            gtdb_genome_rep_idx = header.index('gtdb_genome_representative')

            for line in f:
                tokens = line.strip().split('\t')

                gid = canonical_gid(tokens[genome_idx])
                if gid.startswith('U'):
                    continue

                ncbi_taxonomy = [x.strip() for x in tokens[ncbi_taxonomy_idx].split(';')]
                gtdb_taxonomy = [x.strip() for x in tokens[gtdb_taxonomy_idx].split(';')]

                gtdb_is_rep = False
                if tokens[gtdb_rep_idx] == 't':
                    gtdb_is_rep = True

                gtdb_genome_rep = canonical_gid(tokens[gtdb_genome_rep_idx])
                if gtdb_genome_rep.lower() == 'none':
                    # genome failed QC and not assigned to GTDB species cluster
                    continue

                gtdb_type_strain = False
                if tokens[gtdb_type_designation_idx] == 'type strain of species':
                    gtdb_type_strain = True

                gtdb_type_species = tokens[gtdb_type_species_idx].lower() == 't'

                if gid in untrustworthy_type_gids:
                    # genome flagged as not being trustworthy as type material
                    # in the GTDB. Ideally, this would just be read as a field
                    # from the GTDB and GTDB type material fields would already
                    # reflect these genomes being marked as untrustworthy.
                    gtdb_type_strain = False
                    gtdb_type_species = False

                metadata[gid] = Metadata(gtdb_is_rep=gtdb_is_rep,
                                         gtdb_genome_rep=gtdb_genome_rep,
                                         gtdb_type_strain=gtdb_type_strain,
                                         gtdb_type_species=gtdb_type_species,
                                         gtdb_taxonomy=gtdb_taxonomy,
                                         ncbi_taxonomy=ncbi_taxonomy)

        return metadata

    def create_gtdb_sp_clusters(self, gtdb_metadata):
        """Create GTDB species clusters."""

        clusters = defaultdict(set)
        for gid, metadata in gtdb_metadata.items():
            rid = metadata.gtdb_genome_rep
            clusters[rid].add(gid)

        return clusters

    def identify_sp_synonyms(self, gtdb_metadata, gtdb_sp_clusters, lpsn, genus_synonyms):
        """Find NCBI species names considered synonyms under the GTDB.

        A NCBI species name is consider to be a synonym if:
         1) the type strain of a species is in a GTDB species cluster with a different name
         2) the species representative genome is a type strain and the GTDB and NCBI names differ
        """

        # create dictionary indicating all synonyms specified at LPSN
        lpsn_synonyms = defaultdict(set)
        for sp, synonyms in lpsn.sp_synonyms.items():
            for synonym in synonyms:
                lpsn_synonyms[sp].add(synonym)
                lpsn_synonyms[synonym].add(sp)

            for sp1, sp2 in itertools.permutations(synonyms, 2):
                lpsn_synonyms[sp1].add(sp2)

        # determine all species proposed in the GTDB
        gtdb_species = set()
        for metadata in gtdb_metadata.values():
            gtdb_sp = metadata.gtdb_taxonomy[Taxonomy.SPECIES_INDEX]
            gtdb_species.add(gtdb_sp)

        # determine synonyms
        fout = open(os.path.join(self.output_dir, 'species_synonyms.tsv'), 'w')
        fout.write('GTDB species\tGTDB taxonomy\tGTDB representative\tNCBI synonyms\tSynonyms not recognized at LPSN')
        fout.write('\tSynonym classification\tLPSN recognized name\tLPSN conflicting correct name\tLPSN synonyms\tType genomes\n')

        num_synonym_sp = 0
        for rid, cids in gtdb_sp_clusters.items():
            # check if there are type strain genomes in the cluster
            # with NCBI species assignments that differ from the GTDB
            # species assignment
            gtdb_sp = gtdb_metadata[rid].gtdb_taxonomy[Taxonomy.SPECIES_INDEX]
            if not is_latin_taxon(gtdb_sp):
                # placeholder names have no standing in nomenclature and thus
                # there is no recognized notion of names being synonyms
                continue

            type_gids = [gid for gid in cids if gtdb_metadata[gid].gtdb_type_strain]
            synonym_classification = set()
            ncbi_synonyms_gids = defaultdict(list)
            ncbi_synonyms = set()
            for type_gid in type_gids:
                ncbi_sp = gtdb_metadata[type_gid].ncbi_taxonomy[Taxonomy.SPECIES_INDEX]
                canonical_ncbi_sp = ncbi_sp.replace('Candidatus ', '')
                if canonical_ncbi_sp != gtdb_sp:
                    ncbi_synonyms_gids[ncbi_sp].append(type_gid)
                    ncbi_synonyms.add(ncbi_sp)

                    if canonical_ncbi_sp in gtdb_species:
                        self.logger.error('Invalid synonym: GTDB species = {}; NCBI species = {}; gid = {}'.format(
                            gtdb_sp,
                            ncbi_sp,
                            type_gid))

                    if rid == type_gid:
                        # need to use the generic name here since in some cases the NCBI genus
                        # assignment differs from the generic name
                        gtdb_generic = 'g__' + generic_name(gtdb_sp)
                        ncbi_generic = 'g__' + generic_name(canonical_ncbi_sp)
                        if gtdb_generic != ncbi_generic:
                            if ncbi_generic in genus_synonyms.get(gtdb_generic, set()):
                                synonym_classification.add('REP_MERGE_OF_GENERA')
                            else:
                                synonym_classification.add('REP_TRANSFER_OF_SPECIES')
                        else:
                            synonym_classification.add('REP_SPECIES_EPITHET_CHANGE?')
                    else:
                        synonym_classification.add('CLUSTERED_TYPE_STRAIN_GENOME')

            if len(ncbi_synonyms) > 0:
                num_synonym_sp += 1

                lpsn_correct_name = lpsn.sp_correct_names.get(gtdb_sp, '')
                if lpsn_correct_name == gtdb_sp:
                    lpsn_correct_name = ''

                ncbi_synonyms_str = []
                for ncbi_sp, sp_gids in ncbi_synonyms_gids.items():
                    ncbi_synonyms_str.append('{}: {}'.format(ncbi_sp, ', '.join(sp_gids)))

                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gtdb_sp,
                    '; '.join(gtdb_metadata[rid].gtdb_taxonomy[0:Taxonomy.SPECIES_INDEX]),
                    rid,
                    ', '.join(ncbi_synonyms),
                    ', '.join(ncbi_synonyms - lpsn_synonyms.get(gtdb_sp, set())),
                    ', '.join(synonym_classification),
                    gtdb_sp in lpsn.species,
                    lpsn_correct_name,
                    ', '.join(lpsn_synonyms.get(gtdb_sp, [])),
                    '; '.join(ncbi_synonyms_str)))

        fout.close()

        return num_synonym_sp

    def identify_genus_synonyms(self, gtdb_metadata, gtdb_sp_clusters, lpsn):
        """Find NCBI genera names considered synonyms under the GTDB.

        A NCBI genus name is consider to be a synonym if the type species of a genus is
        assigned to a different genus in the GTDB.
        """

        # create dictionary indicating all synonyms specified at LPSN
        lpsn_synonyms = defaultdict(set)
        for genus, synonyms in lpsn.genus_synonyms.items():
            for synonym in synonyms:
                lpsn_synonyms[genus].add(synonym)
                lpsn_synonyms[synonym].add(genus)

            for g1, g2 in itertools.permutations(synonyms, 2):
                lpsn_synonyms[g1].add(g2)

        # get genomes in each GTDB genus
        gtdb_genera_gids = defaultdict(set)
        for gid, metadata in gtdb_metadata.items():
            gtdb_genus = metadata.gtdb_taxonomy[Taxonomy.GENUS_INDEX]
            gtdb_genera_gids[gtdb_genus].add(gid)

        fout = open(os.path.join(self.output_dir, 'genus_synonyms.tsv'), 'w')
        fout.write('GTDB genus\tGTDB taxonomy\tNCBI synonyms\tSynonyms not recognized at LPSN')
        fout.write('\tLPSN recognized name\tLPSN conflicting correct name\tLPSN synonyms\tType genomes\n')

        synonyms = {}
        for gtdb_genus, gids in gtdb_genera_gids.items():
            # check if there are type species genomes in the GTDB genus
            # with NCBI genus assignments that differ from the GTDB
            # genus assignment
            if not is_latin_taxon(gtdb_genus):
                # placeholder names have no standing in nomenclature and thus
                # there is no recognized notion of names being synonyms
                continue

            type_gids = [gid for gid in gids if gtdb_metadata[gid].gtdb_type_species]

            ncbi_synonyms_gids = defaultdict(list)
            ncbi_synonyms = set()
            for type_gid in type_gids:
                ncbi_genus = gtdb_metadata[type_gid].ncbi_taxonomy[Taxonomy.GENUS_INDEX]
                canonical_ncbi_genus = ncbi_genus.replace('Candidatus ', '')
                ncbi_synonyms_gids[ncbi_genus].append(type_gid)
                if canonical_ncbi_genus != gtdb_genus:
                    ncbi_synonyms.add(ncbi_genus)

                    if canonical_ncbi_genus in gtdb_genera_gids:
                        self.logger.error('Invalid synonym: GTDB genus = {}; NCBI genus = {}; gid = {}'.format(
                            gtdb_genus,
                            ncbi_genus,
                            type_gid))

            if len(ncbi_synonyms) > 0:
                synonyms[gtdb_genus] = ncbi_synonyms

                ncbi_synonyms_str = []
                for ncbi_genus, genus_gids in ncbi_synonyms_gids.items():
                    ncbi_synonyms_str.append('{}: {}'.format(ncbi_genus, ', '.join(genus_gids)))

                lpsn_correct_name = lpsn.genus_correct_names.get(gtdb_genus, '')
                if lpsn_correct_name == gtdb_genus:
                    lpsn_correct_name = ''

                fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gtdb_genus,
                    '; '.join(gtdb_metadata[type_gids[0]].gtdb_taxonomy[0:Taxonomy.GENUS_INDEX]),
                    ', '.join(ncbi_synonyms),
                    ', '.join(ncbi_synonyms - lpsn_synonyms.get(gtdb_genus, set())),
                    gtdb_genus in lpsn.genera,
                    lpsn_correct_name,
                    ', '.join(lpsn_synonyms.get(gtdb_genus, [])),
                    '; '.join(ncbi_synonyms_str)))

        fout.close()

        return synonyms

    def run(self, gtdb_metadata_file, lpsn_gss_metadata, untrustworthy_type_material):
        """Generate files indicating NCBI names considered synonyms under the GTDB taxonomy."""

        self.logger.info('Identifying recognized taxa and synonyms at LPSN.')
        lpsn = LPSN(None, lpsn_gss_metadata)

        self.logger.info('Identifying genomes considered to be untrustworthy as type in GTDB.')
        untrustworthy_type_gids = self.parse_untrustworthy_type_genomes(untrustworthy_type_material)
        self.logger.info(f' - identified {len(untrustworthy_type_gids):,} untrustworthy type genomes.')

        self.logger.info('Reading GTDB metadata.')
        gtdb_metadata = self.read_gtdb_metadata(gtdb_metadata_file, untrustworthy_type_gids)
        self.logger.info(f' - read metadata for {len(gtdb_metadata):,} genomes')

        self.logger.info('Determining GTDB species clusters.')
        gtdb_sp_clusters = self.create_gtdb_sp_clusters(gtdb_metadata)
        self.logger.info(f' - identified {len(gtdb_sp_clusters):,} species clusters')

        self.logger.info('Identifying genus synonyms.')
        genus_synonyms = self.identify_genus_synonyms(gtdb_metadata, gtdb_sp_clusters, lpsn)
        self.logger.info(f' - identified {len(genus_synonyms):,} GTDB genera with NCBI synonyms')

        self.logger.info('Identifying species synonyms.')
        num_synonyms = self.identify_sp_synonyms(gtdb_metadata, gtdb_sp_clusters, lpsn, genus_synonyms)
        self.logger.info(f' - identified {num_synonyms:,} GTDB species with NCBI synonyms')
