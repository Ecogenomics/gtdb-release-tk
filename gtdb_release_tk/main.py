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

import os
import sys
import csv
import re
import logging
from collections import defaultdict

import dendropy

from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists
from biolib.taxonomy import Taxonomy

from gtdb_release_tk.website_data import WebsiteData
from gtdb_release_tk.reps_per_rank import RepsPerRank
from gtdb_release_tk.itol import iTOL
from gtdb_release_tk.lpsn import LPSN
from gtdb_release_tk.bacdive import BacDive
from gtdb_release_tk.strains import Strains


from gtdb_release_tk.plots.genome_category_per_rank import GenomeCategoryPerRank
from gtdb_release_tk.plots.nomenclatural_per_rank import NomenclaturalPerRank
from gtdb_release_tk.plots.genomic_stats import GenomicStats
from gtdb_release_tk.plots.genome_quality import GenomeQuality
from gtdb_release_tk.plots.ncbi_compare import NCBI_Compare
from gtdb_release_tk.plots.species_rep_type import SpeciesRepType

from gtdb_release_tk.tables.taxa_count import TaxaCount
from gtdb_release_tk.tables.top_taxa import TopTaxa


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def taxonomy_files(self, options):
        """Generate taxonomy files for GTDB website."""

        check_file_exists(options.metadata_file)
        check_file_exists(options.gtdb_sp_clusters_file)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.taxonomy_files(options.metadata_file,
                         options.gtdb_sp_clusters_file,
                         options.user_gid_table)

        self.logger.info('Done.')

    def tree_files(self, options):
        """Generate tree files spanning representative genomes."""

        check_file_exists(options.metadata_file)
        check_file_exists(options.bac_tree)
        check_file_exists(options.ar_tree)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.tree_files(options.metadata_file,
                     options.bac_tree,
                     options.ar_tree,
                     options.user_gid_table)

        self.logger.info('Done.')

    def sp_cluster_file(self, options):
        """Generate file indicating GTDB species clusters."""

        check_file_exists(options.metadata_file)
        check_file_exists(options.gtdb_sp_clusters_file)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.sp_cluster_file(options.metadata_file,
                          options.gtdb_sp_clusters_file,
                          options.user_gid_table)

        self.logger.info('Done.')

    def hq_genome_file(self, options):
        """Generate file indicating HQ genomes."""

        check_file_exists(options.metadata_file)
        check_file_exists(options.gtdb_sp_clusters_file)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.hq_genome_file(options.metadata_file,
                         options.gtdb_sp_clusters_file,
                         options.user_gid_table)

        self.logger.info('Done.')

    def ssu_files(self, options):
        """Generate 16S rRNA FASTA files for GTDB website."""

        check_file_exists(options.metadata_file)
        check_file_exists(options.gtdb_sp_clusters_file)
        check_file_exists(options.user_gid_table)
        check_file_exists(options.genome_path_file)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.ssu_files(options.metadata_file,
                    options.gtdb_sp_clusters_file,
                    options.user_gid_table,
                    options.genome_path_file,
                    options.min_16S_ar_len,
                    options.min_16S_bac_len,
                    options.min_contig_len)

        self.logger.info('Done.')

    def marker_files(self, options):
        """Generate marker gene file."""

        check_dir_exists(options.bac120_gene_dir)
        check_dir_exists(options.ar122_gene_dir)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.marker_files(options.bac120_gene_dir,
                       options.ar122_gene_dir,
                       options.user_gid_table)

        self.logger.info('Done.')

    def msa_files(self, options):
        """Generate concatenated MSA files."""

        check_file_exists(options.bac120_msa_file)
        check_file_exists(options.ar122_msa_file)
        check_file_exists(options.metadata_file)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.msa_files(options.bac120_msa_file,
                    options.ar122_msa_file,
                    options.metadata_file,
                    options.user_gid_table)

        self.logger.info('Done.')

    def metadata_files(self, options):
        """Generate GTDB metadata files."""

        check_file_exists(options.metadata_file)
        check_file_exists(options.metadata_fields)
        check_file_exists(options.gtdb_sp_clusters_file)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.metadata_files(options.metadata_file,
                         options.metadata_fields,
                         options.gtdb_sp_clusters_file,
                         options.user_gid_table)

        self.logger.info('Done.')

    def dict_file(self, options):
        """Generate GTDB dictionary file."""

        check_file_exists(options.taxonomy_file)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.dict_file(options.taxonomy_file)

        self.logger.info('Done.')

    def arb_files(self, options):
        """Generate ARB metadata file."""

        self.logger.warning('IN PROGRESS')
        return

        check_file_exists(options.metadata_file)
        check_file_exists(options.metadata_fields)
        check_file_exists(options.user_gid_table)
        check_file_exists(options.bac120_msa_file)
        check_file_exists(options.ar122_msa_file)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.arb_files(options.metadata_file,
                    options.metadata_fields,
                    options.user_gid_table,
                    options.bac120_msa_file,
                    options.ar122_msa_file)

        self.logger.info('Done.')

    def validate(self, options):
        """Perform validation checks on GTDB website files."""

        check_file_exists(options.taxonomy_file)
        check_file_exists(options.tree_file)
        check_file_exists(options.metadata_file)
        check_file_exists(options.msa_file)
        check_file_exists(options.ssu_file)
        check_file_exists(options.sp_clusters_file)
        check_file_exists(options.hq_genome_file)

        p = WebsiteData(None, None)
        p.validate(options.taxonomy_file,
                   options.tree_file,
                   options.metadata_file,
                   options.msa_file,
                   options.ssu_file,
                   options.sp_clusters_file,
                   options.hq_genome_file)

        self.logger.info('Done.')
        
    def add_sp_label(self, options):
        """Generate tree with species labels."""

        check_file_exists(options.taxonomy_file)
        check_file_exists(options.tree_file)
        
        self.logger.info('Reading GTDB taxonomy.')
        gtdb_taxonomy = Taxonomy().read(options.taxonomy_file)
        
        self.logger.info('Reading input tree.')
        tree = dendropy.Tree.get_from_path(options.tree_file, 
                                            schema='newick', 
                                            rooting='force-rooted', 
                                            preserve_underscores=True)
                                            
        self.logger.info('Appending species labels.')
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                gid = node.taxon.label
                species = gtdb_taxonomy[gid][Taxonomy.SPECIES_INDEX].replace('s__', '')
                node.taxon.label += ' | {}'.format(species)
                
        self.logger.info('Writing output tree.')
        tree.write_to_path(options.output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)

        self.logger.info('Done.')

    def tax_comp_files(self, options):
        """Generate comparison files between 2 taxonomy ordered by degree of polyphyly."""
        check_file_exists(options.reference_taxonomy_file)
        check_file_exists(options.new_taxonomy_file)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData('', options.output_dir)
        p.tax_comp_files(options.reference_taxonomy_file,
                         options.new_taxonomy_file,
                         options.user_gid_table,
                         options.changes_only,
                         options.filter_placeholder_name,
                         options.top_change)

        self.logger.info('Done.')

    def json_tree_parser(self, options):
        """generate Json file used to create the tree in http://gtdb.ecogenomic.org/tree"""

        p = WebsiteData(options.release_number, options.output_dir)
        p.json_tree_parser(options.taxonomy_file,
                           options.metadata_file)
        self.logger.info('Done.')

    def genome_category_rank(self, options):
        """Plot number of MAGs, SAGs, and isolates for each taxonomic rank."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        make_sure_path_exists(options.output_dir)

        p = GenomeCategoryPerRank(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file)

        self.logger.info('Done.')

    def nomenclatural_rank(self, options):
        """Plot nomenclatural status of species for each taxonomic rank."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        make_sure_path_exists(options.output_dir)

        p = NomenclaturalPerRank(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file,
              options.domain)

        self.logger.info('Done.')

    def ncbi_compare(self, options):
        """Bar plot comparing GTDB and NCBI taxonomies."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        make_sure_path_exists(options.output_dir)

        p = NCBI_Compare(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file,
              options.all_genomes,
              options.domain)

        self.logger.info('Done.')

    def sp_rep_type(self, options):
        """Pie chart indicating type information for the GTBD species representatives."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        make_sure_path_exists(options.output_dir)

        p = SpeciesRepType(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file,
              options.domain)

        self.logger.info('Done.')

    def itol(self, options):
        """Generate tree and iTOL files for producing iTOL tree image."""

        check_file_exists(options.bac120_tree)
        check_file_exists(options.ar122_tree)
        make_sure_path_exists(options.output_dir)

        p = iTOL(options.release_number, options.output_dir)
        p.run(options.bac120_tree,
              options.ar122_tree)

        self.logger.info('Done.')

    def genomic_stats(self, options):
        """Plot of common genomic statistics."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        make_sure_path_exists(options.output_dir)

        p = GenomicStats(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file,
              options.all_genomes,
              options.width,
              options.height)

        self.logger.info('Done.')

    def genome_quality(self, options):
        """Scatter plot showing quality of GTDB representative genomes."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        make_sure_path_exists(options.output_dir)

        p = GenomeQuality(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file)

        self.logger.info('Done.')

    def reps_per_rank(self, options):
        """Select representative genomes at each taxonomic rank."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        check_file_exists(options.bac120_msa_file)
        check_file_exists(options.ar122_msa_file)
        check_file_exists(options.ssu_fasta_file)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = RepsPerRank(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file,
              options.bac120_msa_file,
              options.ar122_msa_file,
              options.ssu_fasta_file,
              options.user_gid_table,
              options.genomes_per_taxon,
              options.min_ssu_len,
              options.min_msa_perc)

        self.logger.info('Done.')

    def taxa_count(self, options):
        """Create table with number of taxa at each taxonomic rank."""

        check_file_exists(options.bac120_metadata_file)
        check_file_exists(options.ar122_metadata_file)
        make_sure_path_exists(options.output_dir)

        p = TaxaCount(options.release_number, options.output_dir)
        p.run(options.bac120_metadata_file,
              options.ar122_metadata_file)

        self.logger.info('Done.')

    def top_taxa(self, options):
        """Create table with most abundant taxa at each taxonomic rank."""

        check_file_exists(options.gtdb_sp_clusters_file)
        make_sure_path_exists(options.output_dir)

        p = TopTaxa(options.release_number, options.output_dir)
        p.run(options.gtdb_sp_clusters_file,
              options.num_taxa)

        self.logger.info('Done.')

    def full_lpsn_wf(self, options):
        """Full workflow to parse LPSN."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.full_lpsn_wf()

    def pull_html(self, options):
        """Pull all genus.html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.pull_html()

    def parse_html(self, options):
        """Parse all html files."""
        make_sure_path_exists(options.output_dir)
        p = LPSN(options.output_dir)
        p.parse_html(options.input_dir)

    def download_strains(self, options):
        make_sure_path_exists(options.output_dir)
        p = BacDive(options.output_dir, options.username, options.pwd)
        p.download_strains()

    def generate_date_table(self, options):
        p = Strains()
        p.generate_date_table(options.lpsn_species_info,
                              options.dsmz_species_info,
                              options.straininfo_species_info,
                              options.output_file)

    def generate_type_table(self, options):
        p = Strains(options.output_dir, options.cpus)
        p.generate_type_strain_table(options.metadata_file,
                                     options.ncbi_names,
                                     options.ncbi_nodes,
                                     options.lpsn_dir,
                                     options.dsmz_dir,
                                     options.straininfo_dir,
                                     options.year_table,
                                     options.source_strain)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if options.subparser_name == 'taxonomy_files':
            self.taxonomy_files(options)
        elif options.subparser_name == 'tree_files':
            self.tree_files(options)
        elif options.subparser_name == 'sp_cluster_file':
            self.sp_cluster_file(options)
        elif options.subparser_name == 'hq_genome_file':
            self.hq_genome_file(options)
        elif options.subparser_name == 'ssu_files':
            self.ssu_files(options)
        elif options.subparser_name == 'marker_files':
            self.marker_files(options)
        elif options.subparser_name == 'msa_files':
            self.msa_files(options)
        elif options.subparser_name == 'metadata_files':
            self.metadata_files(options)
        elif options.subparser_name == 'dict_file':
            self.dict_file(options)
        elif options.subparser_name == 'arb_files':
            self.arb_files(options)
        elif options.subparser_name == 'validate':
            self.validate(options)
        elif options.subparser_name == 'add_sp_label':
            self.add_sp_label(options)
        elif options.subparser_name == 'genome_category_rank':
            self.genome_category_rank(options)
        elif options.subparser_name == 'nomenclatural_rank':
            self.nomenclatural_rank(options)
        elif options.subparser_name == 'ncbi_compare':
            self.ncbi_compare(options)
        elif options.subparser_name == 'sp_rep_type':
            self.sp_rep_type(options)
        elif options.subparser_name == 'itol':
            self.itol(options)
        elif options.subparser_name == 'tax_comp_files':
            self.tax_comp_files(options)
        elif options.subparser_name == 'json_tree_file':
            self.json_tree_parser(options)
        elif options.subparser_name == 'genomic_stats':
            self.genomic_stats(options)
        elif options.subparser_name == 'genome_quality':
            self.genome_quality(options)
        elif options.subparser_name == 'reps_per_rank':
            self.reps_per_rank(options)
        elif options.subparser_name == 'taxa_count':
            self.taxa_count(options)
        elif options.subparser_name == 'top_taxa':
            self.top_taxa(options)
        elif options.subparser_name == 'lpsn':
            if options.lpsn_subparser_name == 'lpsn_wf':
                self.full_lpsn_wf(options)
            elif options.lpsn_subparser_name == 'parse_html':
                self.parse_html(options)
            elif options.lpsn_subparser_name == 'pull_html':
                self.pull_html(options)
            else:
                self.logger.error('Unknown command: ' +
                                  options.lpsn_subparser_name + '\n')
        elif options.subparser_name == 'bacdive':
            if options.bacdive_subparser_name == 'download_strains':
                self.download_strains(options)
        elif options.subparser_name == 'strains':
            if options.strains_subparser_name == 'date_table':
                self.generate_date_table(options)
            if options.strains_subparser_name == 'type_table':
                self.generate_type_table(options)

        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0
