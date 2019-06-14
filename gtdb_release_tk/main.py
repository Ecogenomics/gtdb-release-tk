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

from biolib.common import check_file_exists, make_sure_path_exists

from gtdb_release_tk.website_data import WebsiteData
from gtdb_release_tk.plots.genome_types_per_rank import GenomeTypesPerRank


class OptionsParser():
    def __init__(self):
        """Initialization"""

        self.logger = logging.getLogger('timestamp')

    def taxonomy_files(self, options):
        """Generate taxonomy files for GTDB website."""

        check_file_exists(options.metadata_file)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.taxonomy_files(options.metadata_file,
                         options.user_gid_table)

        self.logger.info('Done.')

    def tree_files(self, options):
        """Generate tree files spanning representative genomes."""

        check_file_exists(options.bac_tree)
        check_file_exists(options.ar_tree)
        check_file_exists(options.user_gid_table)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData(options.release_number, options.output_dir)
        p.tree_files(options.bac_tree,
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

    def tax_comp_files(self, options):
        """Generate comparison files between 2 taxonomy ordered by degree of polyphyly."""
        check_file_exists(options.reference_taxonomy_file)
        check_file_exists(options.new_taxonomy_file)
        make_sure_path_exists(options.output_dir)

        p = WebsiteData('', options.output_dir)
        p.tax_comp_files(options.reference_taxonomy_file,
                         options.new_taxonomy_file,
                         options.user_gid_table,
                         options.keep_no_change,
                         options.keep_placeholder_name,
                         options.top_change)

        self.logger.info('Done.')

    def json_tree_parser(self, options):
        """generate Json file used to create the tree in http://gtdb.ecogenomic.org/tree"""
        print('in DEV')
        self.logger.info('Done.')

        p = WebsiteData(options.release_number, options.output_dir)
        p.json_tree_parser(options.taxonomy_file,
                           options.metadata_file,
                           options.user_gid_table)

    def genome_types_per_rank(self, options):
        """Plot number of MAGs, SAGs, and isolates for each taxonomic rank."""

        check_file_exists(options.metadata_file)
        make_sure_path_exists(options.output_dir)

        p = GenomeTypesPerRank(options.output_dir)
        p.run(options.metadata_file)

        self.logger.info('Done.')

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if options.subparser_name == 'taxonomy_files':
            self.taxonomy_files(options)
        elif options.subparser_name == 'tree_files':
            self.tree_files(options)
        elif options.subparser_name == 'sp_cluster_file':
            self.sp_cluster_file(options)
        elif options.subparser_name == 'ssu_files':
            self.ssu_files(options)
        elif options.subparser_name == 'genome_types_per_rank':
            self.genome_types_per_rank(options)
        elif options.subparser_name == 'tax_comp_files':
            self.tax_comp_files(options)
        elif options.subparser_name == 'json_tree_parser':
            self.json_tree_parser(options)
        else:
            self.logger.error('Unknown command: ' +
                              options.subparser_name + '\n')
            sys.exit()

        return 0
