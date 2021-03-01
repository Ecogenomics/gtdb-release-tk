###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #r
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2019"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import sys

from biolib.logger import logger_setup

from gtdb_release_tk import version
from gtdb_release_tk.cli import get_main_parser
from gtdb_release_tk.main import OptionsParser


def print_help():
    """Help menu."""

    print('')
    print('                ...::: GTDB Release Toolkit v' + version() + ' :::...''')
    print('''\

    External files for GTDB website:
      taxonomy_files  -> Generate GTDB taxonomy files
      tree_files      -> Generate tree files spanning representative genomes
      sp_cluster_file -> Generate file indicating GTDB species clusters
      hq_genome_file  -> Generate file indicating HQ MIMAG genomes
      ssu_files       -> Generate 16S rRNA FASTA files
      marker_files    -> Generate marker gene files
      protein_files   -> Generate the archive containing all protein (AA) for representative genomes.
      nucleotide_files-> Generate the archive containing all protein (nt) for representative genomes.
      msa_files       -> Generate concatenated MSA files
      gene_files      -> Generate untrimmed ar122/bac120 gene files.
      metadata_files  -> Generate GTDB metadata files
      dict_file       -> Generate dictionary containing all GTDB taxa
      arb_files       -> Generate ARB metadata file
      validate        -> Perform validation checks on GTDB website files
      add_sp_label    -> Generate tree with species labels

    Internal files for GTDB website:
      tax_comp_files  -> Generate comparison files between 2 taxonomy ordered by degree of polyphyly
      json_tree_file  -> Generate JSON file used to populate http://gtdb.ecogenomic.org/tree

    Auxiliary data files:
      reps_per_rank -> Select representative genomes at each taxonomic rank

    Plots:
      genomic_stats         -> Histograms of common genomic statistics (e.g., genome size, coding density)
      genome_quality        -> Scatter plot showing quality of GTDB representative genomes
      genome_category_rank  -> Bar plot of MAGs, SAGs, and isolates within each taxonomic rank
      nomenclatural_rank    -> Bar plot of nomenclatural status of taxa within each taxonomic rank
      ncbi_compare          -> Bar plot comparing GTDB and NCBI taxonomies
      sp_rep_type           -> Pie chart indicating type information for the GTBD species representatives
      itol                  -> Generate tree and iTOL files for producing iTOL tree image

    Tables:
      taxa_count -> Create table with number of taxa at each taxonomic rank
      top_taxa   -> Create table with most abundant taxa at each taxonomic rank

    Meta-database:
      nomenclatural_check  -> List of latin names present/absent from NCBI,LPSN,Bacdive.


  Use: gtdb_release_tk <command> -h for command specific help.

  Feature requests or bug reports can be sent to Donovan Parks (donovan.parks@gmail.com)
    or posted on GitHub (https://github.com/Ecogenomics/gtdb_release_tk).
    ''')


def main():
    # get and check options
    args = None
    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv == '--help'):
        print_help()
        sys.exit(0)
    else:
        args = get_main_parser().parse_args()

    try:
        logger_setup(args.output_dir,
                     'gtdb_release_tk.log',
                     'GTDB Release Tk',
                     version(),
                     args.silent)
    except:
        logger_setup(None,
                     'gtdb_release_tk.log',
                     'GTDB Release Tk',
                     version(),
                     args.silent)

    # do what we came here to do
    try:
        parser = OptionsParser()
        parser.parse_options(args)
    except SystemExit:
        print("\n  Controlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise


if __name__ == '__main__':
    main()
