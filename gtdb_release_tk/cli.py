import argparse
from contextlib import contextmanager

from biolib.misc.custom_help_formatter import CustomHelpFormatter


@contextmanager
def subparser(parser, name, desc):
    yield parser.add_parser(name, conflict_handler='resolve', help=desc,
                            formatter_class=CustomHelpFormatter)


@contextmanager
def mutex_group(parser, required):
    group = parser.add_argument_group(f'mutually exclusive '
                                      f'{"required" if required else "optional"} '
                                      f'arguments')
    yield group.add_mutually_exclusive_group(required=required)


@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)


def __help(group):
    group.add_argument('-h', '--help', action="help", help="show help message")


def __debug(group):
    group.add_argument('--debug', action="store_true", default=False,
                       help='create intermediate files for debugging purposes')


def __gtdb_sp_clusters_file(group):
    group.add_argument('--gtdb_sp_clusters_file', default=None,
                       help='file indicating final GTDB species clusters')


def __user_gid_table(group):
    group.add_argument(f'--user_gid_table', default=None,
                       help='table indicate NCBI accession and UBA IDs for user genomes')


def __silent(group):
    group.add_argument('--silent', action="store_true", default=False, help='suppress output')


def __bac_tree(group):
    group.add_argument('--bac_tree', default=None,
                       help='bacterial tree of representative genomes in Newick format')


def __ar_tree(group):
    group.add_argument('--ar_tree', default=None,
                       help='archaeal tree of representative genomes in Newick format')


def __genome_path_file(group):
    group.add_argument('--genome_path_file', default=None, help='file indicating path to GTDB genomes')


def __min_16s_ar_len(group):
    group.add_argument('--min_16S_ar_len', default=0, type=int,
                       help='minimum length of archaeal 16S rRNA gene to include')


def __min_16s_bac_len(group):
    group.add_argument('--min_16S_bac_len', default=0, type=int,
                       help='minimum length of bacterial 16S rRNA gene to include')


def __min_contig_len(group):
    group.add_argument('--min_contig_len', default=0, type=int,
                       help='minimum contig length')


def __bac120_gene_dir(group):
    group.add_argument('bac120_gene_dir', help='directory with individual bacterial marker gene files')


def __ar122_gene_dir(group):
    group.add_argument('ar122_gene_dir', help='directory with individual archaeal marker gene files')


def __taxonomy_file(group, required=False):
    group.add_argument(f'--taxonomy_file', required=required, help='tsv of gid and taxonomy')


def __genome_dirs(group):
    group.add_argument('--genome_dirs', help='File listing the path to all NCBI genome files.')


def __uba_path(group):
    group.add_argument('--uba_path', help='root path for files not present in NCBI (UBA).')


def __cpus(group):
    group.add_argument('--cpus', type=int, default=1, help='number of threads to use (default: 1)')


def __bac120_msa_file(group):
    group.add_argument('--bac120_msa_file', help='MSA for bacterial representative genomes')


def __ar122_msa_file(group):
    group.add_argument('--ar122_msa_file', help='MSA for archaeal representative genomes')


def __metadata_fields(group):
    group.add_argument('--metadata_fields', help='file indicating GTDB metadata fields to retain')


def __tree_file(group):
    group.add_argument('--tree_file', help='tree spanning representative genomes')


def __metadata_file(group, required=False):
    group.add_argument(f'--metadata_file', required=required,
                       help='file indicating metadata of each genome')


def __msa_file(group):
    group.add_argument('--msa_file', help='file with MSA of representative genomes')


def __ssu_file(group):
    group.add_argument('--ssu_file', help='file with SSU genes for representative genomes')


def __sp_clusters_file(group):
    group.add_argument('--sp_clusters_file', help='file indicating GTDB species clusters')


def __hq_genome_file(group):
    group.add_argument('--hq_genome_file', help='file indicating HQ MIMAG genomes')


def __output_tree(group):
    group.add_argument('--output_tree', help='output tree with appended species labels')


def __reference_taxonomy_file(group):
    group.add_argument('--reference_taxonomy_file', help='Old taxonomy file used as reference.')


def __new_taxonomy_file(group):
    group.add_argument('--new_taxonomy_file', help='New taxonomy file.')


def __output_dir(group, required=False):
    group.add_argument(f'--output_dir', required=required,
                       help='output directory')


def __changes_only(group):
    group.add_argument('--changes_only', action="store_true", default=False,
                       help='only report that differ between the two taxonomies')


def __filter_placeholder_name(group):
    group.add_argument('--filter_placeholder_name', action="store_true", default=False,
                       help='remove taxa that have a placeholder name -> no name mapping between the two taxonomies')


def __top_change(group):
    group.add_argument('--filter_placeholder_name', action="store_true", default=False,
                       help='remove taxa that have a placeholder name -> no name mapping between the two taxonomies')


def __release_number(group, required=False):
    group.add_argument(f'--release_number', required=required, help='Release number (ex. 86)')


def __bac120_metadata_file(group):
    group.add_argument('--bac120_metadata_file', help='metadata for bacterial GTDB genomes.')


def __ar122_metadata_file(group):
    group.add_argument('--ar122_metadata_file', help='metadata for archaeal GTDB genomes')


def __all_genomes(group):
    group.add_argument('--all_genomes', action="store_true", default=False,
                       help='create plot over all genomes instead of all GTDB species representatives')


def __width(group):
    group.add_argument('--width', type=float, default=6.5, help='width of plot (inches)')


def __height(group):
    group.add_argument('--height', type=float, default=4, help='width of plot (inches)')


def __domain(group):
    group.add_argument('--domain', choices=['Both', 'Archaea', 'Bacteria'], default='Both',
                       help='restrict plot to genomes from specified domain')


def __bac120_tree(group):
    group.add_argument('--bac120_tree', help='bacterial GTDB reference tree')


def __ar122_tree(group):
    group.add_argument('--ar122_tree', help='archaeal GTDB reference tree')


def __genomes_per_taxon(group):
    group.add_argument('--genomes_per_taxon', type=int, default=1, help='genomes to sample from each taxon')


def __min_ssu_len(group):
    group.add_argument('--min_ssu_len', type=int, default=900, help='minimum length of 16S rRNA gene')


def __min_msa_perc(group):
    group.add_argument('--min_ssu_len', type=float, default=30, help='minimum percentage of amino acids in MSA')


def __ssu_fasta_file(group):
    group.add_argument('--ssu_fasta_file', help='FASTA file with SSU sequences for all GTDB genomes')


def __num_taxa(group):
    group.add_argument('--num_taxa', type=int, default=10, help='number of taxa to report')


def __ncbi_node_file(group, required=False):
    group.add_argument('--ncbi_node_file', help='NCBI node.dmp file.', required=required)


def __ncbi_name_file(group, required=False):
    group.add_argument('--ncbi_name_file', help='NCBI name.dmp file.', required=required)


def __lpsn_species_file(group, required=False):
    group.add_argument('--lpsn_species_file', help='LPSN species file.', required=required)


def __bacdive_species_file(group, required=False):
    group.add_argument('--bacdive_species_file', help='DSMZ species file.', required=required)


def __gtdb_taxonomy(group, required=False):
    group.add_argument('--gtdb_taxonomy', help='GTDB Taxonomy file.', required=required)


def get_main_parser():
    # Setup the main, and sub parsers.
    main_parser = argparse.ArgumentParser(prog='gtdbtk', add_help=False, conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    # External files for GTDB website

    # taxonomy_files  -> Generate GTDB taxonomy files
    with subparser(sub_parsers, 'taxonomy_files', 'Generate taxonomy files for GTDB website.') as parser:
        with arg_group(parser, 'positional arguments') as grp:
            __metadata_file(grp)
            __gtdb_sp_clusters_file(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    # tree_files -> Generate tree files spanning representative genomes
    with subparser(sub_parsers, 'tree_files', 'Generate tree files spanning representative genomes.') as parser:
        with arg_group(parser, 'positional arguments') as grp:
            __metadata_file(grp)
            __bac_tree(grp)
            __ar_tree(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   sp_cluster_file -> Generate file indicating GTDB species clusters
    with subparser(sub_parsers, 'sp_cluster_file', 'Generate file indicating GTDB species clusters.') as parser:
        with arg_group(parser, 'positional arguments') as grp:
            __metadata_file(grp)
            __gtdb_sp_clusters_file(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   hq_genome_file  -> Generate file indicating HQ MIMAG genomes
    with subparser(sub_parsers, 'hq_genome_file', 'Generate file indicating HQ genomes.') as parser:
        with arg_group(parser, 'positional arguments') as grp:
            __metadata_file(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   ssu_files       -> Generate 16S rRNA FASTA files
    with subparser(sub_parsers, 'ssu_files', 'Generate 16S rRNA FASTA files for GTDB website.') as parser:
        with arg_group(parser, 'positional arguments') as grp:
            __metadata_file(grp)
            __gtdb_sp_clusters_file(grp)
            __user_gid_table(grp)
            __genome_path_file(grp)
            __release_number(grp)
            __output_dir(grp)
            __min_16s_ar_len(grp)
            __min_16s_bac_len(grp)
            __min_contig_len(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   marker_files    -> Generate marker gene files
    with subparser(sub_parsers, 'marker_files', 'Generate marker gene file.') as parser:
        with arg_group(parser, 'positional arguments') as grp:
            __bac120_gene_dir(grp)
            __ar122_gene_dir(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   protein_files   -> Generate the archive containing all protein (AA) for representative genomes.
    with subparser(sub_parsers, 'protein_files',
                   'Generate the archive containing all protein (AA) for representative genomes.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __taxonomy_file(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __genome_dirs(grp)
            __uba_path(grp)
            __output_dir(grp)
            __cpus(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   nucleotide_files-> Generate the archive containing all protein (nt) for representative genomes.
    with subparser(sub_parsers, 'nucleotide_files',
                   'Generate the archive containing all protein (nucl) for representative genomes.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __taxonomy_file(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __genome_dirs(grp)
            __uba_path(grp)
            __output_dir(grp)
            __cpus(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   msa_files       -> Generate concatenated MSA files
    with subparser(sub_parsers, 'msa_files', 'Generate concatenated MSA files.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_msa_file(grp)
            __ar122_msa_file(grp)
            __metadata_file(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   gene_files      -> Generate untrimmed ar122/bac120 gene files.
    with subparser(sub_parsers, 'gene_files', 'Generate concatenated untrimmed gene files.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __taxonomy_file(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __output_dir(grp)
            __cpus(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   metadata_files  -> Generate GTDB metadata files
    with subparser(sub_parsers, 'metadata_files', 'Generate GTDB metadata files.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __metadata_file(grp)
            __metadata_fields(grp)
            __user_gid_table(grp)
            __gtdb_sp_clusters_file(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   dict_file       -> Generate dictionary containing all GTDB taxa
    with subparser(sub_parsers, 'dict_file', 'Create dictionary with all GTDB taxa.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __taxonomy_file(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   arb_files       -> Generate ARB metadata file
    with subparser(sub_parsers, 'arb_files', 'Generate ARB metadata file.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __metadata_file(grp)
            __metadata_fields(grp)
            __user_gid_table(grp)
            __bac120_msa_file(grp)
            __ar122_msa_file(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   validate        -> Perform validation checks on GTDB website files
    with subparser(sub_parsers, 'validate', 'Perform validation checks on GTDB website files') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __taxonomy_file(grp)
            __tree_file(grp)
            __metadata_file(grp)
            __msa_file(grp)
            __ssu_file(grp)
            __sp_clusters_file(grp)
            __hq_genome_file(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   add_sp_label    -> Generate tree with species labels
    with subparser(sub_parsers, 'add_sp_label', 'Generate tree with species labels.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __tree_file(grp)
            __taxonomy_file(grp)
            __output_tree(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    # Internal files for GTDB website:
    #   tax_comp_files  -> Generate comparison files between 2 taxonomy ordered by degree of polyphyly
    with subparser(sub_parsers, 'tax_comp_files',
                   'Generate comparison files between 2 taxonomy ordered by degree of polyphyly.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __reference_taxonomy_file(grp)
            __new_taxonomy_file(grp)
            __output_dir(grp)
            __user_gid_table(grp)
            __changes_only(grp)
            __filter_placeholder_name(grp)
            __top_change(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   json_tree_file  -> Generate JSON file used to populate http://gtdb.ecogenomic.org/tree
    with subparser(sub_parsers, 'json_tree_file',
                   'Generate JSON file used to create the tree in http://gtdb.ecogenomic.org/tree.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __taxonomy_file(grp, required=True)
            __metadata_file(grp, required=True)
            __output_dir(grp, required=True)
            __release_number(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    # Auxiliary data files:
    #   reps_per_rank -> Select representative genomes at each taxonomic rank
    with subparser(sub_parsers, 'reps_per_rank', 'Select representative genomes at each taxonomic rank.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __bac120_msa_file(grp)
            __ar122_msa_file(grp)
            __ssu_fasta_file(grp)
            __user_gid_table(grp)
            __release_number(grp)
            __output_dir(grp)
            __genomes_per_taxon(grp)
            __min_ssu_len(grp)
            __min_msa_perc(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)
    #
    # Plots:
    #   genomic_stats         -> Histograms of common genomic statistics (e.g., genome size, coding density)
    with subparser(sub_parsers, 'genomic_stats', 'Plot of common genomic statistics.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __release_number(grp)
            __output_dir(grp)
            __all_genomes(grp)
            __width(grp)
            __height(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   genome_quality        -> Scatter plot showing quality of GTDB representative genomes
    with subparser(sub_parsers, 'genome_quality',
                   'Scatter plot showing quality of GTDB representative genomes.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   genome_category_rank  -> Bar plot of MAGs, SAGs, and isolates within each taxonomic rank
    with subparser(sub_parsers, 'genome_category_rank',
                   'Plot number of MAGs, SAGs, and isolates for each taxonomic rank.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   nomenclatural_rank    -> Bar plot of nomenclatural status of taxa within each taxonomic rank
    with subparser(sub_parsers, 'nomenclatural_rank',
                   'Plot nomenclatural status of species for each taxonomic rank.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __release_number(grp)
            __output_dir(grp)
            __domain(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   ncbi_compare          -> Bar plot comparing GTDB and NCBI taxonomies
    with subparser(sub_parsers, 'ncbi_compare', 'Bar plot comparing GTDB and NCBI taxonomies.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __release_number(grp)
            __output_dir(grp)
            __all_genomes(grp)
            __domain(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   sp_rep_type           -> Pie chart indicating type information for the GTBD species representatives
    with subparser(sub_parsers, 'sp_rep_type',
                   'Pie chart indicating type information for the GTBD species representatives.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __release_number(grp)
            __output_dir(grp)
            __domain(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   itol                  -> Generate tree and iTOL files for producing iTOL tree image
    with subparser(sub_parsers, 'itol', 'Generate tree and iTOL files for producing iTOL tree image.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_tree(grp)
            __ar122_tree(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #
    # Tables:
    #   taxa_count -> Create table with number of taxa at each taxonomic rank
    with subparser(sub_parsers, 'taxa_count', 'Create table with number of taxa at each taxonomic rank.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __bac120_metadata_file(grp)
            __ar122_metadata_file(grp)
            __release_number(grp)
            __output_dir(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    #   top_taxa   -> Create table with most abundant taxa at each taxonomic rank
    with subparser(sub_parsers, 'top_taxa', 'Create table with most abundant taxa at each taxonomic rank.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __gtdb_sp_clusters_file(grp)
            __release_number(grp)
            __output_dir(grp)
            __num_taxa(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    # Meta-database:
    #   nomenclatural_check  -> List of latin names present/absent from NCBI,LPSN,Bacdive.
    with subparser(sub_parsers, 'nomenclatural_check',
                   'List of latin names present/absent from NCBI,LPSN,Bacdive.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __ncbi_node_file(grp, required=True)
            __ncbi_name_file(grp, required=True)
            __lpsn_species_file(grp, required=True)
            __bacdive_species_file(grp, required=True)
            __output_dir(grp, required=True)
            __gtdb_taxonomy(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __silent(grp)
            __debug(grp)
            __help(grp)

    return main_parser
