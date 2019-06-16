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
import shutil
import logging
from pathlib import Path, PurePath
from collections import defaultdict
from collections import Counter


import dendropy

from numpy import (arange as np_arange,
                   array as np_array)

from biolib.taxonomy import Taxonomy
import biolib.seq_io as seq_io

from gtdb_release_tk.common import (parse_user_gid_table,
                                    parse_rep_genomes,
                                    parse_genomic_path_file,
                                    parse_species_clusters)

from gtdb_release_tk.config import SSU_DIR, MAX_SSU_LEN


class WebsiteData(object):
    """Generate data files for website."""

    def __init__(self, release_number, output_dir):
        """Initialization."""

        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def taxonomy_files(self, metadata_file, user_gid_table):
        """Generate taxonomy files for GTDB website."""

        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file, user_gids)
        self.logger.info(' ...identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))

        # create archaeal and bacterial taxonomy files
        fout_ar = open(
            self.output_dir / 'ar122_taxonomy_r{}.tsv'.format(self.release_number), 'w')
        fout_bac = open(
            self.output_dir / 'bac120_taxonomy_r{}.tsv'.format(self.release_number), 'w')

        for gid, taxa in reps.items():
            assert(not gid.startswith('U_'))

            if taxa[0] == 'd__Archaea':
                fout_ar.write('{}\t{}\n'.format(gid, ';'.join(taxa)))
            elif taxa[0] == 'd__Bacteria':
                fout_bac.write('{}\t{}\n'.format(gid, ';'.join(taxa)))
            else:
                self.logger.error(f'Unrecognized domain: {taxa[0]}')
                sys.exit(-1)

        fout_ar.close()
        fout_bac.close()

    def rep_ssu_files(self, reps, genome_paths, min_16S_ar_len, min_16S_bac_len, min_contig_len):
        """Create archaeal and bacterial SSU FASTA files for representative genomes."""

        self.logger.info(
            'Identifying SSU sequence for each representative genome.')
        fout_ar = open(self.output_dir /
                       'ar122_ssu_r{}.tsv'.format(self.release_number), 'w')
        fout_bac = open(self.output_dir /
                        'bac120_ssu_r{}.tsv'.format(self.release_number), 'w')

        failed_ssu_len = 0
        failed_contig_len = 0
        failed_max_ssu_len = 0
        valid_ssu_seqs = 0
        for gid, taxa in reps.items():
            if taxa[0] == 'd__Archaea':
                fout = fout_ar
                min_gene_len = min_16S_ar_len
            elif taxa[0] == 'd__Bacteria':
                fout = fout_bac
                min_gene_len = min_16S_bac_len
            else:
                self.logger.error(f'Unrecognized domain: {taxa[0]}')
                sys.exit(-1)

            genome_path = Path(genome_paths[gid])
            ssu_dir = genome_path / SSU_DIR
            if not ssu_dir.exists():
                self.logger.warning(f'SSU directory does not exist: {ssu_dir}')

            hmm_summary = ssu_dir / 'ssu.hmm_summary.tsv'
            if not hmm_summary.exists():
                continue

            ssu_fasta = ssu_dir / 'ssu.fna'
            seqs = seq_io.read(ssu_fasta.as_posix())

            best_ssu_len = 0
            best_contig_id = None
            best_hit_str = None
            with open(hmm_summary) as f:
                header = f.readline()
                for line in f:
                    line_split = line.strip().split('\t')
                    contig_id = line_split[0]
                    evalue = float(line_split[2])
                    start_hit = int(line_split[3])
                    end_hit = int(line_split[4])
                    ssu_len = int(line_split[5])
                    contig_len = int(line_split[-1])

                    if ssu_len < min_gene_len:
                        failed_ssu_len += 1

                    if contig_len < min_contig_len:
                        failed_contig_len += 1

                    if contig_len > MAX_SSU_LEN:
                        failed_max_ssu_len += 1

                    if ssu_len >= min_gene_len and contig_len >= min_contig_len:
                        if best_ssu_len < ssu_len:
                            # get longest identified SSU sequence
                            best_ssu_len = ssu_len
                            best_contig_id = contig_id
                            best_hit_str = '>{} {} [locus_tag={}] [location={}..{}] [ssu_len={}] [contig_len={}]\n'.format(
                                gid,
                                ';'.join(reps[gid]),
                                contig_id,
                                start_hit,
                                end_hit,
                                ssu_len,
                                contig_len)

            if best_hit_str:
                fout.write(best_hit_str)
                fout.write(f'{seqs[best_contig_id]}\n')
                valid_ssu_seqs += 1

        fout_ar.close()
        fout_bac.close()

        self.logger.info(f'Filtered {failed_ssu_len:,} SSU sequence based on rRNA gene length.')
        self.logger.info(f'Filtered {failed_contig_len:,} SSU sequence based on contig length.')
        self.logger.info(f'Filtered {failed_max_ssu_len:,} SSU sequence based on a maximum allowed length of {MAX_SSU_LEN} bp.')
        self.logger.info(f'Identified {valid_ssu_seqs:,} genomes with a valid SSU sequence.')

    def gtdb_ssu_file(self, reps, sp_clusters, genome_paths):
        """Create SSU FASTA file spanning all genomes in a GTDB species cluster."""

        self.logger.info(
            'Creating FASTA file containing all identified SSU sequences.')
        fout = open(self.output_dir /
                    'ssu_r{}.tsv'.format(self.release_number), 'w')

        gid_rep = {}
        for rid, cids in sp_clusters.items():
            for cid in cids:
                gid_rep[cid] = rid

        for gids in sp_clusters.values():
            for gid in gids:
                genome_path = Path(genome_paths[gid])
                ssu_dir = genome_path / SSU_DIR
                if not ssu_dir.exists():
                    self.logger.warning(f'SSU directory does not exist: {ssu_dir}')

                hmm_summary = ssu_dir / 'ssu.hmm_summary.tsv'
                if not hmm_summary.exists():
                    continue

                ssu_fasta = ssu_dir / 'ssu.fna'
                seqs = seq_io.read(ssu_fasta.as_posix())

                with open(hmm_summary) as f:
                    header = f.readline()
                    for line in f:
                        line_split = line.strip().split('\t')
                        contig_id = line_split[0]
                        evalue = float(line_split[2])
                        start_hit = int(line_split[3])
                        end_hit = int(line_split[4])
                        ssu_len = int(line_split[5])
                        contig_len = int(line_split[-1])

                        assert(not gid.startswith('U_'))
                        fout.write('>{}~{} {} [location={}..{}] [ssu_len={}] [contig_len={}]\n'.format(
                            gid,
                            contig_id,
                            ';'.join(reps[gid_rep[gid]]),
                            start_hit,
                            end_hit,
                            ssu_len,
                            contig_len))
                        fout.write(f'{seqs[contig_id]}\n')

        fout.close()

    def ssu_files(self, metadata_file,
                  gtdb_sp_clusters_file,
                  user_gid_table,
                  genome_path_file,
                  min_16S_ar_len,
                  min_16S_bac_len,
                  min_contig_len):
        """Generate 16S rRNA FASTA files for GTDB website."""

        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file, user_gids)
        self.logger.info(' ...identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))

        # get GTDB species clusters
        self.logger.info('Parsing all species clusters from GTDB metadata.')
        sp_clusters = parse_species_clusters(gtdb_sp_clusters_file, user_gids)
        self.logger.info(' ...identified {:,} representative genomes spanning {:,} genomes.'.format(
            len(sp_clusters),
            sum([len(cid) for cid in sp_clusters.values()])))

        # read path to data directory for each genome
        self.logger.info('Parsing path to data directory for each genome.')
        genome_paths = parse_genomic_path_file(genome_path_file, user_gids)

        # create archaeal and bacterial SSU FASTA file for representative
        # genomes
        self.rep_ssu_files(reps, genome_paths, min_16S_ar_len,
                           min_16S_bac_len, min_contig_len)

        # create SSU FASTA file spanning all genomes in a GTDB species cluster
        self.gtdb_ssu_file(reps, sp_clusters, genome_paths)

    def sp_cluster_file(self, metadata_file,
                        gtdb_sp_clusters_file,
                        user_gid_table):
        """Generate file indicating GTDB species clusters."""

        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file, user_gids)
        self.logger.info(' ...identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))

        # write out file with GTDB species clusters
        fout = open(self.output_dir /
                    'sp_clusters_r{}.tsv'.format(self.release_number), 'w')
        fout.write(
            'Representative genome\tGTDB species\tGTDB taxonomy\tANI circumscription radius')
        fout.write('\tMean intra-species ANI\tMin intra-species ANI')
        fout.write('\tMean intra-species AF\tMin intra-species AF')
        fout.write('\tNo. clustered genomes\tClustered genomes\n')

        with open(gtdb_sp_clusters_file) as f:
            header = f.readline().strip().split('\t')

            type_genome_index = header.index('Type genome')
            ani_index = header.index('ANI radius')
            cluster_size_index = header.index('No. clustered genomes')
            mean_ani_index = header.index('Mean ANI')
            min_ani_index = header.index('Min ANI')
            mean_af_index = header.index('Mean AF')
            min_ani_index = header.index('Min AF')
            clustered_genomes_index = header.index('Clustered genomes')

            for line in f:
                line_split = line.strip().split('\t')

                type_genome = line_split[type_genome_index]
                type_genome = user_gids.get(type_genome, type_genome)

                ani = line_split[ani_index]
                cluster_size = int(line_split[cluster_size_index])
                mean_ani = line_split[mean_ani_index]
                min_ani = line_split[min_ani_index]
                mean_af = line_split[mean_af_index]
                min_af = line_split[min_ani_index]

                clustered_genomes = set([type_genome])
                if cluster_size > 0:
                    for cid in line_split[clustered_genomes_index].split(','):
                        clustered_genomes.add(user_gids.get(cid, cid))

                fout.write('{}\t{}\t{}\t{}'.format(
                    type_genome,
                    reps[type_genome][6],
                    ';'.join(reps[type_genome][0:6]),
                    ani))
                fout.write('\t{}\t{}'.format(mean_ani, min_ani))
                fout.write('\t{}\t{}'.format(mean_af, min_af))
                fout.write('\t{}\t{}\n'.format(cluster_size + 1,
                                               ','.join(sorted(clustered_genomes))))

                assert(not type_genome.startswith('U_'))
                assert(cluster_size + 1 == len(clustered_genomes))

        fout.close()

    def tree_files(self,
                   bac_tree,
                   ar_tree,
                   user_gid_table):
        """Generate tree files spanning representative genomes."""

        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)

        # ensure proper genome IDs for bacterial tree
        self.logger.info('Reading GTDB bacterial reference tree.')
        bac_tree = dendropy.Tree.get_from_path(bac_tree,
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)
        self.logger.info(' ...tree contains {:,} genomes.'.format(
            sum([1 for leaf in bac_tree.leaf_node_iter()])))

        for leaf in bac_tree.leaf_node_iter():
            gid = leaf.taxon.label
            assert(not gid.startswith('D-'))
            leaf.taxon.label = user_gids.get(gid, gid)

        bac_tree.write_to_path(self.output_dir / f'bac120_r{self.release_number}.tree',
                               schema='newick',
                               suppress_rooting=True,
                               unquoted_underscores=True)

        # ensure proper genome IDs for archaeal tree
        self.logger.info('Reading GTDB archaeal reference tree.')
        ar_tree = dendropy.Tree.get_from_path(ar_tree,
                                              schema='newick',
                                              rooting='force-rooted',
                                              preserve_underscores=True)
        self.logger.info(' ...tree contains {:,} genomes.'.format(
            sum([1 for leaf in ar_tree.leaf_node_iter()])))

        for leaf in ar_tree.leaf_node_iter():
            gid = leaf.taxon.label
            assert(not gid.startswith('D-'))
            leaf.taxon.label = user_gids.get(gid, gid)

        ar_tree.write_to_path(self.output_dir / f'ar122_r{self.release_number}.tree',
                              schema='newick',
                              suppress_rooting=True,
                              unquoted_underscores=True)
                              
    def marker_files(self, bac120_gene_dir,
                            ar122_gene_dir,
                            user_gid_table):
        """Generate marker gene file."""
        
        # copy and rename marker info files
        self.logger.info('Copying and renaming marker info files.')
        shutil.copyfile(PurePath(bac120_gene_dir) / 'bac120_markers_info.tsv',
                        self.output_dir / f'bac120_msa_marker_info_r{self.release_number}.tsv')
        shutil.copyfile(PurePath(ar122_gene_dir) / 'ar122_markers_info.tsv',
                        self.output_dir / f'ar122_msa_marker_info_r{self.release_number}.tsv')
                        
        # compress individual marker gene alignments
        self.logger.info('Creating compressed files containing all marker gene MSAs.')
        for msa_dir, prefix in [(ar122_gene_dir, 'ar122'), (bac120_gene_dir, 'bac120')]:
            self.logger.info(f'Reading {msa_dir}.')
            msa_file = self.output_dir / f'{prefix}_msa_individual_genes_r{self.release_number}.tar.gz'
            
            msa_files = set()
            for f in os.listdir(msa_dir):
                if f.endswith('.faa') and '_concatenated' not in f:
                    msa_files.add(f)
            
            cmd = 'tar -zcvf {}_msa_individual_genes_r{}.tar.gz --directory={} {}'.format(
                    prefix,
                    self.release_number,
                    msa_dir, 
                    ' '.join(msa_files))
            os.system(cmd)
            
    def msa_files(self, bac120_msa_file,
                        ar122_msa_file,
                        metadata_file,
                        user_gid_table):
        """Generate concatenated MSA files."""
        
        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)
        
        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file, user_gids)
        self.logger.info(' ...identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
                            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
                            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))
         
        # create MSA files for representative genomes with modified genome IDs
        self.logger.info('Creating MSA files.')
        for msa_file, prefix in [(ar122_msa_file, 'ar122'), (bac120_msa_file, 'bac120')]:
            out_file = f'{prefix}_msa_r{self.release_number}.faa'
            fout = open(self.output_dir / out_file, 'w')
            count = 0
            for seq_id,seq in seq_io.read_seq(msa_file):
                seq_id = user_gids.get(seq_id, seq_id)
                if seq_id in reps:
                    fout.write(f'>{seq_id}\n')
                    fout.write(f'{seq}\n')
                    count += 1
                    
            self.logger.info(f'MSA file {out_file} contains {count:,} genomes.')
                    
            fout.close()
            
    def metadata_files(self, metadata_file,
                        metadata_fields,
                        gtdb_sp_clusters_file,
                        user_gid_table):
        """Generate GTDB metadata files."""
        
        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)
        
        # get GTDB species clusters
        self.logger.info('Parsing all species clusters from GTDB metadata.')
        sp_clusters = parse_species_clusters(gtdb_sp_clusters_file, user_gids)
        self.logger.info(' ...identified {:,} representative genomes spanning {:,} genomes.'.format(
            len(sp_clusters),
            sum([len(cid) for cid in sp_clusters.values()])))
            
        all_gids = set()
        for gids in sp_clusters.values():
            all_gids.update(gids)
        
        # read metadata fields to retain
        self.logger.info('Reading metadata fields to retain.')
        fields = set()
        with open(metadata_fields) as f:
            f.readline()
            for line in f:
                fields.add(line.split('\t')[0])
        self.logger.info(f' ...identified {len(fields):,} fields to retain.')
        
        # get list of all genome IDs in metadata file
        metadata_gids = set()
        with open(metadata_file, encoding='utf-8') as f:
            f.readline()
            for line in f:
                line_split = line.strip().split('\t')
                metadata_gids.add(line_split[0])
        
        # write out bacterial and archaeal metadata files
        self.logger.info('Creating domain-specific metadata fields.')
        fout_bac = open(self.output_dir / f'bac120_metadata_r{self.release_number}.tsv', 'w',
                        encoding='utf-8')
        fout_ar = open(self.output_dir / f'ar122_metadata_r{self.release_number}.tsv', 'w',
                        encoding='utf-8')
        
        fields_sorted = sorted(fields)
        fields_sorted.remove('accession')
        fout_bac.write('{}\t{}\n'.format('accession', '\t'.join(fields_sorted)))
        fout_ar.write('{}\t{}\n'.format('accession', '\t'.join(fields_sorted)))
        
        count = 0
        with open(metadata_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            
            field_index = {}
            for field in fields_sorted:
                field_index[field] = header.index(field)
                
            for line in f:
                line_split = line.strip().split('\t')
                
                gid = line_split[0]
                new_gid = user_gids.get(gid, gid)
                if gid != new_gid and new_gid in metadata_gids:
                    # skip this record since they GTDB user genome is already
                    # in the metadata file under the NCBI accession number
                    continue
                
                if new_gid in all_gids:
                    count += 1
                    gtdb_taxonomy = line_split[field_index['gtdb_taxonomy']]
                    gtdb_taxa = [t.strip() for t in gtdb_taxonomy.split(';')]
                    gtdb_domain = gtdb_taxa[0]
                    
                    if gtdb_domain == 'd__Bacteria':
                        fout = fout_bac
                    elif gtdb_domain == 'd__Archaea':
                        fout = fout_ar
                    else:
                        self.logger.error(f'Unrecognized GTDB domain: {gtdb_domain} {gid}')
                        
                    fout.write(new_gid)
                    for field in fields_sorted:
                        if field == 'gtdb_genome_representative':
                            gid = line_split[field_index[field]]
                            new_gid = user_gids.get(gid, gid)
                            fout.write('\t{}'.format(new_gid))
                        else:
                            fout.write('\t{}'.format(line_split[field_index[field]]))
                    fout.write('\n')
                         
        fout_bac.close()
        fout_ar.close()
        
        self.logger.info(f'Wrote metadata for {count:,} genomes.')
        
    def validate(self,
                    taxonomy_file,
                    tree_file,
                    metadata_file,
                    msa_file,
                    ssu_file,
                    sp_clusters_file):
        """Perform validation checks on GTDB website files."""
        
        # parse taxonomy file
        self.logger.info('Parsing taxonomy file.')
        taxonomy = {}
        for line in open(taxonomy_file):
            line_split = line.strip().split('\t')
            taxonomy[line_split[0]] = line_split[1]
            
        # get GTDB species clusters
        self.logger.info('Validating species clusters.')
        sp_clusters = parse_species_clusters(sp_clusters_file, {})
        assert(len(set(taxonomy).difference(sp_clusters)) == 0)
            
        # validate that tree spans representatives
        self.logger.info('Validating tree file.')
        tree = dendropy.Tree.get_from_path(tree_file,
                                               schema='newick',
                                               rooting='force-rooted',
                                               preserve_underscores=True)
        tree_gids = set([leaf.taxon.label for leaf in tree.leaf_node_iter()])
        assert(len(tree_gids.symmetric_difference(taxonomy)) == 0)
        
        # validate that metadata file contains all representative genomes
        self.logger.info('Validating metadata file.')
        metadata_gids = set()
        with open(metadata_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            gtdb_taxonomy_index = header.index('gtdb_taxonomy')
            gtdb_rep_index = header.index('gtdb_representative')
            gtdb_rep_genome_index = header.index('gtdb_genome_representative')
            
            for line in f:
                line_split = line.strip().split('\t')
                gid = line_split[0]
                metadata_gids.add(gid)
                
                gtdb_rep = line_split[gtdb_rep_index]
                if ((gid in sp_clusters and gtdb_rep == 'f')
                    or (gid not in sp_clusters and gtdb_rep == 't')):
                    self.logger.info('GTDB species representative information is incorrect: {gid}')
                    sys.exit(-1)
                    
                rep_id = line_split[gtdb_rep_genome_index]
                assert(taxonomy[rep_id] == line_split[gtdb_taxonomy_index])
                
                if gid in taxonomy:
                    if taxonomy[gid] != line_split[gtdb_taxonomy_index]:
                        self.logger.error(f'GTDB taxonomy is incorrect in metadata file: {gid}')
                        sys.exit(-1)
                
        assert(len(set(taxonomy).difference(metadata_gids)) == 0)
        
        # validate that MSA file spans representatives
        self.logger.info('Validating MSA file.')
        msa_seqs = seq_io.read(msa_file)
        assert(len(set(msa_seqs).symmetric_difference(taxonomy)) == 0)
        
        # validate that SSU file spans representatives
        self.logger.info('Validating SSU file.')
        ssu_seqs = set()
        for seq_id, seq, annotations in seq_io.read_seq(ssu_file, keep_annotation=True):
            ssu_taxonomy = annotations[0:annotations.find('[')].strip()
            if ssu_taxonomy != taxonomy[seq_id]:
                self.logger.error(f'GTDB taxonomy is incorrect in SSU file: {gid}')
            ssu_seqs.add(seq_id)
        
        assert(len(ssu_seqs.difference(taxonomy)) == 0)
        
        self.logger.info('Passed validation.')
                         
    def tax_comp_files(self, old_taxonomy, new_taxonomy, gid_table, keep_no_change, keep_placeholder_name, top_change):

        gtogfile = open(os.path.join(
            self.output_dir, 'genomes_summary.tsv'), "w")
        gtogfile.write("genome_id\tfirst_taxonomy\tsecond_taxonomy\tstatus\n")
        for dom in ["d__Archaea", "d__Bacteria"]:
            tempdict = {}  # temporary dictionary for the information in old taxonomy
            olddict = {}
            newdict = {}
            mapdict = {}  # dictionary of uba_mappings
            oldnumber = 0  # number of genes in old_tax
            newnumber = 0  # number of genes in new_tax
            intersection = 0  # number of genes common to both taxonomies

            with open(gid_table, "r") as tax_map:
                for line in tax_map:
                    genes_id = line.strip('\n').split("\t")
                    # IDs that start with U_
                    u_id = genes_id[0]
                    # IDs that start with GCA or GCF
                    gc_id = genes_id[2]
                    # delete GB or RS prefix to ensure id-consistency
                    if gc_id.startswith("GB") or gc_id.startswith("RS"):
                        gc_id = gc_id[3:]
                    # map k, v of all id_mappings
                    mapdict[gc_id] = u_id
                    mapdict[u_id] = gc_id

            with open(new_taxonomy) as taxnew:
                for line in taxnew:
                    newnumber += 1
                    genome_info = line.strip().split("\t")
                    gene_id = genome_info[0]
                    ranks = genome_info[1].split(";")
                    if gene_id.startswith("GB") or gene_id.startswith("RS"):
                        gene_id = gene_id[3:]
                    if gene_id.startswith("GC") or (gene_id.startswith("U") and gene_id in mapdict):
                        # temporary dictionary of IDs of interest in
                        # new_taxonomy
                        tempdict[gene_id] = {"d": ranks[0], "p": ranks[1], "c": ranks[2], "o": ranks[3],
                                             "f": ranks[4], "g": ranks[5], "s": ranks[6], 'full': genome_info[1]}
                        if tempdict.get(gene_id).get("d") != dom:
                            tempdict.pop(gene_id)

            with open(old_taxonomy) as taxold:
                for line in taxold:
                    oldnumber += 1
                    genome_info = line.strip().split("\t")
                    gene_id = genome_info[0]
                    ranks = genome_info[1].split(";")
                    if gene_id.startswith("GB") or gene_id.startswith("RS"):
                        gene_id = gene_id[3:]
                    if gene_id.startswith("GC"):
                        if gene_id in tempdict:
                            intersection += 1
                            newdict[gene_id] = tempdict.get(gene_id)
                            olddict[gene_id] = {"d": ranks[0], "p": ranks[1], "c": ranks[2], "o": ranks[3],
                                                "f": ranks[4], "g": ranks[5], "s": ranks[6], 'full': genome_info[1]}
                        else:
                            if gene_id in mapdict:
                                other_id = mapdict.get(gene_id)
                                if other_id in tempdict:
                                    intersection += 1
                                    newdict[gene_id] = tempdict.get(other_id)
                                    olddict[gene_id] = {"d": ranks[0], "p": ranks[1], "c": ranks[2], "o": ranks[3],
                                                        "f": ranks[4], "g": ranks[5], "s": ranks[6], 'full': genome_info[1]}
                    else:
                        if gene_id in mapdict:
                            other_id = mapdict.get(gene_id)
                            if other_id in tempdict:
                                intersection += 1
                                newdict[gene_id] = tempdict.get(other_id)
                                olddict[gene_id] = {"d": ranks[0], "p": ranks[1], "c": ranks[2], "o": ranks[3],
                                                    "f": ranks[4], "g": ranks[5], "s": ranks[6], 'full': genome_info[1]}
                        else:
                            if gene_id in tempdict:
                                intersection += 1
                                newdict[gene_id] = tempdict.get(gene_id)
                                olddict[gene_id] = {"d": ranks[0], "p": ranks[1], "c": ranks[2], "o": ranks[3],
                                                    "f": ranks[4], "g": ranks[5], "s": ranks[6], 'full': genome_info[1]}

            pref = 'bac'
            if dom == "d__Archaea":
                pref = "arc"

            print("Analysis of {0}".format(dom))
            print("Number of genomes in old taxonomy: {0}".format(oldnumber))
            print("Number of genomes in new taxonomy: {0}".format(newnumber))
            print("Number of genomes present in both taxonomy: {0}".format(
                intersection))
            print("Number of genomes in new taxonomy, unavailable in old taxonomy: {0}".format(
                newnumber - intersection))
            print("Number of genomes in old taxonomy, unavailable in new taxonomy: {0}".format(
                oldnumber - intersection))

            self.genomebygenomecompare(olddict, newdict, gtogfile)

            print("- Analysis of phyla")
            self.rundistrib(olddict, newdict, "p", os.path.join(self.output_dir, pref +
                                                                "_phylum_difference.tsv"), "phylum", keep_placeholder_name, keep_no_change, top_change)
            print("- Analysis of classes")
            self.rundistrib(olddict, newdict, "c", os.path.join(self.output_dir, pref +
                                                                "_class_difference.tsv"), "class", keep_placeholder_name, keep_no_change, top_change)
            print("- Analysis of orders")
            self.rundistrib(olddict, newdict, "o", os.path.join(self.output_dir, pref +
                                                                "_order_difference.tsv"), "order", keep_placeholder_name, keep_no_change, top_change)
            print("- Analysis of families")
            self.rundistrib(olddict, newdict, "f", os.path.join(self.output_dir, pref +
                                                                "_family_difference.tsv"), "family", keep_placeholder_name, keep_no_change, top_change)
            print("- Analysis of genera")
            self.rundistrib(olddict, newdict, "g", os.path.join(self.output_dir, pref +
                                                                "_genus_difference.tsv"), "genus", keep_placeholder_name, keep_no_change, top_change)

        gtogfile.close()
        
    def json_tree_file(self, taxonomy_file, metadata_file, gid_table):
        user_gids = parse_user_gid_table(gid_table)

    def genomebygenomecompare(self, olddict, newdict, gtogfile):
        orderrank = ["domain", "phylum", "class",
                     'order', 'family', 'genus', 'species']
        for k, v in newdict.items():
            if k in olddict:
                if v.get('full') == olddict.get(k).get('full'):
                    gtogfile.write("{0}\t{1}\t{2}\tidentical\n".format(
                        k, olddict.get(k).get('full'), v.get('full')))
                else:
                    for rk in orderrank:
                        if v.get(rk[0:1]) != olddict.get(k).get(rk[0:1]):
                            gtogfile.write("{0}\t{1}\t{2}\t{3} difference\n".format(
                                k, olddict.get(k).get('full'), v.get('full'), rk))
            else:
                gtogfile.write(
                    "{0}\tnot_available\t{1}\t\n".format(k, v.get('full')))
        for k, v in olddict.items():
            if k not in newdict:
                gtogfile.write("{0}\t{1}\tnot_available\t\n".format(
                    k, olddict.get(k).get('full')))

    def rundistrib(self, olddict, newdict, rank, outfile, ref, keep_placeholder_name, keep_no_change, top_change):
        orderrank = ["d", "p", "c", 'o', 'f', 'g', 's']
        different_ranks = []
        # take the different ranks in the dictionaries and then count the appearances of each pair
        # of ranks ie. (p__Fusobacteria, p_Fusobacteriota),200
        for k, v in newdict.items():
            if k in olddict:
                different_ranks.append((olddict.get(k).get(rank), v.get(rank)))
        order_counter = Counter(different_ranks).most_common()
        results_dict = {}
        """create a dictionary with the name of the rank in olddict as key and the number of 
        genomes belonging to that rank encountered, along with different rank names on newdict
        for those genomes and how many times the different names appear."""
        for item in order_counter:
            if item[0][0] in results_dict:
                results_dict[item[0][0]]["nber_g"] = results_dict[item[0]
                                                                  [0]]["nber_g"] + item[1]
                results_dict[item[0][0]].get(
                    "genomes").append((item[0][1], item[1]))
            else:
                results_dict[item[0][0]] = {
                    "nber_g": item[1], "genomes": [(item[0][1], item[1])]}

        to_pop = []  # store taxonomies to take out of the file
        for k, v in results_dict.items():
            if top_change == True:
                top_gen = results_dict[k]['genomes'][0][0]
                if top_gen == k:
                    to_pop.append(k)
                else:
                    results_dict[k]['genomes'] = [
                        results_dict[k]['genomes'][0]]
            if keep_no_change == True and k == results_dict[k]['genomes'][0][0] and results_dict[k]['nber_g'] == results_dict[k]['genomes'][0][1]:
                to_pop.append(k)
            if keep_placeholder_name == True and self.hasNumber(k) == True and results_dict[k]['genomes'][0][0] == rank + "__":
                to_pop.append(k)
            elif keep_placeholder_name == True and self.hasNumber(results_dict[k]['genomes'][0][0]) == True and k == rank + "__":
                to_pop.append(k)
        for k in to_pop:
            results_dict.pop(k)

        outf = open(outfile, "w")
        outf.write(
            "Previous {0}\tNumber of genomes\tNumber of new {0}\tList of new {0}\n".format(ref))
        for k, v in results_dict.items():
            for sk, sv in olddict.items():
                if sv.get(rank) == k:
                    prankold = sv.get(orderrank[orderrank.index(rank) - 1])
                    break
            number_sub = len(v.get("genomes"))

            results = []

            if len(v.get("genomes")) == 1 and v.get("genomes")[0][0] == k:
                continue

            for newg in v.get("genomes"):
                newg_name = newg[0]
                newg_numb = float(newg[1])
                for sk, sv in newdict.items():
                    if sv.get(rank) == newg_name:
                        pranknew = sv.get(orderrank[orderrank.index(rank) - 1])
                        break
                res = newg_name
                if prankold != pranknew:
                    res = res + "(" + pranknew + ")"
                res = "{0} {1}%".format(res, round(
                    (newg_numb / v.get("nber_g")) * 100, 2))
                results.append(res)

            outf.write("{0}\t{1}\t{2}\t{3}\n".format(
                k, v.get("nber_g"), number_sub, ", ".join(results)))
        outf.close()
