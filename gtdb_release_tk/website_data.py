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
        fout_ar = open(self.output_dir / 'ar122_taxonomy_r{}.tsv'.format(self.release_number), 'w')
        fout_bac = open(self.output_dir / 'bac120_taxonomy_r{}.tsv'.format(self.release_number), 'w')
        
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
        
        self.logger.info('Identifying SSU sequence for each representative genome.')
        fout_ar = open(self.output_dir / 'ar122_ssu_r{}.tsv'.format(self.release_number), 'w')
        fout_bac = open(self.output_dir / 'bac120_ssu_r{}.tsv'.format(self.release_number), 'w')

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
        
        self.logger.info('Creating FASTA file containing all identified SSU sequences.')
        fout = open(self.output_dir / 'ssu_r{}.tsv'.format(self.release_number), 'w')
        
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
                            
        # create archaeal and bacterial SSU FASTA file for representative genomes
        self.rep_ssu_files(reps, genome_paths, min_16S_ar_len, min_16S_bac_len, min_contig_len)
        
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
        fout = open(self.output_dir / 'sp_clusters_r{}.tsv'.format(self.release_number), 'w')
        fout.write('Representative genome\tGTDB species\tGTDB taxonomy\tANI circumscription radius')
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
                fout.write('\t{}\t{}\n'.format(cluster_size+1, 
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
