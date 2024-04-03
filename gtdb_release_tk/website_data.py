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

import json
import logging
import ntpath
import os
import shutil
import multiprocessing as mp
import glob

from biolib.common import make_sure_path_exists, canonical_gid
from tqdm import tqdm
import sys
from collections import Counter
from collections import defaultdict

from pathlib import Path, PurePath

import biolib.seq_io as seq_io
import dendropy

from gtdb_release_tk.common import (parse_user_gid_table, parse_canonical_gid_table,
                                    parse_rep_genomes,
                                    parse_genomic_path_file,
                                    parse_species_clusters,
                                    parse_gtdb_metadata, parse_tophit_file, parse_taxonomy_file)
from gtdb_release_tk.config import SSU_DIR, MAX_SSU_LEN, BAC120_MARKERS, AR53_MARKERS


class WebsiteData(object):
    """Generate data files for website."""

    def __init__(self, release_number, output_dir):
        """Initialization."""

        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def taxonomy_files(self,
                       metadata_file,
                       gtdb_sp_clusters_file,
                       only_reps):
        """Generate taxonomy files for GTDB website."""

        # get GTDB species clusters
        self.logger.info('Parsing all species clusters from GTDB metadata.')
        sp_clusters = parse_species_clusters(gtdb_sp_clusters_file)

        genome_rep_id = {}
        for rep_id, gids in sp_clusters.items():
            for gid in gids:
                genome_rep_id[gid] = rep_id

        self.logger.info(' ...identified {:,} representative genomes spanning {:,} genomes.'.format(
            len(sp_clusters),
            len(genome_rep_id)))

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)
        self.logger.info(' ...identified {:,} representative genomes ({:,} bacterial, {:,} archaeal).'.format(
            sum([1 for t in reps.values()]),
            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))

        # create archaeal and bacterial taxonomy files
        if only_reps:
            fout_ar = open(
                self.output_dir / 'ar53_taxonomy_r{}_reps.tsv'.format(self.release_number), 'w')
            fout_bac = open(
                self.output_dir / 'bac120_taxonomy_r{}_reps.tsv'.format(self.release_number), 'w')
        else:
            fout_ar = open(
                self.output_dir / 'ar53_taxonomy_r{}.tsv'.format(self.release_number), 'w')
            fout_bac = open(
                self.output_dir / 'bac120_taxonomy_r{}.tsv'.format(self.release_number), 'w')

        for gid in genome_rep_id:
            assert(not gid.startswith('U_'))

            rep_id = genome_rep_id[gid]
            taxa = reps[rep_id]

            if only_reps and gid not in reps:
                continue

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
                       'ar53_ssu_reps_r{}.fna'.format(self.release_number), 'w')
        fout_bac = open(self.output_dir /
                        'bac120_ssu_reps_r{}.fna'.format(self.release_number), 'w')

        failed_ssu_len = 0
        failed_contig_len = 0
        failed_max_ssu_len = 0
        valid_ssu_seqs = 0
        for gid, taxa in tqdm(reps.items()):
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
                    'ssu_all_r{}.fna'.format(self.release_number), 'w')

        gid_rep = {}
        for rid, cids in sp_clusters.items():
            for cid in cids:
                gid_rep[cid] = rid

        # initialising progress bar objects
        high_loop = tqdm(range(1), bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:100}{r_bar}')
        low_loop = tqdm(range(1), bar_format='{desc:<5.5}{percentage:3.0f}%|{bar:100}{r_bar}')

        high_loop.refresh()
        high_loop.reset(total=len(sp_clusters))

        for gids in sp_clusters.values():
            high_loop.update()

            low_loop.refresh()
            low_loop.reset(total=len(gids))
            for gid in gids:
                low_loop.update()
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
                  genome_path_file,
                  min_16S_ar_len,
                  min_16S_bac_len,
                  min_contig_len):
        """Generate 16S rRNA FASTA files for GTDB website."""

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)
        self.logger.info(' ...identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))

        # get GTDB species clusters
        self.logger.info('Parsing all species clusters from GTDB metadata.')
        sp_clusters = parse_species_clusters(gtdb_sp_clusters_file)
        self.logger.info(' ...identified {:,} representative genomes spanning {:,} genomes.'.format(
            len(sp_clusters),
            sum([len(cid) for cid in sp_clusters.values()])))

        # read path to data directory for each genome
        self.logger.info('Parsing path to data directory for each genome.')
        genome_paths = parse_genomic_path_file(genome_path_file)

        # create archaeal and bacterial SSU FASTA file for representative
        # genomes
        self.rep_ssu_files(reps, genome_paths, min_16S_ar_len,
                           min_16S_bac_len, min_contig_len)

        # create SSU FASTA file spanning all genomes in a GTDB species cluster
        self.gtdb_ssu_file(reps, sp_clusters, genome_paths)

    def sp_cluster_file(self, metadata_file,
                        gtdb_sp_clusters_file,
                        canonical_gid_table):
        """Generate file indicating GTDB species clusters."""

        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')

        # parse canonical ID mapping table
        self.logger.info('Parsing canonical ID mapping table.')
        canonical_gids = parse_canonical_gid_table(canonical_gid_table)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)
        self.logger.info(' - identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
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

            type_genome_index = header.index('Representative')
            ani_index = header.index('ANI radius')
            cluster_size_index = header.index('No. clustered genomes')
            mean_ani_index = header.index('Mean ANI')
            min_ani_index = header.index('Min ANI')
            mean_af_index = header.index('Mean AF')
            min_af_index = header.index('Min AF')
            clustered_genomes_index = header.index('Clustered genomes')

            for line in f:
                line_split = line.strip().split('\t')

                type_genome = line_split[type_genome_index]
                type_genome = canonical_gids.get(type_genome, type_genome)

                ani = line_split[ani_index]
                cluster_size = int(line_split[cluster_size_index])
                mean_ani = line_split[mean_ani_index]
                min_ani = line_split[min_ani_index]
                mean_af = line_split[mean_af_index]
                min_af = line_split[min_af_index]

                clustered_genomes = set([type_genome])
                if cluster_size > 0:
                    for cid in line_split[clustered_genomes_index].split(','):
                        cid = canonical_gids.get(cid, cid)
                        clustered_genomes.add(cid)

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

    def hq_genome_file(self,metadata_file):
        """Generate file indicating HQ genomes."""


        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)
        self.logger.info(' - identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))

        # get metadata required to access genome quality
        metadata = parse_gtdb_metadata(metadata_file, [
            'checkm_completeness',
            'checkm_contamination',
            'ssu_length',
            'lsu_23s_length',
            'lsu_5s_length',
            'trna_aa_count',
            'ncbi_genome_category',
            'mimag_high_quality',
            'mimag_medium_quality',
            'mimag_low_quality',
            'gtdb_taxonomy'])

        fout = open(self.output_dir /
                    'hq_mimag_genomes_r{}.tsv'.format(self.release_number), 'w')
        fout.write('Representative genome')
        fout.write(
            '\tCompleteness (%)\tContamination (%)\t23S rRNA length\t5S rRNA length\t16S rRNA length\ttRNA count')
        fout.write(
            '\tGenome category\tGTDB taxonomy\tGTDB representative genome\n')

        genome_count = 0
        species = set()
        rep_count = 0
        for gid in metadata:
            if gid.startswith('U_'):
                continue

            gtdb_taxa = [t.strip()
                         for t in metadata[gid].gtdb_taxonomy.split(';')]
            gtdb_sp = gtdb_taxa[6]
            gtdb_domain = gtdb_taxa[0]
            if gtdb_sp == 's__':
                continue

            comp = metadata[gid].checkm_completeness
            cont = metadata[gid].checkm_contamination
            lsu_23s_length = metadata[gid].lsu_23s_length
            lsu_5s_length = metadata[gid].lsu_5s_length
            trna_aa_count = metadata[gid].trna_aa_count
            ssu_length = metadata[gid].ssu_length

            if gtdb_domain == 'd__Archaea':
                ssu_length_threshold = 900
            else:
                ssu_length_threshold = 1200

            genome_category = 'isolate'
            category_str = metadata[gid].ncbi_genome_category
            if category_str:
                if ('metagenome' in metadata[gid].ncbi_genome_category.lower()
                        or 'environment' in metadata[gid].ncbi_genome_category.lower()):
                    genome_category = 'MAG'
                elif 'single cell' in metadata[gid].ncbi_genome_category.lower():
                    genome_category = 'SAG'

            if (comp > 90
                and cont < 5
                and lsu_23s_length is not None and lsu_23s_length >= 1900
                and lsu_5s_length is not None and lsu_5s_length>= 80
                and ssu_length is not None and ssu_length >= ssu_length_threshold
                    and trna_aa_count >= 18):
                fout.write('{}'.format(gid))
                fout.write('\t{:.1f}\t{:.2f}\t{}\t{}\t{}\t{}'.format(
                    comp,
                    cont,
                    lsu_23s_length,
                    lsu_5s_length,
                    ssu_length,
                    trna_aa_count))
                fout.write('\t{}\t{}\t{}'.format(
                    genome_category,
                    metadata[gid].gtdb_taxonomy,
                    gid in reps))
                fout.write('\n')
                
                if not metadata[gid].mimag_high_quality:
                    self.logger.warning('This genome should be considered high quality, but the mimag_high_quality flag was set to False: {}'.format(gid))

                genome_count += 1
                species.add(gtdb_sp)
                rep_count += (gid in reps)
        fout.close()

        self.logger.info('Identified {:,} high-quality MIMAG genomes from {:,} species spanning {:,} GTDB representatives.'.format(
            genome_count,
            len(species),
            rep_count))

    def qc_file(self,qc_failed_file,canonical_gid_table):
        # parse canonical ID mapping table
        self.logger.info('Parsing canonical ID mapping table.')
        canonical_gids = parse_canonical_gid_table(canonical_gid_table)

        # write out file with GTDB species clusters
        fout = open(self.output_dir /
                    'qc_failed_r{}.tsv'.format(self.release_number), 'w')

        with open(qc_failed_file) as fin:
            header = fin.readline()
            fout.write(header)
            for line in fin:
                infos = line.strip().split('\t')
                if infos[0] in canonical_gids:
                    infos[0] = canonical_gids.get(infos[0])
                fout.write('\t'.join(infos)+'\n')

    def lpsn_urls(self,lpsn_species_file):
        species_dict = {}

        with open(lpsn_species_file, 'r') as slf:
            for line in slf:
                _letter, raw_spename, url = line.strip().split("\t")
                spename = raw_spename.replace('"', '').replace('Candidatus ', '')
                species_dict['s__' + spename] = url

        with open(os.path.join(self.output_dir,'species_urls_{}.json'.format(self.release_number)), 'w') as file:
            file.write(json.dumps(species_dict))



    def tree_files(self,
                   metadata_file,
                   bac_tree,
                   ar_tree,
                   canonical_gid_table):
        """Generate tree files spanning representative genomes."""

        # parse canonical ID mapping table
        self.logger.info('Parsing canonical ID mapping table.')
        canonical_gids = parse_canonical_gid_table(canonical_gid_table)

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
            gid = canonical_gids.get(gid, gid)
            leaf.taxon.label = gid

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
            gid = canonical_gids.get(gid, gid)
            leaf.taxon.label = gid

        ar_tree.write_to_path(self.output_dir / f'ar53_r{self.release_number}.tree',
                              schema='newick',
                              suppress_rooting=True,
                              unquoted_underscores=True)

        # create trees with species assignments as terminal labels
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        gtdb_reps = parse_rep_genomes(metadata_file)
        self.logger.info(' ...identified {:,} representative genomes ({:,} bacterial, {:,} archaeal).'.format(
            sum([1 for t in gtdb_reps.values()]),
            sum([1 for t in gtdb_reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in gtdb_reps.values() if t[0] == 'd__Archaea'])))

        self.logger.info('Creating trees with species labels.')
        for domain_tree, output_tree in [(bac_tree, f'bac120_r{self.release_number}.sp_labels.tree'),
                                         (ar_tree, f'ar53_r{self.release_number}.sp_labels.tree')]:
            for leaf in domain_tree.leaf_node_iter():
                gid = leaf.taxon.label
                gtdb_sp = gtdb_reps[gid][6]
                leaf.taxon.label = gtdb_sp

            domain_tree.write_to_path(self.output_dir / output_tree,
                                      schema='newick',
                                      suppress_rooting=True,
                                      unquoted_underscores=True)

    def marker_files(self, bac120_gene_dir,
                     ar53_gene_dir):
        """Generate marker gene file."""


        # copy and rename marker info files
        self.logger.info('Copying and renaming marker info files.')
        shutil.copyfile(PurePath(bac120_gene_dir) / 'bac120_markers_info.tsv',
                        self.output_dir / f'bac120_msa_marker_info_r{self.release_number}.tsv')
        shutil.copyfile(PurePath(ar53_gene_dir) / 'ar53_markers_info.tsv',
                        self.output_dir / f'ar53_msa_marker_info_r{self.release_number}.tsv')

        # compress individual marker gene alignments
        self.logger.info(
            'Creating compressed files containing all marker gene MSAs.')
        for msa_dir, prefix in [(ar53_gene_dir, 'ar53'), (bac120_gene_dir, 'bac120')]:
            self.logger.info(f'Reading {msa_dir}.')
            msa_file = self.output_dir / f'{prefix}_msa_marker_genes_r{self.release_number}.tar.gz'

            processed_sequence_path = os.path.join(self.output_dir, prefix)
            make_sure_path_exists(processed_sequence_path)

            msa_files = set()
            for f in os.listdir(msa_dir):
                if f.endswith('.faa') and '_concatenated' not in f:
                    fasta_seqs = seq_io.read_fasta(os.path.join(msa_dir, f))
                    with open(os.path.join(processed_sequence_path, os.path.basename(f)), 'w') as temp_file:
                        for k, seqs in fasta_seqs.items():
                            temp_file.write('>{}\n{}\n'.format(
                                k, seqs))

                    msa_files.add(f)

            cmd = 'tar -zcvf {}_msa_marker_genes_r{}.tar.gz --directory={} {}'.format(
                prefix,
                self.release_number,
                processed_sequence_path,
                ' '.join(msa_files))
            os.system(cmd)

    def msa_files(self, bac120_msa_file,
                  ar53_msa_file,
                  metadata_file,canonical_gid_table):
        """Generate concatenated MSA files."""

        # parse canonical ID mapping table
        self.logger.info('Parsing canonical ID mapping table.')
        canonical_gids = parse_canonical_gid_table(canonical_gid_table)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)
        self.logger.info(' ...identified {:,} bacterial and {:,} archaeal representative genomes.'.format(
            sum([1 for t in reps.values() if t[0] == 'd__Bacteria']),
            sum([1 for t in reps.values() if t[0] == 'd__Archaea'])))

        # create MSA files for representative genomes with modified genome IDs
        self.logger.info('Creating MSA files.')
        for msa_file, prefix in [(ar53_msa_file, 'ar53'), (bac120_msa_file, 'bac120')]:
            out_file = f'{prefix}_msa_r{self.release_number}.faa'
            fout = open(self.output_dir / out_file, 'w')
            count = 0
            for seq_id, seq in seq_io.read_seq(msa_file):
                seq_id = canonical_gids.get(seq_id,seq_id)
                if seq_id in reps:
                    fout.write(f'>{seq_id}\n')
                    fout.write(f'{seq}\n')
                    count += 1

            self.logger.info(f'MSA file {out_file} contains {count:,} genomes.')

            fout.close()

    def _gene_files_worker(self, q_worker, q_writer, q_results):
        while True:
            # Extract the next job
            job = q_worker.get(block=True)
            if job is None:
                break
            gid, domain, dir_base = job
            gid_no_header = gid.replace('RS_', '').replace('GB_', '')
            dir_prodigal = os.path.join(dir_base, 'prodigal')

            try:
                path_pfam_th = os.path.join(f'{dir_prodigal}',f'{gid_no_header}_pfam_lite_tophit.tsv.gz')
            except IndexError as error:
                print(dir_prodigal)
                print(os.path.join(f'{dir_prodigal}',f'{gid_no_header}_pfam_lite_tophit.tsv.gz'))
            path_tigr_th = os.path.join(f'{dir_prodigal}',f'{gid_no_header}_tigrfam_lite_tophit.tsv.gz')

            th_pfam = parse_tophit_file(path_pfam_th)
            th_tigr = parse_tophit_file(path_tigr_th)
            th_all = {**th_pfam, **th_tigr}


            # Determine which marker set to use
            if domain == 'd__Archaea':
                marker_set = AR53_MARKERS
            elif domain == 'd__Bacteria':
                marker_set = BAC120_MARKERS
            else:
                raise Exception('Unknown domain.')

            # Load the NT and AA files
            path_fna = os.path.join(f'{dir_prodigal}',f'{gid_no_header}_protein.fna.gz')
            path_faa = os.path.join(f'{dir_prodigal}',f'{gid_no_header}_protein.faa.gz')

            fna = seq_io.read_fasta(path_fna)
            faa = seq_io.read_fasta(path_faa)

            # Extract the sequences for each marker
            set_thall = set(th_all.keys())
            for marker in marker_set:
                if marker in set_thall:

                    sorted_hits = sorted(
                        th_all[marker], key=lambda x: (-x[2], x[1]))
                    best_gid = sorted_hits[0][0]

                    # Multihit check if seqs are the same, if not blanks
                    if len(th_all[marker]) > 1:

                        multi_seqs = {faa[x[0]] for x in th_all[marker]}

                        # Multiple hits with different genes - don't keep it
                        if len(multi_seqs) > 1:
                            pass
                        else:
                            gene_id, e_val, bit_score = sorted_hits[0]
                            q_results.put(
                                (gid, marker, best_gid, fna[gene_id], faa[gene_id], domain))

                    # Just a single hit - take it
                    else:
                        gene_id, e_val, bit_score = sorted_hits[0]
                        q_results.put(
                            (gid, marker, best_gid, fna[gene_id], faa[gene_id], domain))

            # Done
            q_writer.put(gid)

        return True

    def _gene_files_writer(self, q_writer, n_total):
        with tqdm(total=n_total) as pbar:
            while True:
                job = q_writer.get(block=True)
                if job is None:
                    break
                pbar.update(1)
        return True

    def gene_files(self, path_taxonomy, genome_dirs_file, metadata_file, cpus, only_reps=False):
        cpus = max(1, cpus)

        # Parse the metadata file
        self.logger.info('Parsing taxonomy file.')
        user_tax = parse_taxonomy_file(path_taxonomy)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)

        # Parse the genomic path
        self.logger.info(
            'Parsing genome directory files ( including UBA paths).')
        genome_dirs = parse_genomic_path_file(genome_dirs_file)

        # Setup the multiprocessing objects
        manager = mp.Manager()
        q_worker = manager.Queue()
        q_writer = manager.Queue()
        q_results = manager.Queue()

        # Create a queue of items to be processed
        self.logger.info('Creating a queue of jobs to be processed.')
        iterations = 0
        for gid, dir_base in user_tax.items():
            domain = dir_base.split(';')[0]
            if gid in genome_dirs:
                if only_reps and gid in reps:
                    q_worker.put((gid, domain, genome_dirs.get(gid)))
                    iterations += 1
                elif not only_reps:
                    q_worker.put((gid, domain, genome_dirs.get(gid)))
                    iterations += 1
        [q_worker.put(None) for _ in range(cpus)]



        # Create the workers
        self.logger.info('Starting jobs.')
        p_workers = [mp.Process(target=self._gene_files_worker,
                                args=(q_worker, q_writer, q_results))
                     for _ in range(cpus)]
        p_writer = mp.Process(target=self._gene_files_writer,
                              args=(q_writer, iterations))

        # Start each of the threads.
        try:
            # Start the writer and each processing thread.
            p_writer.start()
            for p_worker in p_workers:
                p_worker.start()

            for p_worker in p_workers:
                p_worker.join()

                # Gracefully terminate the program.
                if p_worker.exitcode != 0:
                    # print error message
                    self.logger.error("Error in worker process.")
                    raise Exception("Error in worker process.")


            # Stop the writer thread.
            q_writer.put(None)
            p_writer.join()

        except Exception:
            for p in p_workers:
                p.terminate()
            p_writer.terminate()
            raise

        # Parse the results
        arc_fna, arc_faa, bac_fna, bac_faa = dict(), dict(), dict(), dict()
        print("Length of q_results: ", q_results.qsize())
        q_results.put(None)
        cur_result = q_results.get()
        while cur_result is not None:
            gid, marker, gene_id, fna, faa, domain = cur_result

            if domain == 'd__Archaea':
                out_fna, out_faa = arc_fna, arc_faa
            elif domain == 'd__Bacteria':
                out_fna, out_faa = bac_fna, bac_faa
            else:
                raise Exception("Unknown domain.")

            if marker not in out_fna:
                out_fna[marker] = list()
            if marker not in out_faa:
                out_faa[marker] = list()

            out_fna[marker].append((gid, gene_id, fna))
            out_faa[marker].append((gid, gene_id, faa))

            cur_result = q_results.get()

        # Write the output to disk
        for out_dict, domain, ext in [(arc_fna, 'ar53', 'fna'),
                                      (arc_faa, 'ar53', 'faa'),
                                      (bac_fna, 'bac120', 'fna'),
                                      (bac_faa, 'bac120', 'faa')]:
            cur_root = os.path.join(self.output_dir, f'{domain}_{self.release_number}_individual_genes', ext)
            make_sure_path_exists(cur_root)

            for marker, seqs in sorted(out_dict.items()):
                cur_path = os.path.join(cur_root, f'{marker}.{ext}')

                with open(cur_path, 'w') as fh:
                    for cur_gid, cur_gene_id, cur_seq in seqs:
                        process_id = cur_gid
                        fh.write(f'>{process_id}\n{cur_seq}\n')
                pass

        self.logger.info('Done.')

    def protein_files(self, raw_taxonomy, genome_dirs_file, metadata_file):

        # Parse the metadata file
        self.logger.info('Parsing taxonomy file.')
        user_tax = parse_taxonomy_file(raw_taxonomy)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)

        domain = None
        make_sure_path_exists(os.path.join(self.output_dir, 'archaea'))
        make_sure_path_exists(os.path.join(self.output_dir, 'bacteria'))

        # Parse the genomic path
        self.logger.info(
            'Parsing genome directory files.')
        genome_dirs = parse_genomic_path_file(genome_dirs_file)

        # keep only the representative genomes in user_tax
        user_tax = {k: v for k, v in user_tax.items() if k in reps}


        for gid, taxonomy in tqdm(user_tax.items()):
            if gid in reps:
                path_of_interest = os.path.join(
                    genome_dirs.get(gid), 'prodigal')
                protein_files = glob.glob(f'{path_of_interest}/*_protein.faa.gz')
                for protein_file in protein_files:
                    if taxonomy.startswith('d__Archaea'):
                        domain = 'archaea'
                    else:
                        domain = 'bacteria'
                    shutil.copy(protein_file, os.path.join(
                        self.output_dir, domain, gid + '_protein.faa.gz'))

    def nucleotide_files(self, raw_taxonomy, genome_dirs_file, metadata_file):

        # Parse the metadata file
        self.logger.info('Parsing taxonomy file.')
        user_tax = parse_taxonomy_file(raw_taxonomy)

        # get representative genome for each GTDB species cluster
        self.logger.info('Parsing representative genomes from GTDB metadata.')
        reps = parse_rep_genomes(metadata_file)

        domain = None
        make_sure_path_exists(os.path.join(self.output_dir, 'archaea'))
        make_sure_path_exists(os.path.join(self.output_dir, 'bacteria'))

        # Parse the genomic path
        self.logger.info(
            'Parsing genome directory files ( including UBA paths).')
        genome_dirs = parse_genomic_path_file(genome_dirs_file)

        # keep only the representative genomes in user_tax
        user_tax = {k: v for k, v in user_tax.items() if k in reps}

        for gid, taxonomy in tqdm(user_tax.items()):
            if gid in reps:
                path_of_interest = os.path.join(
                    genome_dirs.get(gid), 'prodigal')
                protein_files = glob.glob(f'{path_of_interest}/*_protein.fna.gz')
                for protein_file in protein_files:
                    if taxonomy.startswith('d__Archaea'):
                        domain = 'archaea'
                    else:
                        domain = 'bacteria'
                    shutil.copy(protein_file, os.path.join(
                        self.output_dir, domain, gid + '_protein.fna.gz'))

    def metadata_files(self, metadata_file,
                       metadata_fields,
                       gtdb_sp_clusters_file):
        """Generate GTDB metadata files."""

        # get GTDB species clusters
        self.logger.info('Parsing all species clusters from GTDB metadata.')
        sp_clusters = parse_species_clusters(gtdb_sp_clusters_file)
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
        fout_ar = open(self.output_dir / f'ar53_metadata_r{self.release_number}.tsv', 'w',
                       encoding='utf-8')

        fields_sorted = sorted(fields)
        fields_sorted.remove('accession')
        fout_bac.write('{}\t{}\n'.format(
            'accession', '\t'.join(fields_sorted)))
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

                if gid in all_gids:
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

                    fout.write(gid)
                    for field in fields_sorted:
                        if field == 'gtdb_genome_representative':
                            gid_reps = line_split[field_index[field]]
                            fout.write('\t{}'.format(gid_reps))
                        else:
                            fout.write('\t{}'.format(
                                line_split[field_index[field]]))
                    fout.write('\n')

        fout_bac.close()
        fout_ar.close()

        self.logger.info(f'Wrote metadata for {count:,} genomes.')

    def copy_all_genes(self, bac120_taxonomy, ar53_taxonomy, user_gid_table, genome_dirs_file, extension):
        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)

        # parse genome directory files ( including UBA paths)
        self.logger.info(
            'Parsing genome directory files ( including UBA paths).')
        genome_dirs = parse_genomic_path_file(genome_dirs_file, user_gids)

        for domain, taxonomy_file in [('archaea', ar53_taxonomy), ('bacteria', bac120_taxonomy)]:
            count = 0
            with open(taxonomy_file, 'r') as tf:

                for line in tf:
                    gid = line.strip().split('\t')[0]
                    if gid in genome_dirs:
                        if '/GCA/' in genome_dirs.get(gid):
                            prefix = 'GCA'
                        else:
                            prefix = 'GCF'
                    else:
                        self.logger.info(
                            "{} is missing a directory".format(gid))
                        sys.exit()
                    out_dir = '{}/{}/{}'.format(self.output_dir, domain,
                                                '/'.join([prefix, gid[1:4], gid[4:7], gid[7:10]]))
                    #print([prefix, gid[1:4], gid[4:7], gid[7:10]])

                    make_sure_path_exists(out_dir)
                    full_id = '_'.join(os.path.basename(
                        genome_dirs.get(gid)).split('_')[:2])
                    #==========================================================
                    # #print(
                    #     '{}/prodigal/{}{}'.format(genome_dirs.get(gid), full_id, extension))
                    # #print('{}/{}{}'.format(out_dir, gid, extension))
                    #==========================================================
                    shutil.copy('{}/prodigal/{}{}'.format(genome_dirs.get(gid), full_id, extension),
                                '{}/{}{}'.format(out_dir, gid, extension))
                    count += 1
                    print(count, end='\r')

    def dict_file(self, taxonomy_file):
        """Generate GTDB dictionary file."""

        taxa = set()
        with open(taxonomy_file, encoding='utf-8') as f:
            for line in f:
                tokens = line.strip().split('\t')

                for taxon in tokens[1].split(';'):
                    taxa.add(taxon.strip()[3:])
                    taxa.add(taxon.strip())

                    if taxon.startswith('s__'):
                        generic, specific = taxon[3:].split()
                        taxa.add(generic.strip())
                        taxa.add(specific.strip())
                        taxa.add(f's__{generic.strip()}')
                        taxa.add(f'{generic[0]}. {specific}')
                        taxa.add(f's__{generic[0]}. {specific}')

        fout = open(self.output_dir / f'gtdb_r{self.release_number}.dic', 'w', encoding='utf-8')
        for t in sorted(taxa):
            fout.write(f'{t}\n')
        fout.close()

    def arb_files(self, metadata_file,
                  metadata_fields,
                  user_gid_table,
                  bac120_msa_file,
                  ar53_msa_file):
        """Generate ARB metadata files."""

        # parse user genome ID mapping table
        self.logger.info('Parsing user genome ID mapping table.')
        user_gids = parse_user_gid_table(user_gid_table)

        # read metadata fields to retain
        self.logger.info('Reading metadata fields to retain.')
        fields = set()
        with open(metadata_fields) as f:
            f.readline()
            for line in f:
                fields.add(line.strip().split('\t')[0])
        self.logger.info(f' ...identified {len(fields):,} fields to retain.')

        # remove fields that cause problems in ARB
        self.logger.info('Removing fields known to cause issues in ARB.')
        if 'gtdb_clustered_genomes' in fields:
            fields.remove('gtdb_clustered_genomes')

        # read genome metadata
        self.logger.info('Reading genome metadata.')
        metadata = defaultdict(lambda: {})
        with open(metadata_file, encoding='utf-8') as f:
            header = f.readline().strip().split('\t')

            for line in f:
                line_split = line.strip().split('\t')

                gid = line_split[0]
                gid = user_gids.get(gid, gid)

                for idx, field in enumerate(header):
                    if field not in fields:
                        continue

                    value = line_split[idx]
                    metadata[gid][field] = value

        # read MSA files
        seqs = {}
        for seq_file in [bac120_msa_file, ar53_msa_file]:
            for seq_id, seq in seq_io.read_seq(bac120_msa_file):
                seqs[seq_id] = seq

        # write out ARB metadata file for each domain
        for prefix, msa_file in [('bac120', bac120_msa_file), ('ar53', ar53_msa_file)]:
            # read MSA files
            seqs = seq_io.read(msa_file)

            # write out record for each genome
            fout = open(self.output_dir / f'{prefix}_arb_r{self.release_number}.tsv', 'w', encoding='utf-8')
            for gid, msa in seqs.items():
                gid = user_gids.get(gid, gid)

                fout.write("BEGIN\n")
                fout.write("db_name={}\n".format(gid))
                for field, value in metadata[gid].items():
                    # replace equal signs as these are incompatible with the
                    # ARB parser
                    if value:
                        value = value.replace('=', '/')

                    if len(value) > 1000:
                        print(gid, field, value)

                    fout.write("{}={}\n".format(field, value))

                fout.write("multiple_homologs=\n")
                fout.write("aligned_seq={}\n".format(msa))
                fout.write("END\n\n")
            fout.close()

    def validate(self,
                 taxonomy_file,
                 tree_file,
                 metadata_file,
                 msa_file,
                 ssu_file,
                 sp_clusters_file,
                 hq_genome_file):
        """Perform validation checks on GTDB website files."""

        # parse taxonomy file
        self.logger.info('Parsing taxonomy file.')
        taxonomy = {}
        for line in open(taxonomy_file):
            line_split = line.strip().split('\t')
            taxonomy[line_split[0]] = line_split[1]

        # get GTDB species clusters
        self.logger.info('Validating species clusters.')
        sp_clusters = parse_species_clusters(sp_clusters_file)
        all_gids = set()
        for gids in sp_clusters.values():
            all_gids.update(gids)
        assert(len(set(taxonomy).difference(all_gids)) == 0)

        # get HQ MIMAG genomes
        self.logger.info('Validating MIMAG HQ genome file.')
        hq_genomes = set()
        with open(hq_genome_file) as f:
            header = f.readline().strip().split('\t')

            gtdb_taxonomy_index = header.index('GTDB taxonomy')
            gtdb_rep_index = header.index('GTDB representative genome')

            for line in f:
                line_split = line.strip().split('\t')

                gid = line_split[0]
                if gid in taxonomy:
                    is_rep = (line_split[gtdb_rep_index].lower() == 'true')

                    assert(line_split[gtdb_taxonomy_index] == taxonomy[gid])
                    assert(gid in taxonomy)

                    if is_rep:
                        assert(gid in sp_clusters)

        # validate that tree spans representatives
        self.logger.info('Validating tree file.')
        tree = dendropy.Tree.get_from_path(tree_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        tree_gids = set([leaf.taxon.label for leaf in tree.leaf_node_iter()])
        assert(len(tree_gids.difference(sp_clusters)) == 0)

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
                    self.logger.info(
                        'GTDB species representative information is incorrect: {gid}')
                    sys.exit(-1)

                rep_id = line_split[gtdb_rep_genome_index]
                assert(taxonomy[rep_id] == line_split[gtdb_taxonomy_index])
                assert(rep_id in sp_clusters)
                assert(rep_id in tree_gids)

                if gid in taxonomy:
                    if taxonomy[gid] != line_split[gtdb_taxonomy_index]:
                        self.logger.error(f'GTDB taxonomy is incorrect in metadata file: {gid}')
                        sys.exit(-1)

        assert(len(set(taxonomy).symmetric_difference(metadata_gids)) == 0)

        # validate that MSA file spans representatives
        self.logger.info('Validating MSA file.')
        msa_seqs = seq_io.read(msa_file)
        assert(len(set(msa_seqs).symmetric_difference(tree_gids)) == 0)

        # validate that SSU file spans representatives
        self.logger.info('Validating SSU file.')
        ssu_seqs = set()
        for seq_id, seq, annotations in seq_io.read_seq(ssu_file, keep_annotation=True):
            ssu_taxonomy = annotations[0:annotations.find('[')].strip()
            if ssu_taxonomy != taxonomy[seq_id]:
                self.logger.error(f'GTDB taxonomy is incorrect in SSU file: {gid}')
            ssu_seqs.add(seq_id)

        assert(len(ssu_seqs.difference(tree_gids)) == 0)

        self.logger.info('Passed validation.')

    def tax_comp_files(self, old_taxonomy, new_taxonomy, changes_only, filter_placeholder_name, top_change):

        gtogfile = open(os.path.join(
            self.output_dir, 'genomes_summary.tsv'), "w")
        old_taxonomy_name = os.path.splitext(ntpath.basename(old_taxonomy))[0]
        new_taxonomy_name = os.path.splitext(ntpath.basename(new_taxonomy))[0]
        gtogfile.write("genome_id\t{}\t{}\tstatus\n".format(
            old_taxonomy_name, new_taxonomy_name))
        for dom in ["d__Archaea", "d__Bacteria"]:
            tempdict = {}  # temporary dictionary for the information in old taxonomy
            olddict = {}
            newdict = {}
            mapdict = {}  # dictionary of uba_mappings
            oldnumber = 0  # number of genes in old_tax
            newnumber = 0  # number of genes in new_tax
            intersection = 0  # number of genes common to both taxonomies

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
            self.rundistrib(olddict, newdict, "p", os.path.join(self.output_dir,
                                                                pref + "_phylum_difference.tsv"),
                            "Phylum", "Phyla",
                            old_taxonomy_name, new_taxonomy_name,
                            filter_placeholder_name, changes_only, top_change)
            print("- Analysis of classes")
            self.rundistrib(olddict, newdict, "c", os.path.join(self.output_dir,
                                                                pref + "_class_difference.tsv"),
                            "Class", "Classes",
                            old_taxonomy_name, new_taxonomy_name,
                            filter_placeholder_name, changes_only, top_change)
            print("- Analysis of orders")
            self.rundistrib(olddict, newdict, "o", os.path.join(self.output_dir,
                                                                pref + "_order_difference.tsv"),
                            "Order", "Orders",
                            old_taxonomy_name, new_taxonomy_name,
                            filter_placeholder_name, changes_only, top_change)
            print("- Analysis of families")
            self.rundistrib(olddict, newdict, "f", os.path.join(self.output_dir,
                                                                pref + "_family_difference.tsv"),
                            "Family", "Families",
                            old_taxonomy_name, new_taxonomy_name,
                            filter_placeholder_name, changes_only, top_change)
            print("- Analysis of genera")
            self.rundistrib(olddict, newdict, "g", os.path.join(self.output_dir,
                                                                pref + "_genus_difference.tsv"),
                            "Genus", "Genera",
                            old_taxonomy_name, new_taxonomy_name,
                            filter_placeholder_name, changes_only, top_change)

            print("- Analysis of species")
            self.rundistrib(olddict, newdict, "s", os.path.join(self.output_dir,
                                                                pref + "_species_difference.tsv"),
                            "Species", "Species",
                            old_taxonomy_name, new_taxonomy_name,
                            filter_placeholder_name, changes_only, top_change)

        gtogfile.close()

    def clean_empty_children(self, d):
        if "children" in d:
            v = d.get("children")
            if not v:
                del d["children"]
                return d
            else:
                tmp_list = []
                for it in v:
                    dict_tmp = self.clean_empty_children(it)
                    tmp_list.append(dict_tmp)
                d["children"] = tmp_list
                return d
        else:
            return d

    def count_children(self, d):
        if "extra_genomes" in d:
            d["countgen"] = d["extra_genomes"]
            return d
        if "children" in d and "value" not in d:
            v = d.get("children")
            countgen = 0
            for it in v:
                dict_tmp = self.count_children(it)
                countgen += dict_tmp["countgen"]
            d["countgen"] = countgen
            return d

        elif "children" not in d and "value" not in d:
            d["countgen"] = 1
            return d
        return d

    def add_type_species(self, d, type_spe_list):
        if "children" in d:
            v = d.get("children")
            for it in v:
                if it['name'] == 's__Coprobacter secundus':
                    a=1
                if it['name'] in type_spe_list:
                    it['type_material'] = 'type_species'
                else:
                    dict_tmp = self.add_type_species(it, type_spe_list)
        return d

    def conv(self, k, v, meta_information, type_spe_list):
        dict_rank = {"name": k,
                     "type": "d",
                     "children": list({"name": x,
                                       "type": "p",
                                       "children": list({"name": z,
                                                         "type": "c",
                                                         "children": list({"name": o,
                                                                           "type": "o",
                                                                           "children": list({"name": f,
                                                                                             "type": "f",
                                                                                             "children": list({"name": g,
                                                                                                               "type": "g",
                                                                                                               "children": list({"name": s,
                                                                                                                                 "type": "s",
                                                                                                                                 "children": list(self.process_genome(gen, meta_information) for gen in self.sorted_genome(v[x][z][o][f][g][s], meta_information))} for s in self.sort_ranks(v[x][z][o][f][g].keys()))} for g in self.sort_ranks(v[x][z][o][f].keys()))} for f in self.sort_ranks(v[x][z][o].keys()))} for o in self.sort_ranks(v[x][z].keys()))} for z in self.sort_ranks(v[x].keys()))} for x in self.sort_ranks(v.keys()))}
        cleaned_dict = self.clean_empty_children(dict_rank)
        cleaned_dict = self.count_children(cleaned_dict)
        cleaned_dict = self.add_type_species(cleaned_dict, type_spe_list)
        return cleaned_dict

    def sorted_genome(self, list_genomes, meta_information):
        print(list_genomes)
        for g in list_genomes:
            print(g)
            print(meta_information.get(g).get('rep'))
        list_representatives = [
            gen for gen in list_genomes if meta_information.get(gen).get('rep') == 't']
        list_type_species = [gen for gen in list_genomes if meta_information.get(gen).get(
            'type_str') and gen not in list_representatives]

        dict_quality = {gen: meta_information.get(gen).get('comp') - 5 * meta_information.get(gen).get(
            'conta') for gen in list_genomes if gen not in list_representatives + list_type_species}
        list_quality = sorted(dict_quality, key=dict_quality.get)
        full_list = list_representatives + list_type_species + list_quality
        if len(full_list) > 100:
            shorten_list = full_list[0:100]
            remaining_genomes = len(full_list) - len(shorten_list)
            shorten_list.append(f'+ {remaining_genomes:,} additional genomes.')
            return shorten_list
        return full_list

    def is_valid_name(self, rank_name):
        rank_name = rank_name[3:].split('_')[0]
        if rank_name.replace(' ', '').isalpha() and rank_name[0].isupper() and rank_name[1:].islower():
            return True
        return False

    def sort_ranks(self, list_ranks):
        results = []
        valid_names = sorted([x for x in list_ranks if self.is_valid_name(x)])
        placeholder_names = sorted(
            [x for x in list_ranks if x not in valid_names])
        results = valid_names + placeholder_names
        return results

    def process_genome(self, genome, meta_information):
        result_genome = {"name": genome, "type": 'genome'}
        if genome.startswith('+'):
            result_genome["extra_genomes"] = int(
                genome.split(' ')[1].replace(',', ''))
            return result_genome
        if meta_information.get(genome).get('type').startswith('type strain') and meta_information.get(genome).get('type_str') is True:
            result_genome["type_material"] = meta_information.get(
                genome).get('type').replace(' ', '_')
        if meta_information.get(genome).get('rep') == 't':
            result_genome["rep"] = True
        return result_genome

    def json_tree_parser(self, taxonomy_file, metadata_file):
        meta_information, type_spe_list = self.parse_metadata(metadata_file)
        tax_dict = {}
        with open(taxonomy_file, 'r') as f:
            for line in f:
                genome, tax_string = line.strip().split("\t")
                drank, prank, crank, orank, frank, grank, srank = tax_string.split(
                    ';')
                if len(drank) > 3:
                    if drank not in tax_dict:
                        tax_dict[drank] = {}
                else:
                    tax_dict.set_default(
                        'd__Undefined', {}).set_default(genome, {})
                    continue
                if len(prank) > 3:
                    if prank not in tax_dict[drank]:
                        tax_dict[drank][prank] = {}
                if len(crank) > 3:
                    if crank not in tax_dict[drank][prank]:
                        tax_dict[drank][prank][crank] = {}
                if len(orank) > 3:
                    if orank not in tax_dict[drank][prank][crank]:
                        tax_dict[drank][prank][crank][orank] = {}
                if len(frank) > 3:
                    if frank not in tax_dict[drank][prank][crank][orank]:
                        tax_dict[drank][prank][crank][orank][frank] = {}
                if len(grank) > 3:
                    if grank not in tax_dict[drank][prank][crank][orank][frank]:
                        tax_dict[drank][prank][crank][orank][frank][grank] = {}
                if len(srank) > 3:
                    if srank not in tax_dict[drank][prank][crank][orank][frank][grank]:
                        tax_dict[drank][prank][crank][orank][frank][grank][srank] = []
                    tax_dict[drank][prank][crank][orank][frank][grank][srank].append(
                        genome)
        with open('genome_taxonomy_r{}_count.json'.format(self.release_number), 'w') as outfile:
            json.dump({"name": "Domain", "type": "root", "children": list(self.conv(
                k, tax_dict[k], meta_information, type_spe_list) for k in sorted(tax_dict.keys()))}, outfile, check_circular=False)

    def parse_metadata(self, metadatafile):
        result = {}
        type_spe_list = []
        with open(metadatafile) as metafile:
            header = metafile.readline()
            listheaders = header.rstrip().split('\t')
            gr = listheaders.index('gtdb_representative')
            gtd = listheaders.index("gtdb_type_designation_ncbi_taxa")
            gcomp = listheaders.index("checkm_completeness")
            gconta = listheaders.index("checkm_contamination")
            gts = listheaders.index("gtdb_species_type_strain")
            # We based gtsofg on gtdb_genus_type_species (based in GTDB) instead of
            # gtdb_type_species_of_genus ( based on NCBI)
            # gtsofg = listheaders.index("gtdb_type_species_of_genus")
            gtsofg = listheaders.index("gtdb_genus_type_species")
            gtdb_tax = listheaders.index("gtdb_taxonomy")
            for line in metafile:
                infos = line.rstrip().split('\t')
                result[infos[0]] = {'rep': infos[gr],
                                    'type': infos[gtd],
                                    'comp': float(infos[gcomp]),
                                    'conta': float(infos[gconta]),
                                    'type_spe': infos[gtsofg],
                                    'type_str': self.str2bool(infos[gts])}
                if infos[gtsofg] == 't':
                    type_spe_list.append(infos[gtdb_tax].split(';')[6])

        return (result, type_spe_list)

    def str2bool(self, v):
        return v.lower() in ("yes", "true", "t", "1")

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

    def rundistrib(self,
                   olddict, newdict,
                   rank,
                   outfile,
                   rank_label, plural_rank_label,
                   old_taxonomy_name, new_taxonomy_name,
                   filter_placeholder_name,
                   changes_only,
                   top_change):
        """Compare taxa in a pair of taxonomies."""

        orderrank = ["d", "p", "c", 'o', 'f', 'g', 's']

        # take the different ranks in the dictionaries and then count the appearances of each pair
        # of ranks ie. (p__Fusobacteria, p_Fusobacteriota),200
        different_ranks = []
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

            if changes_only and k == results_dict[k]['genomes'][0][0] and results_dict[k]['nber_g'] == results_dict[k]['genomes'][0][1]:
                to_pop.append(k)
            if filter_placeholder_name and self.hasNumber(k) == True and results_dict[k]['genomes'][0][0] == rank + "__":
                to_pop.append(k)
            elif filter_placeholder_name and self.hasNumber(results_dict[k]['genomes'][0][0]) == True and k == rank + "__":
                to_pop.append(k)

        self.logger.info(
            'Removed {} taxa from tables based on filtering parameters.'.format(len(to_pop)))

        for k in to_pop:
            results_dict.pop(k)

        outf = open(outfile, "w")
        outf.write("{0} {1}\tNo. genomes\tNo. {2} {3}\t{2} {3}\n".format(
            old_taxonomy_name,
            rank_label,
            new_taxonomy_name,
            plural_rank_label))
        for k, v in results_dict.items():
            for sk, sv in olddict.items():
                if sv.get(rank) == k:
                    prankold = sv.get(orderrank[orderrank.index(rank) - 1])
                    break
            number_sub = len(v.get("genomes"))

            results = []
            for newg in v.get("genomes"):
                newg_name = newg[0]
                newg_numb = float(newg[1])
                for sk, sv in newdict.items():
                    if sv.get(rank) == newg_name:
                        pranknew = sv.get(orderrank[orderrank.index(rank) - 1])
                        break

                res = newg_name
                if not newg_name.startswith('s__') or newg_name == 's__':
                    # report parent taxa, except at the rank of species unless the species
                    # name (and hence genus) is unspecified
                    if prankold.replace('Candidatus ', '') != pranknew.replace('Candidatus ', ''):
                        res = res + "(" + pranknew + ")"
                res = "{0} {1}%".format(res, round(
                    (newg_numb / v.get("nber_g")) * 100, 2))
                results.append(res)

            outf.write("{0}\t{1}\t{2}\t{3}\n".format(
                k, v.get("nber_g"), number_sub, ", ".join(results)))
        outf.close()
