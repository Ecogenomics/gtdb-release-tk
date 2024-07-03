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
__copyright__ = 'Copyright 2019'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import logging
import operator
import gzip
from pathlib import PurePath
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy
from biolib.seq_io import read_seq

   
class RepsPerRank(object):
    """Select representative genomes at each taxonomic rank."""

    def __init__(self, release_number, output_dir):
        
        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def read_ssu_file(self, ssu_fasta_file):
        """Read length of SSU sequences for genomes."""
        
        ssu_length = {}
        for seq_id, seq in read_seq(ssu_fasta_file):
            gid = seq_id.split('~')[0]
            
            if gid in ssu_length and len(seq) < ssu_length[gid]:
                continue
                
            ssu_length[gid] = len(seq) - seq.upper().count('N')
            
        return ssu_length
        
    def read_msa_file(self, msa_file):
        """Determine percentage of amino acids for each genome in MSA file."""
        
        msa_perc = {}
        for seq_id, seq in read_seq(msa_file):
            seq = seq.upper()
            aa = len(seq) - seq.count('-') - seq.count('_') - seq.count('*')
            msa_perc[seq_id] = aa*100.0 / len(seq)
        
        return msa_perc
        
    def read_metadata(self, bac120_metadata_file, ar120_metadata_file):
        """Read metadata for genomes."""
        
        Metadata = namedtuple('Metadata', ['gtdb_representative',
                                            'gtdb_type_strain',
                                            'genome_quality',
                                            'genome_comp',
                                            'genome_cont',
                                            'gtdb_taxonomy',
                                            'ncbi_taxonomy'])

        metadata = {}
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            maybe_gzip_open = open
            if mf.endswith('.gz'):
                maybe_gzip_open = gzip.open

            with maybe_gzip_open(mf, mode = 'rt', encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                genome_index = header.index('accession')
                comp_index = header.index('checkm_completeness')
                cont_index = header.index('checkm_contamination')
                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                ncbi_taxonomy_index = header.index('ncbi_taxonomy')
                gtdb_rep_index = header.index('gtdb_representative')
                gtdb_type_index = header.index('gtdb_type_designation_ncbi_taxa')
                    
                for line in f:
                    line_split = line.strip().split('\t')
                    
                    gid = line_split[genome_index]

                    comp = float(line_split[comp_index])
                    cont = float(line_split[cont_index])
                    quality = comp-5*cont
                        
                    ncbi_taxonomy = line_split[ncbi_taxonomy_index]
                    gtdb_taxonomy = [x.strip() for x in line_split[gtdb_taxonomy_index].split(';')]

                    gtdb_representative = False
                    if line_split[gtdb_rep_index] == 't':
                        gtdb_representative = True
                       
                    gtdb_type_strain = False
                    if line_split[gtdb_type_index] == 'type strain of species':
                        gtdb_type_strain = True

                    metadata[gid] = Metadata(gtdb_representative=gtdb_representative,
                                                    gtdb_type_strain=gtdb_type_strain,
                                                    genome_quality=quality,
                                                    genome_comp=comp,
                                                    genome_cont=cont,
                                                    gtdb_taxonomy=gtdb_taxonomy,
                                                    ncbi_taxonomy=ncbi_taxonomy)

        return metadata
   
    def get_genomes_in_rank(self, domain, metadata):
        """Get genomes in each taxa at each rank."""
        
        genomes_in_rank = defaultdict(lambda: defaultdict(set))
        for gid, m in metadata.items():
            taxa = m.gtdb_taxonomy
            
            if ((domain == 'archaea' and taxa[0] == 'd__Archaea')
                or (domain == 'bacteria' and taxa[0] == 'd__Bacteria')):

                for rank_index, taxon in enumerate(taxa):
                    genomes_in_rank[rank_index][taxon].add(gid)
        
        return genomes_in_rank
        
    def filter_genomes(self, 
                        genome_ids,
                        metadata,
                        ar122_msa_per, 
                        bac120_msa_per, 
                        min_msa_perc):
        """Filter genomes based on selection criteria."""
        
        quality_genome_ids = {}
        filtered_desc = {}
        for gid in genome_ids:
            qual = metadata[gid].genome_quality
            
            if max(ar122_msa_per.get(gid, 0), bac120_msa_per.get(gid, 0)) < min_msa_perc:
                filtered_desc[gid] = 'Failed MSA test: %.1f' % max(ar122_msa_per.get(gid, 0), bac120_msa_per.get(gid, 0))
                continue
                 
            quality_genome_ids[gid] = qual
            
        return quality_genome_ids, filtered_desc
        
    def ssu_test(self, 
                    genome_ids, 
                    ssu_length, 
                    min_ssu_len):
        """Filter genomes based on selection criteria."""
        
        passed_ssu_test = set()
        for gid in genome_ids:
            if ssu_length.get(gid, -1) < min_ssu_len:
                continue

            passed_ssu_test.add(gid)
            
        return passed_ssu_test
        
    def qualitative_metadata(self, genome_ids, metadata):
        """Get qualitative annotation information about genomes."""
        
        gtdb_rep = set()
        gtdb_type_strain = set()
        for gid in genome_ids:
            m = metadata[gid]
            
            if m.gtdb_representative:
                gtdb_rep.add(gid)
                
            if m.gtdb_type_strain:
                gtdb_type_strain.add(gid)
                
        return gtdb_rep, gtdb_type_strain
        
    def sample_genomes(self,
                        genomes_per_taxon,
                        domain, 
                        rank_to_sample, 
                        genomes_in_rank, 
                        metadata,
                        ssu_length, min_ssu_len,
                        ar122_msa_per, bac120_msa_per, min_msa_perc):
        """Sample genomes from taxa."""
        
        rank_label = Taxonomy.rank_labels[rank_to_sample]

        # sample genomes
        output_file = os.path.join(self.output_dir, 'gtdb_r{}_selected_genomes.{}.{}.tsv'.format(
                                                        self.release_number,
                                                        domain, 
                                                        rank_label))
                                                        
        fout = open(output_file, 'w')
        fout.write('#Genome ID\t{}\tGTDB taxonomy\tNCBI taxonomy\tCompleteness\tContamination\tQuality\tDescription\n'.format(
                        rank_label.capitalize()))
        
        missing_file = os.path.join(self.output_dir, 'gtdb_r{}_unrepresented_taxa.{}.{}.tsv'.format(
                                    self.release_number,
                                    domain, 
                                    rank_label))
        fout_missing = open(missing_file, 'w')
        
        no_taxon_rep = set()
        ssu_status = {}
        for taxon, genome_ids in genomes_in_rank[rank_to_sample].items():
            # filter genomes based on general selection criteria
            quality_genome_ids, filtered_desc = self.filter_genomes(genome_ids,
                                                                    metadata,
                                                                    ar122_msa_per, 
                                                                    bac120_msa_per, 
                                                                    min_msa_perc)
                                                        
            # filter genomes based on 16S rRNA
            pass_ssu_test = self.ssu_test(quality_genome_ids.keys(),
                                            ssu_length, 
                                            min_ssu_len)
                                            
            if pass_ssu_test:
                # restrict genomes to those that passed filtering and have a SSU sequences
                failed_ssu_test = set(quality_genome_ids).difference(pass_ssu_test)
                for genome_id in failed_ssu_test:
                    quality_genome_ids.pop(genome_id)
                                  
                ssu_status[taxon] = True
            else:
                # based selection of genomes passing general filtering criteria and accept
                # that this taxon doesn't have a quality genome with a SSU sequence
                if len(quality_genome_ids): # otherwise the failure is due to other filtering
                    ssu_status[taxon] = False 
                                                        
            # get qualitative information about genomes
            gtdb_rep, gtdb_type_strain = self.qualitative_metadata(quality_genome_ids,
                                                                metadata)
            
            gtdb_type_strain_ids = gtdb_type_strain.intersection(quality_genome_ids)
            gtdb_rep_ids = gtdb_rep.intersection(quality_genome_ids)

            sampled_genomes = []
            if gtdb_type_strain_ids:
                quality = {gid:quality_genome_ids[gid] for gid in gtdb_type_strain_ids}
                quality_sorted = sorted(quality.items(), key=operator.itemgetter(1,0), reverse=True)
                sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])

            if gtdb_rep_ids and len(sampled_genomes) != genomes_per_taxon:
                quality = {gid:quality_genome_ids[gid] for gid in gtdb_rep_ids if gid not in sampled_genomes}
                quality_sorted = sorted(quality.items(), key=operator.itemgetter(1,0), reverse=True)
                sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])

            if len(sampled_genomes) != genomes_per_taxon:
                remaining_genomes = set(quality_genome_ids.keys()).difference(sampled_genomes)

                if remaining_genomes:
                    quality = {gid:quality_genome_ids[gid] for gid in remaining_genomes}
                    quality_sorted = sorted(quality.items(), key=operator.itemgetter(1,0), reverse=True)
                    sampled_genomes.extend([x[0] for x in quality_sorted[0:genomes_per_taxon - len(sampled_genomes)]])
                
            for gid in sampled_genomes:
                desc = 'Selected based on genome quality'
                if gid in gtdb_type_strain_ids:
                    desc = 'Type strain of species'
                elif gid in gtdb_rep_ids:
                    desc = 'GTDB representative genome'
                
                fout.write('{}\t{}\t{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}\t{}\n'.format(
                                gid,
                                taxon, 
                                ';'.join(metadata[gid].gtdb_taxonomy),
                                metadata[gid].ncbi_taxonomy, 
                                metadata[gid].genome_comp, 
                                metadata[gid].genome_cont,
                                metadata[gid].genome_quality,
                                desc))
                                                                            
            if not sampled_genomes:
                fout_missing.write('{}\n'.format(taxon))
                for gid in genome_ids:
                    fout_missing.write('\t{}\t{}\n'.format(
                                        gid, 
                                        filtered_desc[gid]))
                no_taxon_rep.add(taxon)
                
        fout.close()
        fout_missing.close()
        
        if no_taxon_rep:
            self.logger.info('  {}: no reps for {:,} taxa'.format(rank_label, len(no_taxon_rep)))
            
        self.logger.info('  {}: {:,} of {:,} taxa are missing a suitable 16S rRNA gene'.format(
                            rank_label, 
                            sum([1 for status in ssu_status.values() if not status]), 
                            len(ssu_status)))
        
    def run(self, 
                bac120_metadata_file,
                ar53_metadata_file,
                bac120_msa_file,
                ar53_msa_file,
                ssu_fasta_file,
                genomes_per_taxon,
                min_ssu_len,
                min_msa_perc):
        """Select representative genomes at each taxonomic rank."""
                
        # get genome metadata
        self.logger.info('Reading GTDB metadata.')
        metadata = self.read_metadata(bac120_metadata_file, ar53_metadata_file)
        
        # get percentage of amino acids in MSAs
        self.logger.info('Determining percentage of amino acids in MSAs.')
        ar122_msa_per = self.read_msa_file(ar53_msa_file)
        bac120_msa_per = self.read_msa_file(bac120_msa_file)

        # get length of SSU sequences
        self.logger.info('Reading SSU length.')
        ssu_length = self.read_ssu_file(ssu_fasta_file)
            
        # sample genomes per taxonomic group at the ranks of phylum to genus
        self.logger.info('Selecting representatives.')
        for domain in ['archaea', 'bacteria']:
            self.logger.info('Processing {} taxa:'.format(domain.capitalize()))
            
            # get genomes for each taxon
            genomes_in_rank = self.get_genomes_in_rank(domain, metadata)

            # dereplicate from each taxon at each rank
            for rank_to_sample in range(1, 7):
                self.sample_genomes(genomes_per_taxon,
                                        domain, 
                                        rank_to_sample, 
                                        genomes_in_rank, 
                                        metadata,
                                        ssu_length, min_ssu_len,
                                        ar122_msa_per, bac120_msa_per, min_msa_perc)
