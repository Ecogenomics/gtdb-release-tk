from gtdb_release_tk.common import assert_file_exists, optional_int, optional_float

from typing import Dict

from gtdb_release_tk.models.taxonomy_string import TaxonomyString


class MetadataFile(object):

    def __init__(self, path, header, rows):
        self.path: str = path
        self.header: Dict[str, int] = header
        self.rows: Dict[str, MetadataRow] = rows

    @classmethod
    def read(cls, path):
        """Instantiates a new MetadataFile from disk."""
        assert_file_exists(path)
        rows = dict()
        headers = dict()
        with open(path, 'r') as f:
            f_headers = f.readline().rstrip()
            for i, f_header in enumerate(f_headers.split()):
                headers[f_header] = i
            for line in f.readlines():
                cur_row = MetadataRow(headers, line.rstrip().split('\t'))
                rows[cur_row.accession] = cur_row
        return cls(path, headers, rows)

    def get_rep_tax(self):
        out = dict()
        for gid, row in self.rows.items():
            if row.gtdb_representative.startswith('t'):
                out[gid] = row.gtdb_taxonomy
        return out

class MetadataRow(object):

    def __init__(self, header: Dict[str, int], row):
        self.accession = row[header['accession']]
        self.ambiguous_bases = row[header['ambiguous_bases']]
        self.checkm_completeness = optional_float(row[header['checkm_completeness']])
        self.checkm_contamination = optional_float(row[header['checkm_contamination']])
        self.checkm_marker_count = row[header['checkm_marker_count']]
        self.checkm_marker_lineage = row[header['checkm_marker_lineage']]
        self.checkm_marker_set_count = row[header['checkm_marker_set_count']]
        self.checkm_strain_heterogeneity = row[header['checkm_strain_heterogeneity']]
        self.coding_bases = row[header['coding_bases']]
        self.coding_density = row[header['coding_density']]
        self.contig_count = optional_int(row[header['contig_count']])
        self.gc_count = optional_int(row[header['gc_count']])
        self.gc_percentage = optional_float(row[header['gc_percentage']])
        self.genome_size = optional_int(row[header['genome_size']])
        self.gtdb_genome_representative = row[header['gtdb_genome_representative']]
        self.gtdb_representative = row[header['gtdb_representative']]
        self.gtdb_taxonomy = TaxonomyString(row[header['gtdb_taxonomy']])
        self.gtdb_type_designation = row[header['gtdb_type_designation']]
        self.gtdb_type_designation_sources = row[header['gtdb_type_designation_sources']]
        self.gtdb_type_species_of_genus = row[header['gtdb_type_species_of_genus']]
        self.l50_contigs = row[header['l50_contigs']]
        self.l50_scaffolds = optional_int(row[header['l50_scaffolds']])
        self.longest_contig = row[header['longest_contig']]
        self.longest_scaffold = optional_int(row[header['longest_scaffold']])
        self.lsu_23s_contig_len = row[header['lsu_23s_contig_len']]
        self.lsu_23s_count = row[header['lsu_23s_count']]
        self.lsu_23s_length = row[header['lsu_23s_length']]
        self.lsu_23s_query_id = row[header['lsu_23s_query_id']]
        self.lsu_5s_contig_len = row[header['lsu_5s_contig_len']]
        self.lsu_5s_count = row[header['lsu_5s_count']]
        self.lsu_5s_length = row[header['lsu_5s_length']]
        self.lsu_5s_query_id = row[header['lsu_5s_query_id']]
        self.lsu_silva_23s_blast_align_len = row[header['lsu_silva_23s_blast_align_len']]
        self.lsu_silva_23s_blast_bitscore = row[header['lsu_silva_23s_blast_bitscore']]
        self.lsu_silva_23s_blast_evalue = row[header['lsu_silva_23s_blast_evalue']]
        self.lsu_silva_23s_blast_perc_identity = row[header['lsu_silva_23s_blast_perc_identity']]
        self.lsu_silva_23s_blast_subject_id = row[header['lsu_silva_23s_blast_subject_id']]
        self.lsu_silva_23s_taxonomy = row[header['lsu_silva_23s_taxonomy']]
        self.mean_contig_length = row[header['mean_contig_length']]
        self.mean_scaffold_length = row[header['mean_scaffold_length']]
        self.mimag_high_quality = row[header['mimag_high_quality']]
        self.mimag_low_quality = row[header['mimag_low_quality']]
        self.mimag_medium_quality = row[header['mimag_medium_quality']]
        self.n50_contigs = optional_int(row[header['n50_contigs']])
        self.n50_scaffolds = optional_int(row[header['n50_scaffolds']])
        self.ncbi_assembly_level = row[header['ncbi_assembly_level']]
        self.ncbi_assembly_name = row[header['ncbi_assembly_name']]
        self.ncbi_assembly_type = row[header['ncbi_assembly_type']]
        self.ncbi_bioproject = row[header['ncbi_bioproject']]
        self.ncbi_biosample = row[header['ncbi_biosample']]
        self.ncbi_contig_count = row[header['ncbi_contig_count']]
        self.ncbi_contig_n50 = row[header['ncbi_contig_n50']]
        self.ncbi_country = row[header['ncbi_country']]
        self.ncbi_date = row[header['ncbi_date']]
        self.ncbi_genbank_assembly_accession = row[header['ncbi_genbank_assembly_accession']]
        self.ncbi_genome_category = row[header['ncbi_genome_category']]
        self.ncbi_genome_representation = row[header['ncbi_genome_representation']]
        self.ncbi_isolate = row[header['ncbi_isolate']]
        self.ncbi_isolation_source = row[header['ncbi_isolation_source']]
        self.ncbi_lat_lon = row[header['ncbi_lat_lon']]
        self.ncbi_molecule_count = row[header['ncbi_molecule_count']]
        self.ncbi_ncrna_count = row[header['ncbi_ncrna_count']]
        self.ncbi_organism_name = row[header['ncbi_organism_name']]
        self.ncbi_protein_count = row[header['ncbi_protein_count']]
        self.ncbi_refseq_category = row[header['ncbi_refseq_category']]
        self.ncbi_rrna_count = row[header['ncbi_rrna_count']]
        self.ncbi_scaffold_count = row[header['ncbi_scaffold_count']]
        self.ncbi_scaffold_l50 = row[header['ncbi_scaffold_l50']]
        self.ncbi_scaffold_n50 = row[header['ncbi_scaffold_n50']]
        self.ncbi_scaffold_n75 = row[header['ncbi_scaffold_n75']]
        self.ncbi_scaffold_n90 = row[header['ncbi_scaffold_n90']]
        self.ncbi_seq_rel_date = row[header['ncbi_seq_rel_date']]
        self.ncbi_spanned_gaps = row[header['ncbi_spanned_gaps']]
        self.ncbi_species_taxid = row[header['ncbi_species_taxid']]
        self.ncbi_ssu_count = row[header['ncbi_ssu_count']]
        self.ncbi_strain_identifiers = row[header['ncbi_strain_identifiers']]
        self.ncbi_submitter = row[header['ncbi_submitter']]
        self.ncbi_taxid = row[header['ncbi_taxid']]
        self.ncbi_taxonomy = row[header['ncbi_taxonomy']]
        self.ncbi_taxonomy_unfiltered = row[header['ncbi_taxonomy_unfiltered']]
        self.ncbi_total_gap_length = row[header['ncbi_total_gap_length']]
        self.ncbi_total_length = optional_int(row[header['ncbi_total_length']])
        self.ncbi_translation_table = row[header['ncbi_translation_table']]
        self.ncbi_trna_count = row[header['ncbi_trna_count']]
        self.ncbi_type_material_designation = row[header['ncbi_type_material_designation']]
        self.ncbi_ungapped_length = row[header['ncbi_ungapped_length']]
        self.ncbi_unspanned_gaps = row[header['ncbi_unspanned_gaps']]
        self.ncbi_wgs_master = row[header['ncbi_wgs_master']]
        self.protein_count = row[header['protein_count']]
        self.scaffold_count = optional_int(row[header['scaffold_count']])
        self.ssu_contig_len = row[header['ssu_contig_len']]
        self.ssu_count = row[header['ssu_count']]
        self.ssu_gg_blast_align_len = row[header['ssu_gg_blast_align_len']]
        self.ssu_gg_blast_bitscore = row[header['ssu_gg_blast_bitscore']]
        self.ssu_gg_blast_evalue = row[header['ssu_gg_blast_evalue']]
        self.ssu_gg_blast_perc_identity = row[header['ssu_gg_blast_perc_identity']]
        self.ssu_gg_blast_subject_id = row[header['ssu_gg_blast_subject_id']]
        self.ssu_gg_taxonomy = row[header['ssu_gg_taxonomy']]
        self.ssu_length = row[header['ssu_length']]
        self.ssu_query_id = row[header['ssu_query_id']]
        self.ssu_silva_blast_align_len = row[header['ssu_silva_blast_align_len']]
        self.ssu_silva_blast_bitscore = row[header['ssu_silva_blast_bitscore']]
        self.ssu_silva_blast_evalue = row[header['ssu_silva_blast_evalue']]
        self.ssu_silva_blast_perc_identity = row[header['ssu_silva_blast_perc_identity']]
        self.ssu_silva_blast_subject_id = row[header['ssu_silva_blast_subject_id']]
        self.ssu_silva_taxonomy = row[header['ssu_silva_taxonomy']]
        self.total_gap_length = row[header['total_gap_length']]
        self.trna_aa_count = row[header['trna_aa_count']]
        self.trna_count = row[header['trna_count']]
        self.trna_selenocysteine_count = row[header['trna_selenocysteine_count']]

    def __repr__(self):
        return f'{self.accession}'
