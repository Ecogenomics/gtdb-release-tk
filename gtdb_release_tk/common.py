import gzip
import os
import sys
from collections import defaultdict, namedtuple
from biolib.common import canonical_gid


ENV_CATEGORIES = set(['derived from single cell',
                      'derived from metagenome',
                      'derived from environmental_sample',
                      'derived from environmental sample'])


def canonical_taxon_name(taxon):
    """Get canonical version of taxon."""

    if '_' in taxon[3:]:
        taxon = taxon[0:taxon.rfind('_')]

    return taxon


def sp_cluster_type_category(gids, genome_category):
    """Determine genome types in each species cluster."""

    sp_genome_types = set()
    for gid in gids:
        if gid.startswith('UBA') or genome_category[gid] in ENV_CATEGORIES:
            sp_genome_types.add('ENV')
        elif genome_category[gid] == 'none':
            sp_genome_types.add('ISOLATE')
        else:
            print(f'Unrecognized genome category: {genome_category[gid]}')
            sys.exit(-1)

    if len(sp_genome_types) == 2:
        category = 'BOTH'
    elif list(sp_genome_types)[0] == 'ENV':
        category = 'ENV'
    elif list(sp_genome_types)[0] == 'ISOLATE':
        category = 'ISOLATE'
    else:
        print(f'Error in sp_cluster_type_category: {sp_genome_types}')
        sys.exit(-1)

    return category


def parse_user_gid_table(user_gid_table):
    """Parse user genome ID table."""

    user_gids = {}
    with open(user_gid_table) as f:
        for line in f:
            line_split = line.strip().split('\t')

            if len(line_split) == 3:
                user_gids[line_split[0]] = line_split[2]
                # Temporary fix for migration from UBA to GCA R95
                user_gids[line_split[1]] = line_split[2]
            elif len(line_split) == 2:
                user_gids[line_split[0]] = line_split[1]
            else:
                print('Error parsing user genome ID table.')
                sys.exit(-1)

    return user_gids


def parse_canonical_gid_table(canonical_gid_table):
    """Parse user canonical ID table."""

    canonical_gids = {}
    with open(canonical_gid_table) as f:
        for line in f:
            canonical_id, ncbi_id = line.strip().split('\t')
            canonical_gids[canonical_id] = ncbi_id
    return canonical_gids


def parse_rep_genomes(gtdb_metadata_file):
    """Parse species representative genomes and their GTDB taxonomy."""

    reps = {}
    with open(gtdb_metadata_file, encoding='utf-8') as f:
        header = f.readline().strip().split('\t')

        gtdb_taxonomy_index = header.index('gtdb_taxonomy')
        gtdb_rep_index = header.index('gtdb_representative')

        for line in f:
            line_split = line.strip().split('\t')

            gid = line_split[0]
            gtdb_rep = line_split[gtdb_rep_index]
            if gtdb_rep.lower().startswith('f'):
                continue

            gtdb_taxa = [t.strip()
                         for t in line_split[gtdb_taxonomy_index].split(';')]
            reps[gid] = gtdb_taxa

    return reps


def parse_genomic_path_file(genome_path_file):
    """Parse path to data directory for each genome."""

    genome_paths = {}
    with open(genome_path_file) as f:
        for line in f:
            line_split = line.strip().split('\t')

            gid = line_split[0]
            genome_paths[gid] = line_split[1]

    return genome_paths


def parse_species_clusters(gtdb_sp_clusters_file):
    """Parse GTDB species clusters."""

    sp_clusters = {}
    with open(gtdb_sp_clusters_file) as f:
        header = f.readline().strip().split('\t')

        if 'Type genome' in header:
            type_genome_index = header.index('Type genome')
        else:
            type_genome_index = header.index('Representative genome')
        cluster_size_index = header.index('No. clustered genomes')
        clustered_genomes_index = header.index('Clustered genomes')

        for line in f:
            line_split = line.strip().split('\t')

            rid = line_split[type_genome_index]
            sp_clusters[rid] = set([rid])

            cluster_size = int(line_split[cluster_size_index])
            if cluster_size > 0:
                for cid in line_split[clustered_genomes_index].split(','):
                    sp_clusters[rid].add(cid)

    return sp_clusters


def parse_gtdb_metadata(metadata_file, fields):
    """Parse genome quality from GTDB metadata.

    Parameters
    ----------
    metadata_file : str
        Metadata for all genomes in CSV file.
    fields : iterable
        Fields  to read.

    Return
    ------
    dict : d[genome_id] -> namedtuple
        Value for fields indicted by genome IDs.
    """

    gtdb_metadata = namedtuple('gtdb_metadata', ' '.join(fields))
    m = {}

    with open(metadata_file, encoding='utf-8') as f:
        headers = f.readline().strip().split('\t')

        genome_index = headers.index('accession')

        indices = []
        for field in fields:
            indices.append(headers.index(field))

        for line in f:
            line_split = line.strip().split('\t')
            gid = line_split[genome_index]

            values = []
            for i in indices:
                # save values as floats or strings
                v = line_split[i]
                try:
                    values.append(float(v))
                except ValueError:
                    if v is None or v == '' or v == 'none':
                        values.append(None)
                    elif v == 'f' or v.lower() == 'false':
                        values.append(False)
                    elif v == 't' or v.lower() == 'true':
                        values.append(True)
                    else:
                        values.append(v)
            m[gid] = gtdb_metadata._make(values)

    return m


def parse_tophit_file(path):
    out = dict()

    open_file = open
    if path.endswith('.gz'):
        open_file = gzip.open

    with open_file(path, 'rt') as fh:
        fh.readline()
        for line in fh.readlines():
            gid, hits = line.strip().split('\t')
            for hit in hits.split(';'):
                fam_id, e_val, bitscore = hit.split(',')
                e_val, bitscore = float(e_val), float(bitscore)
                if fam_id not in out:
                    out[fam_id] = list()
                out[fam_id].append((gid, e_val, bitscore))

    return out


def parse_taxonomy_file(path):
    out = dict()
    with open(path, 'r') as fh:
        for line in fh.readlines():
            gid, tax = line.strip().split('\t')
            out[gid] = tax
            #out[gid.replace('RS_', '').replace('GB_', '')] = tax
    return out
