
import os
import sys


def parse_user_gid_table(user_gid_table):
    """Parse user genome ID table."""
    
    user_gids = {}
    with open(user_gid_table) as f:
        for line in f:
            line_split = line.strip().split('\t')
            
            if len(line_split) == 3:
                user_gids[line_split[0]] = line_split[2]
            elif len(line_split) == 2:
                user_gids[line_split[0]] = line_split[1]
            else:
                print('Error parsing user genome ID table.')
                sys.exit(-1)
    
    return user_gids
    
    
def parse_rep_genomes(gtdb_metadata_file, user_gids):
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
                
            gid = user_gids.get(gid, gid)
            gtdb_taxa = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')]
            reps[gid] = gtdb_taxa
            
    return reps
    
    
def parse_genomic_path_file(genome_path_file, user_gids):
    """Parse path to data directory for each genome."""
    
    genome_paths = {}
    with open(genome_path_file) as f:
        for line in f:
            line_split = line.strip().split('\t')
            
            gid = line_split[0]
            gid = user_gids.get(gid, gid)
            
            genome_paths[gid] = line_split[1]
            
    return genome_paths
    
    
def parse_species_clusters(gtdb_sp_clusters_file, user_gids):
    """Parse GTDB species clusters."""
    
    sp_clusters = {}
    with open(gtdb_sp_clusters_file) as f:
        header = f.readline().strip().split('\t')
        
        type_genome_index = header.index('Type genome')
        cluster_size_index = header.index('No. clustered genomes')
        clustered_genomes_index = header.index('Clustered genomes')
        
        for line in f:
            line_split = line.strip().split('\t')
            
            rid = line_split[type_genome_index]
            rid = user_gids.get(rid, rid)
            sp_clusters[rid] = set([rid])
            
            cluster_size = int(line_split[cluster_size_index])
            if cluster_size > 0:
                for cid in line_split[clustered_genomes_index].split(','):
                    cid = user_gids.get(cid, cid)
                    sp_clusters[rid].add(cid)
            
    return sp_clusters