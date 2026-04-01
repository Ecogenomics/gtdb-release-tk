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

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2025'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'

import dendropy
import pandas as pd
import os

def convert_canonical_to_full_ids(input_tree,canonical_ids,output_tree):

    # we load the input tree with dendropy

    tree = dendropy.Tree.get(path=input_tree, schema="newick")

    # we create a dictionary with the canonical ids as keys and the full ids as values
    with open(canonical_ids, 'r') as f:
        for line in f:
            full_id, canonical_id = line.strip().split('\t')
            tree.replace_node_with_label(canonical_id, full_id)

    # we write the output tree with dendropy
    tree.write(path=output_tree, schema="newick")

def create_ena_xreference_table(bac120_metadata,ar53_metadata,release_number,output_dir):

    # load the bac120 metadata as a dataframe
    bac120_df = pd.read_csv(bac120_metadata,sep='\t')
    ar53_df = pd.read_csv(ar53_metadata,sep='\t')

    #only keep accession and ncbi_genbank_assembly_accession columns
    bac120_df = bac120_df[['accession','ncbi_genbank_assembly_accession']]
    ar53_df = ar53_df[['accession','ncbi_genbank_assembly_accession']]

    # combine the two dataframes
    combined_df = pd.concat([bac120_df,ar53_df])

    # in accession replace RS_ and GB_ with nothing
    combined_df['accession'] = combined_df['accession'].str.replace('RS_','').str.replace('GB_','')

    # rename the columns accession-> SOURCE_PRIMARY_ID and ncbi_genbank_assembly_accession-> TARGET_PRIMARY_ACC
    combined_df = combined_df.rename(columns={'accession':'SOURCE_PRIMARY_ID','ncbi_genbank_assembly_accession':'TARGET_PRIMARY_ACC'})

    # save the dataframe as a tsv file
    combined_df.to_csv(os.path.join(output_dir,f'gtdb_genomes_R{release_number}.tsv'),sep='\t',index=False)


