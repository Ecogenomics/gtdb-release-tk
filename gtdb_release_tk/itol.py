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

import logging
from collections import deque
from pathlib import Path

import dendropy
import pandas as pd
import seaborn as sns
from biolib.newick import parse_label
from gtdblib.file.itol.collapse import iTolCollapseFile
from gtdblib.file.itol.dataset_color_strip import iTolDatasetColorStripFile
from gtdblib.file.itol.label import iTolLabelFile
from gtdblib.file.itol.popup import iTolPopupFile
from gtdblib.file.itol.tree_colors import iTolTreeColorsFile
from gtdblib.util.color import rgb_to_hex, TABLEAU_20


class iTOL:
    """Generate tree and iTOL files for producing iTOL tree image."""

    def __init__(self, output_dir):
        if output_dir:
            self.output_dir = Path(output_dir)

        self.logger = logging.getLogger('timestamp')

        self.colors = ['#a6cee3',
                       '#1f78b4',
                       '#b2df8a',
                       '#33a02c',
                       '#fb9a99',
                       '#e31a1c',
                       '#fdbf6f',
                       '#ff7f00',
                       '#cab2d6',
                       '#6a3d9a',
                       '#ffff99',
                       '#b15928']

    def run(self, tree_file, metadata_file):
        """Generate tree and iTOL files for producing iTOL tree image."""
        self.logger.info('Creating trees with iTOL labels.')

        self.logger.info(f'Reading GTDB reference tree: {tree_file}')
        tree = dendropy.Tree.get_from_path(tree_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        self.logger.info(f'Found {len(tree.leaf_nodes()):,} leaves in tree.')

        self.logger.info(f'Reading GTDB metadata file: {metadata_file}')
        metadata = pd.read_csv(metadata_file, sep='\t', index_col=0, low_memory=False)
        metadata['domain'] = metadata['gtdb_taxonomy'].apply(lambda x: x.split(';')[0])
        metadata['phylum'] = metadata['gtdb_taxonomy'].apply(lambda x: x.split(';')[1])
        metadata['class'] = metadata['gtdb_taxonomy'].apply(lambda x: x.split(';')[2])
        metadata['order'] = metadata['gtdb_taxonomy'].apply(lambda x: x.split(';')[3])
        metadata['family'] = metadata['gtdb_taxonomy'].apply(lambda x: x.split(';')[4])
        metadata['genus'] = metadata['gtdb_taxonomy'].apply(lambda x: x.split(';')[5])
        metadata['species'] = metadata['gtdb_taxonomy'].apply(lambda x: x.split(';')[6])
        self.logger.info(f'Found {len(metadata):,} genomes in metadata file.')

        self.logger.info(f'Getting the last common ancestor of each phylum.')
        d_phylum_to_lca = get_phylum_to_lca(tree, metadata)
        d_phylum_palette = get_phylum_colours(tree, metadata)
        write_color_datastrip(d_phylum_to_lca, d_phylum_palette, self.output_dir / 'dataset_strip_phylum.txt')

        self.logger.info('Pulling internal node labels.')
        d_int_label_to_lca = get_internal_nodes_with_labels(tree)
        self.logger.info(f'Found {len(d_int_label_to_lca):,} internal nodes with labels')

        self.logger.info('Writing node label file.')
        write_internal_node_labels(tree, d_int_label_to_lca, metadata, self.output_dir / 'labels.txt')

        self.logger.info('Writing popup information file')
        write_popup_file(tree, d_int_label_to_lca, metadata, self.output_dir / 'popup.txt')

        self.logger.info('Getting tree colour palette')
        write_tree_colours(tree, metadata, d_int_label_to_lca, self.output_dir / 'tree_colours.txt', d_phylum_palette)

        self.logger.info(f'Finding the LCA of all taxa')
        d_taxon_to_lca = get_taxon_to_lca(tree, metadata)
        write_collapse_file(d_taxon_to_lca, self.output_dir / 'collapse_phylum.txt', 'p__')
        write_collapse_file(d_taxon_to_lca, self.output_dir / 'collapse_class.txt', 'c__')
        write_collapse_file(d_taxon_to_lca, self.output_dir / 'collapse_order.txt', 'o__')
        write_collapse_file(d_taxon_to_lca, self.output_dir / 'collapse_family.txt', 'f__')
        write_collapse_file(d_taxon_to_lca, self.output_dir / 'collapse_genus.txt', 'g__')

        self.logger.info('Making tree compatible with iTOL (stripping labels)')
        strip_tree_labels(tree)
        path_tree_out = self.output_dir / 'gtdb_itol.tree'
        tree.write_to_path(path_tree_out, schema='newick', suppress_rooting=True, unquoted_underscores=True)
        return


def get_nodes_from_taxon_mrca(tree: dendropy.Tree, metadata: pd.DataFrame, rank):
    d_taxon_to_node = dict()

    # 1. Shortcut by using the internal node labels to get the MRCA of the rank
    for node in tree.preorder_node_iter():
        _, label, _ = parse_label(node.label)
        if label:
            for taxon in label.split(';'):
                if taxon.startswith(f'{rank[0]}__'):
                    d_taxon_to_node[taxon] = node
                    break

    # 2. Find any remaining nodes that weren't grabbed in the first step
    taxa_to_get = set(metadata[rank].values) - set(d_taxon_to_node.keys())
    for taxon in taxa_to_get:
        df_subset = metadata[(metadata[rank] == taxon) & (metadata['gtdb_representative'] == 't')]
        gids = set(df_subset.index)
        d_taxon_to_node[taxon] = tree.mrca(taxon_labels=gids)

    # 3. Calculate the LCA, this will most likely be used in writing files
    d_taxon_to_lca = dict()
    for taxon, node in d_taxon_to_node.items():
        d_taxon_to_lca[taxon] = get_lca(node)
    return d_taxon_to_node, d_taxon_to_lca


def get_taxon_to_lca(tree: dendropy.Tree, metadata: pd.DataFrame):
    d_taxon_to_lca = dict()

    d_prefix_to_rank = {
        'd__': 'domain',
        'p__': 'phylum',
        'c__': 'class',
        'o__': 'order',
        'f__': 'family',
        'g__': 'genus',
        's__': 'species'
    }
    all_taxa = set()
    for cur_rank in ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        all_taxa.update(set(metadata[cur_rank].unique()))

    # 1. Shortcut by using the internal node labels to get the MRCA of the rank
    for node in tree.preorder_node_iter():
        _, label, _ = parse_label(node.label)
        if label and label != 'None':
            lca = get_lca(node)
            for taxon in label.split(';'):
                d_taxon_to_lca[taxon] = lca

    # 2. Find any remaining nodes that weren't grabbed in the first step
    taxa_to_get = all_taxa - set(d_taxon_to_lca.keys())
    for taxon in taxa_to_get:
        df_subset = metadata[
            (metadata[d_prefix_to_rank[taxon[0:3]]] == taxon) & (metadata['gtdb_representative'] == 't')]
        gids = set(df_subset.index)
        if len(gids) == 1:
            d_taxon_to_lca[taxon] = [list(gids)[0]]
        else:
            d_taxon_to_lca[taxon] = get_lca(tree.mrca(taxon_labels=gids))

    return d_taxon_to_lca


def strip_tree_labels(tree: dendropy.Tree):
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        support, taxon, auxiliary_info = parse_label(node.label)
        if taxon:
            node.label = str(support)
    return


def get_lca(node: dendropy.Node):
    if node.is_leaf():
        return [node.taxon.label]
    if len(node.child_nodes()) < 2:
        return [node.child_nodes()[0].leaf_nodes()[0]]
    left = node.child_nodes()[0].leaf_nodes()[0]
    right = node.child_nodes()[1].leaf_nodes()[0]
    return [left.taxon.label, right.taxon.label]


def get_lca_str(node: dendropy.Node):
    if node.is_leaf():
        return [node.taxon.label]
    if len(node.child_nodes()) < 2:
        return [node.child_nodes()[0].leaf_nodes()[0]]
    left = node.child_nodes()[0].leaf_nodes()[0]
    right = node.child_nodes()[1].leaf_nodes()[0]
    return f'{left.taxon.label}|{right.taxon.label}'


def write_collapse_file(d_taxon_to_lca, path, taxon_prefix):
    file = iTolCollapseFile(path)
    for taxon, lr_nodes in d_taxon_to_lca.items():
        if len(lr_nodes) == 2 and taxon.startswith(taxon_prefix):
            file.insert("|".join(lr_nodes))
    file.write()


def get_colours_for_phyla(phyla):
    out = dict()
    pal = sns.color_palette("husl", len(phyla))
    for phylum, cur_pal in zip(phyla, pal):
        out[phylum] = rgb_to_hex(*[round(x * 255) for x in cur_pal])
    return out


def write_tree_colours(tree, metadata, d_int_label_to_lca, path, d_color_palette):
    file = iTolTreeColorsFile(path)

    rank_idxs = ('p__', 'c__', 'o__', 'f__', 'g__', 's__')

    # 1. Map each taxon to it's higher rank
    d_taxon_to_phylum = dict()
    for _, row in metadata.iterrows():
        for taxon in row['gtdb_taxonomy'].split(';'):
            d_taxon_to_phylum[taxon] = row['phylum']

    # 2. Colour the internal node labels by the highest degree of tax novelty
    for node in tree.preorder_node_iter():
        _, label, _ = parse_label(node.label)
        if label and label != 'None':
            if label.startswith('d__'):
                continue
            highest_rank = label.split('; ')[0]
            if highest_rank.startswith('p__'):
                colour = d_color_palette[highest_rank][0]
            else:
                colour = d_color_palette[d_taxon_to_phylum[highest_rank]][rank_idxs.index(highest_rank[0:3])]
            file.insert_range('|'.join(d_int_label_to_lca[label]), colour, label)

    # 3. Add leaf node labels (note: go up if this is a singleton)
    for leaf_node in tree.leaf_node_iter():
        gid = leaf_node.taxon.label

        parent = leaf_node.parent_node
        highest_taxon = None
        while parent is not None:
            _, label, _ = parse_label(parent.label)
            if label is not None and label != 'None':
                highest_taxon = label.split('; ')[-1]
                break
            parent = parent.parent_node

        if highest_taxon is None:
            raise Exception(f'Unable to determine highest taxon for {gid}')
        if highest_taxon.startswith('s__'):
            continue

        # Append the taxon string up to this rank
        if highest_taxon.startswith('d__'):
            highest_taxon = metadata.loc[gid]['phylum']
            colour = d_color_palette[highest_taxon][0]
        else:
            colour = d_color_palette[d_taxon_to_phylum[highest_taxon]][rank_idxs.index(highest_taxon[0:3])]
        file.insert_range(gid, colour, highest_taxon)

    file.write()


def write_internal_node_labels(tree: dendropy.Tree, d_label_to_lca, metadata: pd.DataFrame, path: Path):
    file = iTolLabelFile(path)

    # 1. Add internal node labels (trivial)
    for label, lr_nodes in d_label_to_lca.items():
        file.insert("|".join(lr_nodes), label)

    # 2. Add leaf node labels (note: go up if this is a singleton)
    df_sp_reps = metadata[metadata['gtdb_representative'] == 't']
    for leaf_node in tree.leaf_node_iter():
        gid = leaf_node.taxon.label

        parent = leaf_node.parent_node
        highest_taxon = None
        while parent is not None:
            _, label, _ = parse_label(parent.label)
            if label is not None and label != 'None':
                highest_taxon = label.split('; ')[-1]
                break
            parent = parent.parent_node

        if highest_taxon is None:
            raise Exception(f'Unable to determine highest taxon for {gid}')

        # Append the taxon string up to this rank
        taxonomy = df_sp_reps.loc[leaf_node.taxon.label, 'gtdb_taxonomy'].split(';')
        taxonomy_string = '; '.join(taxonomy[taxonomy.index(highest_taxon) + 1:])
        file.insert(gid, taxonomy_string)

    file.write()


def get_internal_nodes_with_labels(tree):
    out = dict()
    for node in tree.preorder_node_iter():
        _, label, _ = parse_label(node.label)
        if label and label != 'None':
            lca = get_lca(node)
            out[label] = lca
    return out


def get_n_genomes_under_each_taxon(metadata: pd.DataFrame):
    out = dict()
    for gid, row in metadata.iterrows():
        for taxon in row['gtdb_taxonomy'].split(';'):
            if taxon not in out:
                out[taxon] = 0
            out[taxon] += 1
    return out


def write_popup_file(tree: dendropy.Tree, d_int_label_to_lca, metadata: pd.DataFrame, path: Path):
    file = iTolPopupFile(path)

    taxon_to_cnt = get_n_genomes_under_each_taxon(metadata)

    for leaf_node in tree.leaf_node_iter():
        gid = leaf_node.taxon.label
        gid_short = gid[3:]
        row = metadata.loc[gid]
        branch_length = leaf_node.edge_length

        lines = [
            '<b>Representative</b>: ' + f'<a href="https://gtdb.ecogenomic.org/genome?gid={gid_short}" target="_blank">{gid}</a>',
            '<br><b>Genomes in species</b>: ' + f'{taxon_to_cnt[row["species"]]:,}',
            '<br><b>Branch Length</b>: ' + str(branch_length),
            '<br><b>Domain</b>: ' + row['domain'],
            '<br><b>Phylum</b>: ' + row['phylum'],
            '<br><b>Class</b>: ' + row['class'],
            '<br><b>Order</b>: ' + row['order'],
            '<br><b>Family</b>: ' + row['family'],
            '<br><b>Genus</b>: ' + row['genus'],
            '<br><b>Species</b>: ' + row['species'],
        ]
        file.insert(gid, row['species'], ''.join(lines))

    # # Label the internal nodes
    # for node in tree.preorder_node_iter():
    #     _, label, _ = parse_label(node.label)
    #     if label and label != 'None':
    #         lca = d_int_label_to_lca[label]
    #         taxonomy = label.split('; ')
    #         lines = list()
    #         for taxon in taxonomy:
    #             lines.append(f'<b>{taxon}</b>: {taxon_to_cnt[taxon]:,}')
    #         lines.extend([
    #             '<br><b>Branch Length</b>: ' + str(node.edge_length),
    #         ])
    #         file.insert("|".join(lca), taxonomy[0], ''.join(lines))

    file.write()
    return


def get_phylum_palette(metadata: pd.DataFrame):
    # 1. Assign a colour to each phylum
    phyla = sorted(metadata['phylum'].unique())
    cmap = sns.color_palette('Spectral', len(phyla))

    # 2. For each phylum, create a colour map to species
    out = dict()
    for phylum, phylum_color in zip(phyla, cmap):
        phylum_cmap = sns.light_palette(phylum_color, 6, reverse=True)
        out[phylum] = [rgb_to_hex(*[round(y * 255) for y in x]) for x in phylum_cmap]
    return out


def get_tree_colour_palette(tree: dendropy.Tree, metadata: pd.DataFrame):
    # 1. Assign a colour to each phylum
    phyla = sorted(metadata['phylum'].unique())
    cmap = sns.color_palette('Spectral', len(phyla))

    d_phyla_to_colour = dict()
    for phylum, cur_cmap in zip(phyla, cmap):
        d_phyla_to_colour[phylum] = dict()

        cur_hex = rgb_to_hex(*[round(x * 255) for x in cur_cmap])
        sub_ranks = ('p__', 'c__', 'o__', 'f__', 'g__')
        sub_rank_cmap = sns.light_palette(cur_hex, len(sub_ranks), reverse=True)
        for sub_rank, cur_sub_cmap in zip(sub_ranks, sub_rank_cmap):
            d_phyla_to_colour[phylum][sub_rank] = rgb_to_hex(*[round(x * 255) for x in cur_sub_cmap])

    # 2. Map each taxon to it's higher rank
    d_taxon_to_phylum = dict()
    for _, row in metadata.iterrows():
        for taxon in row['gtdb_taxonomy'].split(';'):
            d_taxon_to_phylum[taxon] = row['phylum']

    # Label the internal nodes
    out = dict()
    for node in tree.preorder_node_iter():
        _, label, _ = parse_label(node.label)
        if label and label != 'None':
            if label.startswith('d__'):
                continue
            highest_rank = label.split('; ')[0]
            if highest_rank.startswith('p__'):
                out[label] = d_phyla_to_colour[highest_rank]['p__']
            else:
                out[label] = d_phyla_to_colour[d_taxon_to_phylum[highest_rank]][highest_rank[:3]]
    return out


def get_phylum_colours(tree: dendropy.Tree, metadata: pd.DataFrame):
    """Generate a list of colours for each phylum. Order is inferred through a
    depth first search, to ensure no clade colours are side by side.
    """
    colours = TABLEAU_20

    d_phylum_to_colour = dict()
    queue = deque(tree.seed_node.child_nodes())
    while len(queue) > 0:
        node = queue.pop()

        # Case 1. It's a leaf node, i.e. novel phylum
        if node.is_leaf():
            cur_phylum = metadata.loc[node.taxon.label]['phylum']
            cur_colour = colours[len(d_phylum_to_colour) % len(colours)]
            d_phylum_to_colour[cur_phylum] = cur_colour

        # Case 2. There is a label?
        else:
            _, label, _ = parse_label(node.label)

            # Use the label, add the phylum and stop adding children
            if label and label.startswith('p__'):
                cur_phylum = label.split('; ')[0]
                cur_colour = colours[len(d_phylum_to_colour) % len(colours)]
                d_phylum_to_colour[cur_phylum] = cur_colour

            # Otherwise, add the children and keep searching
            else:
                for child in node.child_nodes():
                    queue.append(child)

    # Generate a colour palette for each phylum (increasing brightness)
    out = dict()
    for phylum, colour in d_phylum_to_colour.items():
        cur_pal = sns.light_palette(colour, 6, reverse=True)
        out[phylum] = [rgb_to_hex(*[round(y * 255) for y in x]) for x in cur_pal]
    return out


def get_phylum_to_lca(tree: dendropy.Tree, metadata: pd.DataFrame):
    """Calculate the LCA for each phylum. Considers singletons."""
    out = dict()

    queue = deque(tree.seed_node.child_nodes())
    while len(queue) > 0:
        node = queue.pop()

        # Case 1. It's a leaf node, i.e. novel phylum
        if node.is_leaf():
            cur_phylum = metadata.loc[node.taxon.label]['phylum']
            out[cur_phylum] = node.taxon.label

        # Case 2. There is a label?
        else:
            _, label, _ = parse_label(node.label)

            # Use the label, add the phylum and stop adding children
            if label and label.startswith('p__'):
                cur_phylum = label.split('; ')[0]
                out[cur_phylum] = get_lca_str(node)

            # Otherwise, add the children and keep searching
            else:
                for child in reversed(node.child_nodes()):
                    queue.append(child)
    return out


def write_color_datastrip(d_phylum_to_lca, d_phylum_palette, path):
    file = iTolDatasetColorStripFile(path, 'Phylum Labels', '#000000', strip_width=100, show_internal=True)

    for phylum, node_id in d_phylum_to_lca.items():
        colour = d_phylum_palette[phylum][0]

        file.insert(node_id, colour, phylum)

    file.write()
    return
