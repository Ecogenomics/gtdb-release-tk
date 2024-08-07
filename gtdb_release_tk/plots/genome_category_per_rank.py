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
import logging
from pathlib import PurePath
from collections import defaultdict

from biolib.taxonomy import Taxonomy
from biolib.plots.abstract_plot import AbstractPlot

from numpy import (arange as np_arange,
                   array as np_array)

from matplotlib import gridspec

from gtdb_release_tk.common import (ENV_CATEGORIES,
                                    sp_cluster_type_category,
                                    summarise_file)
from gtdb_release_tk.taxon_utils import canonical_taxon
from gtdb_release_tk.plots.palette import Palette, DEFAULT_PALETTE


class GenomeCateogryPerRankPlot(AbstractPlot):
    """Box and whisker plot."""

    def __init__(self, options):
        """Initialize."""
        
        AbstractPlot.__init__(self, options)

        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        self.spec = gridspec.GridSpec(
            ncols=2, nrows=1, width_ratios=[8, 1], wspace=0.5)

    def plot(self, plot_num, both, isolate, env, xticklabels, ylabel, palette: Palette):
        """Create stacked bar plot."""

        axis = self.fig.add_subplot(self.spec[plot_num-1])

        ind = np_arange(len(both))
        width = 0.7

        both = np_array(both)
        isolate = np_array(isolate)
        env = np_array(env)

        p1 = axis.bar(ind, both, width, color=palette.colour1)
        p2 = axis.bar(ind, isolate, width, bottom=both, color=palette.colour2)
        p3 = axis.bar(ind, env, width, bottom=both+isolate, color=palette.colour3)

        axis.set_ylim([0, 100])
        axis.set_yticks(range(0, 101, 10))
        axis.set_ylabel(ylabel)

        axis.set_xticks(ind)
        axis.set_xticklabels(xticklabels)

        axis.yaxis.grid(True,
                        linestyle='-',
                        which='major',
                        color='lightgrey',
                        alpha=0.7,
                        zorder=1)
        axis.set_axisbelow(True)

        self.prettify(axis)

        if plot_num == 2:
            axis.legend((p3[0], p2[0], p1[0]), ('Exclusively MAGs and/or SAGs',
                                                'Exclusively isolates',
                                                'Isolate and environmental genomes'),
                        fontsize=self.options.tick_font_size,
                        loc='upper left',
                        bbox_to_anchor=(1, 1),
                        frameon=False)

        #self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()


class GenomeCategoryPerRank(object):
    """Plot number of MAGs, SAGs, and isolates for each taxonomic rank."""

    def __init__(self, release_number, output_dir):
        """Initialization."""

        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def run(self, bac120_metadata_file, ar120_metadata_file, ar_only, bac_only, palette: Palette = DEFAULT_PALETTE):
        """Plot number of MAGs, SAGs, and isolates for each taxonomic rank."""

        # parse GTDB metadata file to determine genomes in each species clusters
        # and the type of these genomes (MAG, SAG, or isolate)
        self.logger.info('Reading GTDB metadata.')
        gtdb_taxonomy = {}
        sp_clusters = defaultdict(set)
        genome_category = {}

        metadata_files = [bac120_metadata_file, ar120_metadata_file]
        if ar_only and not bac_only:
            metadata_files = [ar120_metadata_file]
        if not ar_only and bac_only:
            metadata_files = [bac120_metadata_file]

        for mf in metadata_files:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                gtdb_rep_index = header.index('gtdb_genome_representative')
                genome_category_index = header.index('ncbi_genome_category')

                for line in f:
                    line_split = line.strip().split('\t')

                    gid = line_split[0]

                    taxonomy = line_split[gtdb_taxonomy_index]
                    gtdb_rep = line_split[gtdb_rep_index]
                    if taxonomy == 'none' or gtdb_rep == 'none':
                        continue

                    gtdb_taxa = [t.strip() for t in taxonomy.split(';')]
                    gtdb_taxonomy[gid] = gtdb_taxa

                    sp_clusters[gtdb_taxa[6]].add(gid)

                    genome_category[gid] = line_split[genome_category_index]

        self.logger.info(' ...identified {:,} species clusters spanning {:,} genomes.'.format(
            len(sp_clusters),
            sum([len(gids) for gids in sp_clusters.values()])))

        # determine genome types in each species cluster
        sp_genome_types = {}
        for sp, gids in sp_clusters.items():
            sp_genome_types[sp] = sp_cluster_type_category(
                gids, genome_category)

        # get species in each taxa
        sp_in_taxa = defaultdict(lambda: defaultdict(set))
        for taxa in gtdb_taxonomy.values():
            for rank_index in range(1, 7):
                cur_taxon = taxa[rank_index]
                if rank_index < 5:
                    # canonicalize names above genus
                    cur_taxon = canonical_taxon(cur_taxon)
                sp_in_taxa[rank_index][cur_taxon].add(taxa[6])

        # tabulate number of genome types at each rank
        out_prefix = f'gtdb_r{self.release_number}_genome_category_per_rank'
        self.logger.info('Tabulating genomes types at each rank.')
        fout_count = open(self.output_dir / f'{out_prefix}.tsv', 'w')
        fout_count.write('Rank\tNo. taxa\tBoth\tIsolate\tEnvironmental\n')

        fout_taxa = open(self.output_dir / f'{out_prefix}.taxa.tsv', 'w')
        fout_taxa.write('Rank\tNo. taxa\tBoth\tIsolate\tEnvironmental\n')

        plot_both = []
        plot_isolate = []
        plot_env = []
        plot_labels = []
        for rank_index in range(1, 7):
            fout_count.write(Taxonomy.rank_labels[rank_index])
            fout_taxa.write(Taxonomy.rank_labels[rank_index])

            both = set()
            env = set()
            isolate = set()
            for taxon in sp_in_taxa[rank_index]:
                taxon_categories = set()
                for sp in sp_in_taxa[rank_index][taxon]:
                    taxon_categories.add(sp_genome_types[sp])

                if (('ENV' in taxon_categories and 'ISOLATE' in taxon_categories)
                        or 'BOTH' in taxon_categories):
                    both.add(taxon)
                elif 'ISOLATE' in taxon_categories:
                    isolate.add(taxon)
                elif 'ENV' in taxon_categories:
                    env.add(taxon)
                else:
                    self.logger.error(
                        f'Genomes in species have an unassigned category: {taxon_categories}')
                    sys.exit(-1)

            total_taxa = len(both) + len(isolate) + len(env)
            fout_count.write(
                f'\t{total_taxa}\t{len(both)}\t{len(isolate)}\t{len(env)}\n')

            fout_taxa.write('\t{}\t{}\t{}\t{}\n'.format(
                            total_taxa,
                            ', '.join(sorted(both)),
                            ', '.join(sorted(isolate)),
                            ', '.join(sorted(env))))

            plot_both.append(len(both)*100.0/total_taxa)
            plot_isolate.append(len(isolate)*100.0/total_taxa)
            plot_env.append(len(env)*100.0/total_taxa)
            plot_labels.append('{}\n{:,}'.format(
                Taxonomy.rank_labels[rank_index].capitalize(),
                total_taxa))

        # create plot for taxa at each rank
        self.logger.info('Creating plot.')
        options = AbstractPlot.Options(width=4,
                                       height=3,
                                       label_font_size=7,
                                       tick_font_size=6,
                                       dpi=600)
        plot = GenomeCateogryPerRankPlot(options)
        plot.plot(1, plot_both, plot_isolate, plot_env,
                  plot_labels, ylabel='Taxa (%)', palette=palette)

        # create seperate plot with genomes
        isolate_genomes = sum(
            [1 for c in genome_category.values() if c not in ENV_CATEGORIES])
        env_genomes = sum(
            [1 for c in genome_category.values() if c in ENV_CATEGORIES])
        fout_count.write('Genomes\t{}\t{}\t{}\t{}\n'.format(
            len(genome_category),
            0,
            isolate_genomes,
            env_genomes))
        fout_count.close()
        fout_taxa.close()

        plot_both = [0]
        plot_isolate = [isolate_genomes*100.0/len(genome_category)]
        plot_env = [env_genomes*100.0/len(genome_category)]
        plot_label = ['{}\n{:,}'.format('Genomes', len(genome_category))]
        plot.plot(2, plot_both, plot_isolate, plot_env,
                  plot_label, ylabel='Genomes (%)', palette=palette)

        for ext in ('.png', '.svg'):
            path = os.path.join(self.output_dir, f'{out_prefix}{ext}')
            plot.save_plot(path, dpi=600)
            self.logger.info(summarise_file(path))
