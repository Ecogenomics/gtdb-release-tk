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
import logging
from pathlib import PurePath
from collections import defaultdict

from biolib.taxonomy import Taxonomy
from biolib.plots.abstract_plot import AbstractPlot

from numpy import (arange as np_arange,
                   array as np_array)

from gtdb_release_tk.common import summarise_file
from gtdb_release_tk.plots.palette import Palette, DEFAULT_PALETTE


class NomenclaturalPerRankPlot(AbstractPlot):
    """Box and whisker plot."""

    def __init__(self, options):
        """Initialize."""
        
        AbstractPlot.__init__(self, options)

    def plot(self, plot_latinized, plot_placeholder, xticklabels, palette: Palette = DEFAULT_PALETTE):
        """Create stacked bar plot."""

        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axis = self.fig.add_subplot(111)

        ind = np_arange(len(plot_latinized))
        width = 0.7

        plot_latinized = np_array(plot_latinized)
        plot_placeholder = np_array(plot_placeholder)

        p1 = axis.bar(ind, plot_latinized, width, color=palette.colour1)
        p2 = axis.bar(ind, plot_placeholder, width, bottom=plot_latinized, color=palette.colour2)

        axis.set_ylim([0, 100])
        axis.set_yticks(range(0, 101, 10))
        axis.set_ylabel('Taxa (%)')

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

        axis.legend((p2[0], p1[0]), ('Placeholder', 'Latinized'),
                    fontsize=self.options.tick_font_size,
                    loc='upper left',
                    bbox_to_anchor=(1, 1),
                    frameon=False)

        #self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()


class NomenclaturalPerRank(object):
    """Plot nomenclatural status of species for each taxonomic rank."""

    def __init__(self, release_number, output_dir):
        """Initialization."""

        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def latinized_taxon(self, taxon):
        """Check if taxon is a Latinized name."""

        if '__' in taxon:
            taxon = taxon[3:]

        return (taxon[0].isupper()
                and all(c.islower() for c in taxon[1:]))

    def latinized_species(self, species_name):
        """Check if species name is Latinized."""

        generic, specific = species_name.split()
        generic = generic.replace('s__', '')

        if not self.latinized_taxon(generic):
            return False

        return all(c.islower() for c in specific)

    def run(self, bac120_metadata_file, ar120_metadata_file, domain, palette: Palette = DEFAULT_PALETTE):
        """Plot nomenclatural status of species for each taxonomic rank."""

        # parse GTDB metadata file to determine genomes in each species clusters
        self.logger.info('Reading GTDB metadata.')
        gtdb_taxonomy = {}
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                gtdb_rep_index = header.index('gtdb_representative')

                for line in f:
                    line_split = line.strip().split('\t')

                    gid = line_split[0]

                    gtdb_rep = line_split[gtdb_rep_index]
                    if gtdb_rep != 't':
                        continue

                    taxonomy = line_split[gtdb_taxonomy_index]
                    gtdb_taxa = [t.strip() for t in taxonomy.split(';')]

                    gtdb_domain = gtdb_taxa[0]
                    if domain == 'Both' or domain in gtdb_domain:
                        gtdb_taxonomy[gid] = gtdb_taxa

        self.logger.info(
            ' ...identified {:,} representative genomes.'.format(len(gtdb_taxonomy)))

        # get GTDB taxa at each rank
        gtdb_taxa_at_rank = defaultdict(set)
        for taxa in gtdb_taxonomy.values():
            for rank, taxon in enumerate(taxa):
                gtdb_taxa_at_rank[rank].add(taxon)

        # determine nomenclatural category for each taxa at each rank
        out_prefix = f'gtdb_r{self.release_number}_nomenclatural_per_rank'
        self.logger.info('Determining nomenclatural type of taxa.')
        path_tsv = os.path.join(self.output_dir, f'{out_prefix}.species.tsv')
        with open(path_tsv, 'w') as fout:
            plot_latinized = []
            plot_placeholder = []
            plot_labels = []
            for rank_index in range(1, 7):
                fout.write(Taxonomy.rank_labels[rank_index])

                latinized = 0
                placeholder = 0
                for taxon in gtdb_taxa_at_rank[rank_index]:
                    if rank_index != 6:
                        if self.latinized_taxon(taxon):
                            latinized += 1
                        else:
                            placeholder += 1
                    else:
                        if self.latinized_species(taxon):
                            latinized += 1
                        else:
                            placeholder += 1

                total_taxa = latinized + placeholder
                fout.write(f'\t{total_taxa}\t{latinized}\t{placeholder}\n')

                plot_latinized.append(latinized*100.0/total_taxa)
                plot_placeholder.append(placeholder*100.0/total_taxa)
                plot_labels.append('{}\n{:,}'.format(
                                    Taxonomy.rank_labels[rank_index].capitalize(),
                                    total_taxa))
        self.logger.info(summarise_file(path_tsv))

        # create plot
        self.logger.info('Creating plot.')
        options = AbstractPlot.Options(width=4,
                                       height=3,
                                       label_font_size=7,
                                       tick_font_size=6,
                                       dpi=600)
        plot = NomenclaturalPerRankPlot(options)
        plot.plot(plot_latinized, plot_placeholder, plot_labels, palette)

        for ext in ('.png', '.svg'):
            path = os.path.join(self.output_dir, f'{out_prefix}{ext}')
            plot.save_plot(path, dpi=600)
            self.logger.info(summarise_file(path))
