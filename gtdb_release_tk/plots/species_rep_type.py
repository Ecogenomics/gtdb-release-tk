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

import logging
import os
import sys
from collections import defaultdict
from pathlib import PurePath

from biolib.plots.abstract_plot import AbstractPlot

from gtdb_release_tk.common import summarise_file
from gtdb_release_tk.files.metadata import MetadataFile


class SpeciesRepTypePlot(AbstractPlot):
    """Pie plot."""

    def __init__(self, options):
        """Initialize."""
        AbstractPlot.__init__(self, options)

    def plot(self, type_strain_categories,
             latinized_categories,
             placeholder_categories):
        """Create pie chart."""

        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axis = self.fig.add_subplot(111)

        # plot inner donut with number of species comprised of
        # just isolates, just MAGs/SAGs, or a mix of isolates 
        # and MAGs/SAGs
        sizes = []
        colors = []
        for c in [type_strain_categories, latinized_categories, placeholder_categories]:
            sizes.append(c.get('ISOLATE', 0))
            colors.append('#fdae6b')

            sizes.append(c.get('MAG', 0))
            colors.append('#b3de69')

            sizes.append(c.get('SAG', 0))
            colors.append('#80b1d3')

        wedgeprops = {'linewidth': 1,
                      'edgecolor': 'white',
                      'width': 0.6}
        patches, texts = axis.pie(sizes,
                                  shadow=False,
                                  startangle=90,
                                  colors=colors,
                                  wedgeprops=wedgeprops)
        axis.axis('equal')

        axis.legend(patches,
                    ('Isolate',
                     'MAG',
                     'SAG'),
                    fontsize=self.options.label_font_size,
                    loc='upper left',
                    bbox_to_anchor=(1, 1),
                    frameon=False)

        # plot donut with number of type strain, latinized, and
        # placeholder species clusters
        num_type_stain = sum(type_strain_categories.values())
        num_lainized = sum(latinized_categories.values())
        num_placeholder = sum(placeholder_categories.values())

        labels = [f'Type strain\n({num_type_stain:,})',
                  f'Latin, no type strain\n({num_lainized:,})',
                  f'Placeholder\n ({num_placeholder:,})']
        sizes = [num_type_stain,
                 num_lainized,
                 num_placeholder]
        colors = ['#fb8072', '#bebada', '#8dd3c7']

        wedgeprops = {'linewidth': 1,
                      'edgecolor': 'white',
                      'width': 0.4}
        textprops = {'fontsize': self.options.label_font_size}
        axis.pie(sizes,
                 labels=labels,
                 autopct='%1.1f%%',
                 pctdistance=0.8,
                 shadow=False,
                 startangle=90,
                 colors=colors,
                 wedgeprops=wedgeprops,
                 textprops=textprops)
        axis.axis('equal')

        self.draw()


class SpeciesRepType(object):
    """Bar plot comparing GTDB and NCBI taxonomies."""

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

    def run(self, bac120_metadata_file: MetadataFile, ar120_metadata_file: MetadataFile, domain):
        """Pie graph comparing GTDB and NCBI taxonomies."""

        # parse GTDB metadata file to determine genomes in each species clusters
        self.logger.info('Reading GTDB metadata.')
        gtdb_taxonomy = dict()
        type_strain = set()
        genome_category = dict()
        sp_clusters = defaultdict(set)
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            for gid, row in mf.rows.items():
                gtdb_genome_rep = row.gtdb_genome_representative
                sp_clusters[gtdb_genome_rep].add(gid)
                genome_category[gid] = row.ncbi_genome_category

                gtdb_rep = row.gtdb_representative
                if gtdb_rep != 't':
                    continue

                taxonomy = row.gtdb_taxonomy
                gtdb_domain = taxonomy.d.rank
                if domain == 'Both' or domain in gtdb_domain:
                    gtdb_taxonomy[gid] = taxonomy

                type_designation = row.gtdb_type_designation
                if type_designation == 'type strain of species':
                    type_strain.add(gid)

        self.logger.info(' ...identified {:,} GTDB species representatives.'.format(len(gtdb_taxonomy)))

        # determine genome types in each species cluster
        sp_genome_types = {}
        for rid in sp_clusters:
            if rid.startswith('UBA') or genome_category[rid] == 'derived from metagenome':
                sp_genome_types[rid] = 'MAG'
            elif genome_category[rid] in ['derived from environmental sample', 'derived from environmental_sample']:
                sp_genome_types[rid] = 'MAG'
            elif genome_category[rid] == 'derived from single cell':
                sp_genome_types[rid] = 'SAG'
            elif genome_category[rid] == 'none':
                sp_genome_types[rid] = 'ISOLATE'
            else:
                print(f'Unrecognized genome category: {genome_category[rid]}')
                sys.exit(1)

        # determine type information for GTDB representatives,
        # and genome type category (isolate, MAG/SAG)
        out_prefix = f'gtdb_r{self.release_number}_sp_rep_type'
        if domain == 'Bacteria':
            out_prefix += '.bacteria'
        elif domain == 'Archaea':
            out_prefix += '.archaea'

        self.logger.info('Determining type information for GTDB representatives.')
        path_tsv = os.path.join(self.output_dir, f'{out_prefix}.tsv')
        with open(path_tsv, 'w') as fout:
            fout.write('Genome ID\tGTDB taxonomy\tGTDB species\tClassification\n')

            type_strain_categories = defaultdict(int)
            latinized_categories = defaultdict(int)
            placeholder_categories = defaultdict(int)
            for gid, taxa in gtdb_taxonomy.items():
                gtdb_sp = taxa.s.rank
                fout.write(f'{gid}\t{taxa}\t{gtdb_sp}')

                if gid in type_strain:
                    fout.write('\ttype strain of species')
                    type_strain_categories[sp_genome_types[gid]] += 1
                elif self.latinized_species(gtdb_sp):
                    fout.write('\tLatinized, not type strain')
                    latinized_categories[sp_genome_types[gid]] += 1
                else:
                    fout.write('\tplaceholder')
                    placeholder_categories[sp_genome_types[gid]] += 1

                fout.write('\n')
        self.logger.info(summarise_file(path_tsv))

        # create plot
        self.logger.info('Creating plot.')
        options = AbstractPlot.Options(width=3,
                                       height=3,
                                       label_font_size=7,
                                       tick_font_size=6,
                                       dpi=600)
        plot = SpeciesRepTypePlot(options)
        plot.plot(type_strain_categories,
                  latinized_categories,
                  placeholder_categories)

        for ext in ('.png', '.svg'):
            path = os.path.join(self.output_dir, f'{out_prefix}{ext}')
            plot.save_plot(path, dpi=600)
            self.logger.info(summarise_file(path))
