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

class GenomeTypesPerRankPlot(AbstractPlot):
    """Box and whisker plot."""

    def __init__(self, options):
        """Initialize."""
        AbstractPlot.__init__(self, options)

    def plot(self, both, isolate, env, xticklabels):
        """Create stacked bar plot."""
        
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axis = self.fig.add_subplot(111)
        
        ind = np_arange(len(both))
        width = 0.7
        
        both = np_array(both)
        isolate = np_array(isolate)
        env = np_array(env)

        p1 = axis.bar(ind, both, width)
        p2 = axis.bar(ind, isolate, width, bottom=both)
        p3 = axis.bar(ind, env, width, bottom=both+isolate)

        axis.set_ylim([0, 100])
        axis.set_yticks(range(0, 101, 10))
        axis.set_ylabel('Taxa (%)')
        
        axis.set_xticks(ind)
        axis.set_xticklabels(xticklabels)
        
        axis.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
        axis.set_axisbelow(True)

        self.prettify(axis)
        
        axis.legend((p3[0], p2[0], p1[0]), ('Exclusively MAGs and/or SAGs', 
                                            'Exclusively isolates', 
                                            'Isolate and environmental genomes'),
                                            fontsize=self.options.tick_font_size)

        #self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()


class GenomeTypesPerRank(object):
    """Plot number of MAGs, SAGs, and isolates for each taxonomic rank."""
    
    def __init__(self, output_dir):
        """Initialization."""
        
        self.output_dir = PurePath(output_dir)
        
        self.logger = logging.getLogger('timestamp')
        
    def run(self, metadata_file):
        """Plot number of MAGs, SAGs, and isolates for each taxonomic rank."""
        
        # parse GTDB metadata file to determine genomes in each species clusters
        # and the type of these genomes (MAG, SAG, or isolate)
        self.logger.info('Reading GTDB metadata.')
        gtdb_taxonomy = {}
        sp_clusters = defaultdict(set)
        genome_category = {}
        with open(metadata_file, encoding='utf-8') as f:
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
        env_categories = set(['derived from single cell',
                                'derived from metagenome', 
                                'derived from environmental_sample'])
        sp_genome_types = defaultdict(set)
        for sp, gids in sp_clusters.items():
            for gid in gids:
                if genome_category[gid] in env_categories:
                    sp_genome_types[sp].add('ENV')
                elif genome_category[gid] == 'none':
                    sp_genome_types[sp].add('ISOLATE')
                else:
                    self.logger.error('Unrecognized genome category: {genome_category[gid]}')
                    sys.exit(-1)
                    
        # get species in each taxa
        sp_in_taxa = defaultdict(lambda: defaultdict(set))
        for taxa in gtdb_taxonomy.values():
            for rank_index in range(1, 7):
                sp_in_taxa[rank_index][taxa[rank_index]].add(taxa[6])
                    
        # tabulate number of genome types at each rank
        self.logger.info('Tabulating genomes types at each rank.')
        fout_count = open(self.output_dir / 'genome_types_per_rank.counts.tsv', 'w')
        fout_count.write('Rank\tNo. taxa\tBoth\tIsolate\tEnvironmental\n')
        
        fout_taxa = open(self.output_dir / 'genome_types_per_rank.taxa.tsv', 'w')
        fout_taxa.write('Rank\tNo. taxa\tBoth\tIsolate\tEnvironmental\n')
        
        plot_both = []
        plot_isolate = []
        plot_end = []
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
                    taxon_categories.update(sp_genome_types[sp])
                    
                if 'ENV' in taxon_categories and 'ISOLATE' in taxon_categories:
                    both.add(taxon)
                elif 'ISOLATE' in taxon_categories:
                    isolate.add(taxon)
                elif 'ENV' in taxon_categories:
                    env.add(taxon)
                else:
                    self.logger.error('Genomes in species have an unassigned category.')
                    sys.exit(-1)
                
            total_taxa = len(both) + len(isolate) + len(env)
            fout_count.write(f'\t{total_taxa}\t{len(both)}\t{len(isolate)}\t{len(env)}\n')
            
            fout_taxa.write('\t{}\t{}\t{}\t{}\n'.format(
                            total_taxa,
                            ', '.join(sorted(both)),
                            ', '.join(sorted(isolate)),
                            ', '.join(sorted(env))))
                            
            plot_both.append(len(both)*100.0/total_taxa)
            plot_isolate.append(len(isolate)*100.0/total_taxa)
            plot_end.append(len(env)*100.0/total_taxa)
            plot_labels.append('{}\n{:,}'.format(
                                Taxonomy.rank_labels[rank_index].capitalize(), 
                                total_taxa))
        
        isolate_genomes = sum([1 for c in genome_category.values() if c not in env_categories])
        env_genomes = sum([1 for c in genome_category.values() if c in env_categories])
        fout_count.write('Genomes\t{}\t{}\t{}\t{}\n'.format(
                            len(genome_category),
                            0,
                            isolate_genomes,
                            env_genomes))
                    
        plot_both.append(0)
        plot_isolate.append(isolate_genomes*100.0/len(genome_category))
        plot_end.append(env_genomes*100.0/len(genome_category))
        plot_labels.append('{}\n{:,}'.format(
                            'Genomes', 
                            len(genome_category)))
            
        fout_count.close()
        fout_taxa.close()
        
        # create plot
        self.logger.info('Creating plot.')
        options = AbstractPlot.Options(width=4.5, 
                                        height=3, 
                                        label_font_size=7, 
                                        tick_font_size=6, 
                                        dpi=600)
        plot = GenomeTypesPerRankPlot(options)
        plot.plot(plot_both, plot_isolate, plot_end, plot_labels)
        
        plot.save_plot(self.output_dir / 'genome_types_per_rank.png', dpi=600)
        plot.save_plot(self.output_dir / 'genome_types_per_rank.svg', dpi=600)