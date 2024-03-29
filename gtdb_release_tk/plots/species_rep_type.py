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

import sys
import logging
from pathlib import PurePath
from collections import defaultdict

from biolib.taxonomy import Taxonomy
from biolib.plots.abstract_plot import AbstractPlot

from numpy import (arange as np_arange,
                    array as np_array)
                    
from gtdb_release_tk.common import ENV_CATEGORIES


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
        textprops  = {'fontsize': self.options.label_font_size}
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
        
    def run(self, 
                bac120_metadata_file, 
                ar120_metadata_file,
                domain):
        """Bar plot comparing GTDB and NCBI taxonomies."""
        
        # parse GTDB metadata file to determine genomes in each species clusters
        self.logger.info('Reading GTDB metadata.')
        gtdb_taxonomy = {}
        type_strain = set()
        genome_category = {}
        sp_clusters = defaultdict(set)
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')
                
                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                gtdb_type_index = header.index('gtdb_type_designation')
                gtdb_rep_index = header.index('gtdb_representative')
                gtdb_genome_rep_index = header.index('gtdb_genome_representative')
                genome_category_index = header.index('ncbi_genome_category')

                for line in f:
                    line_split = line.strip().split('\t')
                    
                    gid = line_split[0]
                    
                    gtdb_genome_rep = line_split[gtdb_genome_rep_index]
                    sp_clusters[gtdb_genome_rep].add(gid)
                    genome_category[gid] = line_split[genome_category_index]

                    gtdb_rep = line_split[gtdb_rep_index]
                    if gtdb_rep != 't':
                        continue
                        
                    taxonomy = line_split[gtdb_taxonomy_index]
                    gtdb_taxa = [t.strip() for t in taxonomy.split(';')]
                    gtdb_domain = gtdb_taxa[0]
                    if domain == 'Both' or domain in gtdb_domain:
                        gtdb_taxonomy[gid] = gtdb_taxa
                        
                    type_designation = line_split[gtdb_type_index]
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
                sys.exit(-1)
                            
        # determine type information for GTDB representatives,
        # and genome type category (isolate, MAG/SAG)
        out_prefix = f'gtdb_r{self.release_number}_sp_rep_type'
        if domain == 'Bacteria':
            out_prefix += '.bacteria'
        elif domain == 'Archaea':
            out_prefix += '.archaea'
            
        self.logger.info('Determining type information for GTDB representatives.')
        fout = open(self.output_dir / f'{out_prefix}.tsv', 'w')
        fout.write('Genome ID\tGTDB taxonomy\tGTDB species\tClassification\n')

        type_strain_categories = defaultdict(int)
        latinized_categories = defaultdict(int)
        placeholder_categories = defaultdict(int)
        for gid, taxa in gtdb_taxonomy.items():
            gtdb_sp = taxa[6]
            fout.write('{}\t{}\t{}'.format(gid, ';'.join(taxa), gtdb_sp))
            
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

        fout.close()
        
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
        
        plot.save_plot(self.output_dir / f'{out_prefix}.png', dpi=600)
        plot.save_plot(self.output_dir / f'{out_prefix}.svg', dpi=600)