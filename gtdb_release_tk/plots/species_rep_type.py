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

class SpeciesRepTypePlot(AbstractPlot):
    """Pie plot."""

    def __init__(self, options):
        """Initialize."""
        AbstractPlot.__init__(self, options)

    def plot(self, type_strain, latinized, placeholder):
        """Create pie chart."""
        
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axis = self.fig.add_subplot(111)
        
        labels = ['Type strain', 'Latin name\n(no type strain)', 'Placeholder']
        sizes = [type_strain, latinized, placeholder]
        colors = ['#80b1d3', '#fdae6b', '#b3de69']

        wedgeprops = {'linewidth': 1, 'edgecolor': 'white'}
        textprops  = {'fontsize': self.options.label_font_size}
        axis.pie(sizes, 
                labels=labels, 
                autopct='%1.1f%%', 
                shadow=False, 
                startangle=90,
                colors=colors,
                wedgeprops=wedgeprops)
        axis.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

        #self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
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
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')
                
                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                gtdb_type_index = header.index('gtdb_type_designation')
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
                        
                    type_designation = line_split[gtdb_type_index]
                    if type_designation == 'type strain of species':
                        type_strain.add(gid)

        self.logger.info(' ...identified {:,} GTDB species representatives.'.format(len(gtdb_taxonomy)))
                            
        # determing type information for GTDB representatives
        out_prefix = f'gtdb_r{self.release_number}_sp_rep_type'
        if domain == 'Bacteria':
            out_prefix += '.bacteria'
        elif domain == 'Archaea':
            out_prefix += '.archaea'
            
        self.logger.info('Determining type information for GTDB representatives.')
        fout = open(self.output_dir / f'{out_prefix}.tsv', 'w')
        fout.write('Genome ID\tGTDB taxonomy\tGTDB species\tClassification\n')

        num_type_strain = 0
        num_latinized = 0
        num_placeholder = 0
        for gid, taxa in gtdb_taxonomy.items():
            gtdb_sp = taxa[6]
            fout.write('{}\t{}\t{}'.format(gid, ';'.join(taxa), gtdb_sp))
            
            if gid in type_strain:
                num_type_strain += 1
                fout.write('\ttype strain of species')
            elif self.latinized_species(gtdb_sp):
                num_latinized += 1
                fout.write('\tLatinized, not type strain')
            else:
                num_placeholder += 1
                fout.write('\tplaceholder')
                
            fout.write('\n')

        fout.close()
        
        # create plot
        self.logger.info('Creating plot.')
        options = AbstractPlot.Options(width=2.5, 
                                        height=2.5, 
                                        label_font_size=7, 
                                        tick_font_size=6, 
                                        dpi=600)
        plot = SpeciesRepTypePlot(options)
        plot.plot(num_type_strain, num_latinized, num_placeholder)
        
        plot.save_plot(self.output_dir / f'{out_prefix}.png', dpi=600)
        plot.save_plot(self.output_dir / f'{out_prefix}.svg', dpi=600)