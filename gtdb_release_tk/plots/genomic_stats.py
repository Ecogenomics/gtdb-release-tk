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

import os
import sys
import logging
import argparse
import ntpath
import csv
import random
import operator
from pathlib import PurePath
from collections import defaultdict, namedtuple

from biolib.plots.abstract_plot import AbstractPlot

from numpy import (median as np_median,
                    mean as np_mean,
                    std as np_std,
                    percentile as np_percentile,
                    arange as np_arange,
                    array as np_array,
                    ones_like as np_ones_like,
                    bincount as np_bincount)
                    
from matplotlib.ticker import FuncFormatter
from pylab import MaxNLocator
from matplotlib.offsetbox import AnchoredText

import gtdb_release_tk.HTML as HTML

class GenomicStatsPlot(AbstractPlot):
    """Create histogram of common genomic statistics."""

    def __init__(self, options, rows, cols):
        """Initialize."""
        
        AbstractPlot.__init__(self, options)
        
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        self.rows = rows
        self.cols = cols

    def plot(self, plot_num, data, xlabel, ylabel, txt_position=None):
        """Create histogram for statistic."""
        
        self.axis = self.fig.add_subplot(self.rows, self.cols, plot_num)

        weights = np_ones_like(data)/float(len(data))
        counts, bins, patches = self.axis.hist(data, 
                                                bins=min(20, len(set(data))-1), 
                                                rwidth=0.9, 
                                                weights=weights, 
                                                color='#fdae6b',
                                                align='mid')

        self.axis.set_xlabel(xlabel)
        self.axis.set_ylabel(ylabel)
        
        self.axis.xaxis.set_major_locator(MaxNLocator(integer=True))
        self.axis.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.1%}'.format(y)))
        
        # report summary statistics
        stat_txt = f'median = {np_median(data):.1f}\n'
        stat_txt += f'mean = {np_mean(data):.1f}\n'
        stat_txt += f'std = {np_std(data):.1f}'
        if txt_position == 'left':
            self.axis.text(0.05, 0.95, 
                            stat_txt, 
                            transform=self.axis.transAxes,
                            fontsize=self.options.tick_font_size,
                            verticalalignment='top')
        elif txt_position == 'right':
            self.axis.text(0.95, 0.95, 
                            stat_txt, 
                            transform=self.axis.transAxes,
                            fontsize=self.options.tick_font_size,
                            verticalalignment='top',
                            horizontalalignment='right')

        self.prettify(self.axis)
        for loc, spine in self.axis.spines.items():
            if loc in ['right', 'top']:
                spine.set_color('none')

        self.fig.tight_layout(pad=0.1, w_pad=1.0, h_pad=1.0)
        self.draw()
   
class GenomicStats(object):
    """Plot of common genomic statistics."""

    def __init__(self, release_number, output_dir):
        
        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')
        
    def read_metadata(self, bac120_metadata_file, ar120_metadata_file):
        """Read metadata for genomes."""
        
        Metadata = namedtuple('Metadata', ['gtdb_representative',
                                            'gtdb_type_strain',
                                            'genome_quality',
                                            'genome_comp',
                                            'genome_cont',
                                            'genome_size',
                                            'gc_percentage',
                                            'contig_count',
                                            'protein_count',
                                            'coding_density',
                                            'ssu_count',
                                            'gtdb_taxonomy',
                                            'ncbi_taxonomy'])

        metadata = {}
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                genome_index = header.index('accession')
                comp_index = header.index('checkm_completeness')
                cont_index = header.index('checkm_contamination')
                genome_size_index = header.index('genome_size')
                gc_percentage_index = header.index('gc_percentage')
                contig_count_index = header.index('contig_count')
                protein_count_index = header.index('protein_count')
                coding_density_index = header.index('coding_density')
                ssu_count_index = header.index('ssu_count')
                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                ncbi_taxonomy_index = header.index('ncbi_taxonomy')
                gtdb_rep_index = header.index('gtdb_representative')
                gtdb_type_index = header.index('gtdb_type_designation')
                    
                for line in f:
                    line_split = line.strip().split('\t')
                    
                    gid = line_split[genome_index]

                    comp = float(line_split[comp_index])
                    cont = float(line_split[cont_index])
                    quality = comp-5*cont
                        
                    ncbi_taxonomy = line_split[ncbi_taxonomy_index]
                    gtdb_taxonomy = [x.strip() for x in line_split[gtdb_taxonomy_index].split(';')]

                    gtdb_representative = False
                    if line_split[gtdb_rep_index] == 't':
                        gtdb_representative = True
                       
                    gtdb_type_strain = False
                    if line_split[gtdb_type_index] == 'type strain of species':
                        gtdb_type_strain = True
                        
                    genome_size = int(line_split[genome_size_index])
                    gc_percentage = float(line_split[gc_percentage_index])
                    contig_count = int(line_split[contig_count_index])
                    
                    protein_count = int(line_split[protein_count_index])
                    coding_density = float(line_split[coding_density_index])
                    
                    ssu_count = int(line_split[ssu_count_index])

                    metadata[gid] = Metadata(gtdb_representative=gtdb_representative,
                                                    gtdb_type_strain=gtdb_type_strain,
                                                    genome_quality=quality,
                                                    genome_comp=comp,
                                                    genome_cont=cont,
                                                    genome_size=genome_size,
                                                    gc_percentage=gc_percentage,
                                                    contig_count=contig_count,
                                                    protein_count=protein_count,
                                                    coding_density=coding_density,
                                                    ssu_count=ssu_count,
                                                    gtdb_taxonomy=gtdb_taxonomy,
                                                    ncbi_taxonomy=ncbi_taxonomy)

        return metadata
        
    def get_stat_data(self, 
                        metadata, 
                        field, 
                        scale_factor, 
                        min_value,
                        max_value,
                        all_genomes):
        """Get statistic across desired genomes."""
        
        data = []
        for gid, m in metadata.items():
            if all_genomes or m.gtdb_representative:
                v = getattr(m, field)/scale_factor
                
                if v < min_value or v > max_value:
                    continue
                
                data.append(v)
                
        return data
 
    def run(self, 
                bac120_metadata_file,
                ar120_metadata_file,
                all_genomes):
        """Plot of common genomic statistics."""

        # get genome metadata
        self.logger.info('Reading GTDB metadata.')
        metadata = self.read_metadata(bac120_metadata_file, ar120_metadata_file)
        
        # create plot for each genomic statistic
        options = AbstractPlot.Options(width=6.5, 
                                        height=4, 
                                        label_font_size=7, 
                                        tick_font_size=6, 
                                        dpi=600)
        plot = GenomicStatsPlot(options, 2, 3)
        
        if all_genomes:
            ylabel = 'Genomes'
            out_prefix = f'gtdb_r{self.release_number}_genomic_stats.genomes'
        else:
            ylabel = 'Species'
            out_prefix = f'gtdb_r{self.release_number}_genomic_stats.species'
            
        
        table_data = []
        table_data.append(('', 
                            'Median', 
                            'Mean', 
                            'Std. deviation',
                            'Min.',
                            'Max.',
                            '5th percentile',
                            '95th percentile'))
        
        # Genome size panel
        self.logger.info(f'Creating genome size plot.')
        data = self.get_stat_data(metadata, 'genome_size', 1e6, 0, 20, all_genomes)
        plot.plot(1, data, 'Genome size (Mb)', f'{ylabel} ({len(data):,})')
        table_data.append(('Genome size (Mb)', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))
        
        # GC-content panel
        self.logger.info(f'Creating GC content plot.')
        data = self.get_stat_data(metadata, 'gc_percentage', 1, 20, 80, all_genomes)
        plot.plot(2, data, 'GC content (%)', f'{ylabel} ({len(data):,})')
        table_data.append(('GC content (%)', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))
        
        # SSU count panel
        self.logger.info(f'Creating SSU count plot.')
        data = self.get_stat_data(metadata, 'ssu_count', 1, -1, 12, all_genomes)
        plot.plot(3, data, 'No. SSU genes', f'{ylabel} ({len(data):,})')
        table_data.append(('No. SSU genes', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))
        
        # CDS count panel
        self.logger.info(f'Creating CDS count plot.')
        data = self.get_stat_data(metadata, 'protein_count', 1, 100, 10000, all_genomes)
        plot.plot(4, data, 'No. CDSs', f'{ylabel} ({len(data):,})')
        table_data.append(('No. CDSs', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))
        
        # Coding density panel
        self.logger.info(f'Creating coding density plot.')
        data = self.get_stat_data(metadata, 'coding_density', 1, 70, 100, all_genomes)
        plot.plot(5, data, 'Coding density (%)', f'{ylabel} ({len(data):,})')
        table_data.append(('Coding density (%)', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))
        
        # Contig count panel
        self.logger.info(f'Creating contig count plot.')
        data = self.get_stat_data(metadata, 'contig_count', 1, 1, 1000, all_genomes)
        plot.plot(6, data, 'No. contigs', f'{ylabel} ({len(data):,})')
        table_data.append(('No. contigs', 
                            f'{np_median(data):.1f}',
                            f'{np_mean(data):.1f}',
                            f'{np_std(data):.1f}',
                            f'{min(data):.1f}',
                            f'{max(data):.1f}',
                            f'{np_percentile(data, 5):.1f}',
                            f'{np_percentile(data, 95):.1f}'))

        plot.save_plot(self.output_dir / f'{out_prefix}.png', dpi=600)
        plot.save_plot(self.output_dir / f'{out_prefix}.svg', dpi=600)
        
        # write out table
        fout = open(self.output_dir / f'{out_prefix}.tsv','w')
        for row in table_data:
            fout.write('{}\n'.format('\t'.join(table_data[0])))
        fout.close()
        
        fout = open(self.output_dir / f'{out_prefix}.html','w')
        htmlcode = HTML.table(table_data,
                                col_align=['left'] + ['center']*(len(table_data[0])-1),
                                col_styles=['font-size: small']*len(table_data[0]),
                                cellpadding=5)
        fout.write(htmlcode)
        fout.close()
