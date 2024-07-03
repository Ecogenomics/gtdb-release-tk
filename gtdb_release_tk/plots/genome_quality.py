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
from pathlib import PurePath
from collections import namedtuple

from biolib.plots.abstract_plot import AbstractPlot
                   
from matplotlib.ticker import NullFormatter
import matplotlib.lines as mlines

from gtdb_release_tk.plots.palette import DEFAULT_PALETTE


class GenomeQualityPlot(AbstractPlot):
    """Create histogram of common genomic statistics."""

    def __init__(self, options):
        """Initialize."""
        
        AbstractPlot.__init__(self, options)
        
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

    def plot(self, comp, cont, 
                    mimag_category, exception_count,
                    xlabel, ylabel, 
                    num_bins, 
                    xlim=None, 
                    ylim=None,
                    palette=DEFAULT_PALETTE):
        """Create plot."""
        
        # setup histograms
        self.histogram_size = 0.5
        histogram_sizeX = self.histogram_size / self.options.width
        histogram_sizeY = self.histogram_size / self.options.height

        padding = 0.1   # inches
        xoffset_figspace = (0.4 + padding)/self.options.width
        yoffset_figspace = (0.3 + padding)/self.options.height
        axes_scatter = self.fig.add_axes([xoffset_figspace, 
                                            yoffset_figspace,
                                            1.0 - xoffset_figspace - histogram_sizeX - (2*padding)/self.options.width, 
                                            1.0 - yoffset_figspace - histogram_sizeY - (2*padding)/self.options.height])

        axes_top_histogram = self.fig.add_axes([xoffset_figspace, 
                                                    1.0 - histogram_sizeY - padding/self.options.height,
                                                    1.0 - xoffset_figspace - histogram_sizeX - (2*padding)/self.options.width, 
                                                    histogram_sizeY])

        axes_right_histogram = self.fig.add_axes([1.0 - histogram_sizeX - padding/self.options.width, 
                                                    yoffset_figspace,
                                                    histogram_sizeX, 
                                                    1.0 - yoffset_figspace - histogram_sizeY - (2*padding)/self.options.height])

        # scatter plot
        hq_x = []
        hq_y = []
        hq_c = []
        mq_x = []
        mq_y = []
        mq_c = []
        p_x = []
        p_y = []
        p_c = []

        for cur_comp, cur_cont, category in zip(comp, cont, mimag_category):
            if category == 'hq':
                hq_x.append(cur_comp)
                hq_y.append(cur_cont)
                hq_c.append(palette.colour3)
            elif category == 'mq':
                mq_x.append(cur_comp)
                mq_y.append(cur_cont)
                mq_c.append(palette.colour2)
            elif category == 'lq':
                p_x.append(cur_comp)
                p_y.append(cur_cont)
                p_c.append(palette.colour1)
            else:
                assert(False)

        mq = axes_scatter.scatter(mq_x, mq_y, c=mq_c, s=15, lw=0.5, edgecolor=(0.2,0.2,0.2))
        hq = axes_scatter.scatter(hq_x, hq_y, c=hq_c, s=15, lw=0.5, edgecolor=(0.2,0.2,0.2))
        p = axes_scatter.scatter(p_x, p_y, c=p_c, s=15, lw=0.5, edgecolor=(0.2,0.2,0.2))
        
        hq_count = len(hq_x)
        mq_count = len(mq_x)
        p_count = len(p_x)
        
        total = hq_count + mq_count + p_count + exception_count
        hq_perc = hq_count * 100.0 / total
        mq_perc = mq_count * 100.0 / total
        p_perc = p_count * 100.0 / total
        exception_perc = exception_count * 100.0 / total
        
        red_square = mlines.Line2D([], [], color='black', marker='x', linestyle='None', markersize=3)
        lgnd = axes_scatter.legend((hq, mq, p, red_square),
                                   ('High quality ({:,} genomes; {:.1f}%)'.format(hq_count, hq_perc), 
                                   'Medium quality ({:,} genomes; {:.1f}%)'.format(mq_count, mq_perc),
                                   'Low quality ({:,} genomes; {:.1f}%)'.format(p_count, p_perc),
                                   'Exceptions ({:,} genomes; {:.1f}%)'.format(exception_count, exception_perc)),
                                   scatterpoints=1,
                                   loc='upper left',
                                   ncol=1,
                                   fontsize=self.options.tick_font_size,
                                   handletextpad=0,
                                   markerscale=1.5,
                                   frameon=False)

        axes_scatter.set_xlabel(xlabel)
        axes_scatter.set_ylabel(ylabel)
        
        if xlim:
            axes_scatter.set_xlim(xlim)
        if ylim:
            axes_scatter.set_ylim(ylim)
        
        # *** Prettify scatter plot
        for line in axes_scatter.yaxis.get_ticklines(): 
            line.set_color(self.axes_colour)
                
        for line in axes_scatter.xaxis.get_ticklines(): 
            line.set_color(self.axes_colour)
            
        for loc, spine in axes_scatter.spines.items():
            spine.set_color(self.axes_colour)

        # get bin weights for percentage histogram plot
        x = [p_x, mq_x, hq_x]
        num_entries = sum([len(c) for c in x])
        weights = []
        for c in range(0, len(x)):
            w = [100.0/num_entries for _ in range(0, len(x[c]))]
            weights.append(w)
    
        # plot top histogram
        axes_top_histogram.xaxis.set_major_formatter(NullFormatter())
        pdf, _bins, _patches = axes_top_histogram.hist(x, 
                                                    bins=num_bins, 
                                                    weights=weights,
                                                    color=[palette.colour1, palette.colour2, palette.colour3],
                                                    edgecolor=(0.2,0.2,0.2),
                                                    lw=0.5, 
                                                    histtype='bar',
                                                    rwidth=0.9,                                                    
                                                    stacked=True)
                                                                                       
        max_y = max(pdf[len(x)-1])
        axes_top_histogram.yaxis.set_label_position('right')
        axes_top_histogram.set_xlim(axes_scatter.get_xlim())
        axes_top_histogram.set_yticks([int(max_y)])
        axes_top_histogram.set_yticklabels(['%d%%' % int(max_y)])
        axes_top_histogram.set_ylim([0, max_y*1.05])
        
        axes_top_histogram.yaxis.tick_right()
        

        # plot right histogram
        y = [p_y, mq_y, hq_y]
        axes_right_histogram.yaxis.set_major_formatter(NullFormatter())
        pdf, _bins, _patches = axes_right_histogram.hist(y, 
                                                        bins=num_bins, 
                                                        orientation='horizontal',
                                                        weights=weights,
                                                        color=[palette.colour1, palette.colour2, palette.colour3],
                                                        edgecolor=(0.2,0.2,0.2),
                                                        lw=0.5, 
                                                        histtype='bar',
                                                        rwidth=0.9,
                                                        stacked=True)
                                                        
        max_x = max(pdf[len(y)-1])
        axes_right_histogram.set_ylim(axes_scatter.get_ylim())
        axes_right_histogram.set_xticks([int(max_x)])
        axes_right_histogram.set_xticklabels(['%d%%' % int(max_x)])
        axes_right_histogram.set_xlim([0, max_x*1.05])

        # *** Prettify histogram plot
        for a in axes_top_histogram.yaxis.majorTicks:
                a.tick1On=False
                a.tick2On=True
            
        for a in axes_top_histogram.xaxis.majorTicks:
                a.tick1On=True
                a.tick2On=False
                
        for line in axes_top_histogram.yaxis.get_ticklines(): 
            line.set_color(self.axes_colour)
        
        for line in axes_top_histogram.xaxis.get_ticklines(): 
            line.set_color(self.axes_colour)

        for loc, spine in axes_top_histogram.spines.items():
                if loc in ['left','top']:
                        spine.set_color('none')
                else:
                    spine.set_color(self.axes_colour)

        for a in axes_right_histogram.yaxis.majorTicks:
                a.tick1On=True
                a.tick2On=False
            
        for a in axes_right_histogram.xaxis.majorTicks:
                a.tick1On=True
                a.tick2On=False
                
        for line in axes_right_histogram.yaxis.get_ticklines(): 
            line.set_color(self.axes_colour)
        
        for line in axes_right_histogram.xaxis.get_ticklines(): 
            line.set_color(self.axes_colour)

        for loc, spine in axes_right_histogram.spines.items():
                if loc in ['right','top']:
                        spine.set_color('none') 
                else:
                    spine.set_color(self.axes_colour)

        self.draw()
   
class GenomeQuality(object):
    """Scatter plot showing quality of GTDB representative genomes."""

    def __init__(self, release_number, output_dir):
        
        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')
        
    def read_metadata(self, bac120_metadata_file, ar120_metadata_file):
        """Read metadata for genomes."""
        
        Metadata = namedtuple('Metadata', ['genome_comp',
                                            'genome_cont',
                                            'mimag_hq',
                                            'mimag_mq',
                                            'mimag_lq'])

        metadata = {}
        mimag_hq_count = 0
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')
                
                genome_index = header.index('accession')
                comp_index = header.index('checkm_completeness')
                cont_index = header.index('checkm_contamination')
                hq_index = header.index('mimag_high_quality')
                mq_index = header.index('mimag_medium_quality')
                lq_index = header.index('mimag_low_quality')
                gtdb_rep_index = header.index('gtdb_representative')
                

                for line in f:
                    line_split = line.strip().split('\t')
                    
                    if line_split[gtdb_rep_index] != 't':
                        continue
                    
                    gid = line_split[genome_index]

                    comp = float(line_split[comp_index])
                    cont = float(line_split[cont_index])
                    
                    mimag_hq = line_split[hq_index] == 't'
                    mimag_mq = line_split[mq_index] == 't'
                    mimag_lq = line_split[lq_index] == 't'

                    metadata[gid] = Metadata(genome_comp=comp,
                                                genome_cont=cont,
                                                mimag_hq=mimag_hq,
                                                mimag_mq=mimag_mq,
                                                mimag_lq=mimag_lq)

        return metadata
 
    def run(self, 
                bac120_metadata_file,
                ar120_metadata_file,
                palette):
        """Scatter plot showing quality of GTDB representative genomes."""

        # get genome metadata
        self.logger.info('Reading GTDB metadata.')
        metadata = self.read_metadata(bac120_metadata_file, ar120_metadata_file)
        self.logger.info(f' - read metadata for {len(metadata):,} representative genomes.')
        
        # get completenss, contamination, and MIMAG quality of each 
        # GTDB species representative
        comp = []
        cont = []
        mimag_category = []
        exception_count = 0
        for gid, m in metadata.items(): 
            if m.mimag_hq:
                mimag_category.append('hq')
                comp.append(m.genome_comp)
                cont.append(m.genome_cont)
            elif m.mimag_mq:
                mimag_category.append('mq')
                comp.append(m.genome_comp)
                cont.append(m.genome_cont)
            elif m.mimag_lq:
                mimag_category.append('lq')
                comp.append(m.genome_comp)
                cont.append(m.genome_cont)
            else:
                #***self.logger.error(f'GTDB representative has no assigned MIMAG quality: {gid}')
                # this occurs becomes some representatives have >10% contamination
                exception_count += 1
                
        self.logger.info(' - HQ = {:,}, MQ = {:,}, LQ = {:,}, exceptions = {:,}'.format(
                            mimag_category.count('hq'),
                            mimag_category.count('mq'),
                            mimag_category.count('lq'),
                            exception_count))
        
        # create plot for each genomic statistic
        options = AbstractPlot.Options(width=4, 
                                        height=4, 
                                        label_font_size=7, 
                                        tick_font_size=6, 
                                        dpi=600)
        plot = GenomeQualityPlot(options)
        
        plot.plot(comp, cont, 
                    mimag_category,
                    exception_count,
                    'Completeness (%)', 
                    'Contamination (%)', 
                    num_bins=25,
                    xlim=(49, 101), 
                    ylim=(-0.2, 10.2),
                    palette=palette)

        out_prefix = f'gtdb_r{self.release_number}_genome_quality.species'
        plot.save_plot(self.output_dir / f'{out_prefix}.png', dpi=600)
        plot.save_plot(self.output_dir / f'{out_prefix}.svg', dpi=600)
