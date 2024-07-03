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
from pathlib import PurePath

from biolib.taxonomy import Taxonomy
from biolib.plots.abstract_plot import AbstractPlot

from numpy import (arange as np_arange,
                    array as np_array)

from gtdb_release_tk.plots.palette import COLOR_BLIND_PALETTE


class NCBI_ComparePlot(AbstractPlot):
    """Box and whisker plot."""

    def __init__(self, options):
        """Initialize."""
        AbstractPlot.__init__(self, options)

    def plot(self, unchanged, passive, active, xticklabels, ylabel, palette):
        """Create stacked bar plot."""
        
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axis = self.fig.add_subplot(111)
        
        ind = np_arange(len(unchanged))
        width = 0.7
        
        unchanged = np_array(unchanged)
        passive = np_array(passive)
        active = np_array(active)

        p1 = axis.bar(ind, unchanged, width, color=palette.colour1)
        p2 = axis.bar(ind, passive, width, bottom=unchanged, color=palette.colour2)
        p3 = axis.bar(ind, active, width, bottom=unchanged+passive, color=palette.colour3)

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
        
        axis.legend((p3[0], p2[0], p1[0]), ('Active change', 'Passive change', 'Unchanged'),
                        fontsize=self.options.tick_font_size,
                        loc='upper left', 
                        bbox_to_anchor=(1, 1),
                        frameon=False)

        #self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()


class NCBI_Compare(object):
    """Bar plot comparing GTDB and NCBI taxonomies."""
    
    def __init__(self, release_number, output_dir):
        """Initialization."""
        
        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)
        
        self.logger = logging.getLogger('timestamp')
        
    def run(self, 
                bac120_metadata_file, 
                ar120_metadata_file,
                all_genomes,
                domain,
                palette):
        """Bar plot comparing GTDB and NCBI taxonomies."""
        
        # parse GTDB metadata file to determine genomes in each species clusters
        self.logger.info('Reading GTDB metadata.')
        gtdb_taxonomy = {}
        ncbi_taxonomy = {}
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')
                
                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                ncbi_taxonomy_index = header.index('ncbi_taxonomy')
                gtdb_rep_index = header.index('gtdb_representative')

                for line in f:
                    line_split = line.strip().split('\t')
                    
                    gid = line_split[0]

                    gtdb_rep = line_split[gtdb_rep_index]
                    if not (all_genomes or gtdb_rep == 't'):
                        continue
                        
                    taxonomy = line_split[gtdb_taxonomy_index]
                    gtdb_taxa = [t.strip() for t in taxonomy.split(';')]
                    gtdb_domain = gtdb_taxa[0]
                    if domain == 'Both' or domain in gtdb_domain:
                        gtdb_taxonomy[gid] = gtdb_taxa
                        
                        taxonomy = line_split[ncbi_taxonomy_index]
                        if taxonomy != 'none':
                            ncbi_taxa = [t.strip() for t in taxonomy.split(';')]
                            ncbi_taxonomy[gid] = ncbi_taxa

        self.logger.info(' ...identified {:,} GTDB genomes and {:,} NCBI genomes.'.format(
                            len(gtdb_taxonomy),
                            len(ncbi_taxonomy)))
                            
        # compare NCBI and GTDB taxa at each rank
        out_prefix = f'gtdb_r{self.release_number}_ncbi_compare'
        if all_genomes:
            out_prefix += '.genomes'
        else:
            out_prefix += '.species'
            
        if domain == 'Bacteria':
            out_prefix += '.bacteria'
        elif domain == 'Archaea':
            out_prefix += '.archaea'
            
        self.logger.info('Comparing NCBI and GTDB taxa.')
        fout = open(self.output_dir / f'{out_prefix}.summary.tsv', 'w')
        fout.write('Rank\tUnchanged\tPassive change\tActive change\n')
        
        fout_taxon = open(self.output_dir / f'{out_prefix}.taxon.tsv', 'w')
        fout_taxon.write('Genome ID\tNCBI taxon\tGTDB taxon\tClassification\n')

        plot_unchanged = []
        plot_passive = []
        plot_active = []
        plot_labels = []
        active_change = set()
        active_change_no_species = set()
        active_or_passive_change = set()
        for rank_index in range(1, 7):
            fout.write(Taxonomy.rank_labels[rank_index])
            
            unchanged = 0
            passive = 0
            active = 0
            for gid, ncbi_taxa in ncbi_taxonomy.items():
                ncbi_taxon = ncbi_taxa[rank_index].replace('Candidatus ', '')
                gtdb_taxon = gtdb_taxonomy[gid][rank_index]
                classification = None
                if ncbi_taxon == gtdb_taxon:
                    unchanged += 1
                    classification = 'unchanged'
                elif ncbi_taxon == Taxonomy.rank_prefixes[rank_index]:
                    passive += 1
                    classification = 'passive change'
                    active_or_passive_change.add(gid)
                else:
                    active += 1
                    classification = 'active change'
                    active_change.add(gid)
                    active_or_passive_change.add(gid)
                    
                    if rank_index != 6:
                        active_change_no_species.add(gid)
                    
                fout_taxon.write(f'{gid}\t{ncbi_taxon}\t{gtdb_taxon}\t{classification}\n')

            total_taxa = unchanged + passive + active
            fout.write(f'\t{total_taxa}\t{unchanged}\t{active}\n')
            
            plot_unchanged.append(unchanged*100.0/total_taxa)
            plot_passive.append(passive*100.0/total_taxa)
            plot_active.append(active*100.0/total_taxa)
            plot_labels.append('{}'.format(
                                Taxonomy.rank_labels[rank_index].capitalize()))
            
            if all_genomes:
                ylabel = 'Genomes (%)\n({:,})'.format(total_taxa)
            else:
                ylabel = 'Representative Genomes (%)\n({:,})'.format(total_taxa)
            
        fout.close()
        fout_taxon.close()
        
        self.logger.info(f'Identified {len(active_change_no_species):,} genomes with one or more active changes above the rank of species.')
        self.logger.info(f'Identified {len(active_change):,} genomes with one or more active changes.')
        self.logger.info(f'Identified {len(active_or_passive_change):,} genomes with one or more active or passive changes.')
        
        # create plot
        self.logger.info('Creating plot.')
        options = AbstractPlot.Options(width=4, 
                                        height=3, 
                                        label_font_size=7, 
                                        tick_font_size=6, 
                                        dpi=600)
        plot = NCBI_ComparePlot(options)
        plot.plot(plot_unchanged, plot_passive, plot_active, plot_labels, ylabel, palette)
        
        plot.save_plot(self.output_dir / f'{out_prefix}.png', dpi=600)
        plot.save_plot(self.output_dir / f'{out_prefix}.svg', dpi=600)