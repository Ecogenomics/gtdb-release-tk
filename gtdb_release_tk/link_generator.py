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

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2019'
__credits__ = ['Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import os

import pandas as pd


class LinkGenerator(object):

    def __init__(self, release_number, output_dir):
        self.release_number = release_number
        self.output_dir = output_dir

    def generate_sandpiper_links(self, web_table_path):
        # Determine the file extension and read accordingly
        file_ext = os.path.splitext(web_table_path)[-1].lower()
        sep = '\t' if file_ext == '.tsv' else ','

        df = pd.read_csv(web_table_path, sep=sep)

        # Filter out rows where 'type' is 'genome'
        df = df[df['type'] != 'genome']
        df = df[df['id'] != 0]

        # Replace spaces in taxon names and create URLs
        df['url'] = df['taxon'].str.replace(" ", "%20").apply(
            lambda taxon: f"https://sandpiper.qut.edu.au/taxonomy/{taxon}"
        )

        # Keep only rows with URLs
        df = df[['id', 'url']].dropna()

        output_path = os.path.join(
            self.output_dir, f'sandpiper_links_{self.release_number}.tsv'
        )

        df.to_csv(output_path, index=False, sep='\t')

        print(f"Sandpiper links generated for release {self.release_number} and saved to {output_path}")





