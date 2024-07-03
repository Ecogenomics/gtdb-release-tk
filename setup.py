import os

from setuptools import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'gtdb_release_tk', 'VERSION'), 'r') as f:
        return f.readline().strip()


setup(
    name='gtdb_release_tk',
    python_requires='>=3.6',
    version=version(),
    packages=['gtdb_release_tk', 'gtdb_release_tk.plots', 'gtdb_release_tk.tables'],
    scripts=['bin/gtdb_release_tk'],
    package_data={'gtdb_release_tk': ['VERSION']},
    url='https://github.com/Ecogenomics/gtdb-release-tk',
    license='GPL3',
    description='The GTDB Release Toolkit provides functionality for updating the GTDB to the next release and generating data files for the GTDB website.',
    install_requires=['biolib>=0.1.0', 'dendropy', 'unidecode', 'requests', 'matplotlib',
                      'numpy', 'tqdm', 'pandas', 'seaborn', 'gtdblib>=1.1.0,<1.2.0'],
)
