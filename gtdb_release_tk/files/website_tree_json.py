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
import json

from gtdb_release_tk.files.metadata import MetadataFile
from gtdb_release_tk.files.taxonomy import TaxonomyFile


class WebsiteTreeJsonFile(object):

    def __init__(self, data):
        self.data = data

    @classmethod
    def read(cls, path):
        with open(path) as f:
            data = json.loads(f.read())
        return cls(data)

    @classmethod
    def create(cls, taxonomy_file: TaxonomyFile, metadata_file: MetadataFile):
        tax_dict = dict()

        meta_information = dict()
        type_spe_list = list()

        for gid, tax in taxonomy_file.data.items():

            # Create metadata
            cur_meta = metadata_file.rows[gid]
            meta_information[cur_meta.accession] = {'rep': cur_meta.gtdb_representative,
                                                    'type': cur_meta.gtdb_type_designation,
                                                    'comp': cur_meta.checkm_completeness,
                                                    'conta': cur_meta.checkm_contamination,
                                                    'type_spe': cur_meta.gtdb_type_species_of_genus}

            if cur_meta.gtdb_type_species_of_genus == 't':
                type_spe_list.append(cur_meta.gtdb_taxonomy.s)

            # Create the tree
            if not tax.d.is_empty():
                if tax.d.rank not in tax_dict:
                    tax_dict[tax.d.rank] = {}
            else:
                tax_dict.set_default('d__Undefined', {}).set_default(gid, {})
                continue
            if not tax.p.is_empty():
                if tax.p.rank not in tax_dict[tax.d.rank]:
                    tax_dict[tax.d.rank][tax.p.rank] = {}
            if not tax.c.is_empty():
                if tax.c.rank not in tax_dict[tax.d.rank][tax.p.rank]:
                    tax_dict[tax.d.rank][tax.p.rank][tax.c.rank] = {}
            if not tax.o.is_empty():
                if tax.o.rank not in tax_dict[tax.d.rank][tax.p.rank][tax.c.rank]:
                    tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank] = {}
            if not tax.f.is_empty():
                if tax.f.rank not in tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank]:
                    tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank][tax.f.rank] = {}
            if not tax.g.is_empty():
                if tax.g.rank not in tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank][tax.f.rank]:
                    tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank][tax.f.rank][tax.g.rank] = {}
            if not tax.s.is_empty():
                if tax.s.rank not in tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank][tax.f.rank][tax.g.rank]:
                    tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank][tax.f.rank][tax.g.rank][tax.s.rank] = []
                tax_dict[tax.d.rank][tax.p.rank][tax.c.rank][tax.o.rank][tax.f.rank][tax.g.rank][tax.s.rank].append(gid)

        to_dump = {"name": "Domain", "type": "root", "children": list(WebsiteTreeJsonFile.conv(
            k, tax_dict[k], meta_information, type_spe_list) for k in sorted(tax_dict.keys()))}

        return cls(to_dump)

    def write(self, path):
        with open(path, 'w') as f:
            json.dump(self.data, f, check_circular=False)

    @staticmethod
    def conv(k, v, meta_information, type_spe_list):
        dict_rank = {"name": k,
                     "type": "d",
                     "children": list({"name": x,
                                       "type": "p",
                                       "children": list({"name": z,
                                                         "type": "c",
                                                         "children": list({"name": o,
                                                                           "type": "o",
                                                                           "children": list({"name": f,
                                                                                             "type": "f",
                                                                                             "children": list(
                                                                                                 {"name": g,
                                                                                                  "type": "g",
                                                                                                  "children": list(
                                                                                                      {"name": s,
                                                                                                       "type": "s",
                                                                                                       "children": list(
                                                                                                           WebsiteTreeJsonFile.process_genome(
                                                                                                               gen,
                                                                                                               meta_information)
                                                                                                           for gen in
                                                                                                           WebsiteTreeJsonFile.sorted_genome(
                                                                                                               v[x][z][
                                                                                                                   o][
                                                                                                                   f][
                                                                                                                   g][
                                                                                                                   s],
                                                                                                               meta_information))}
                                                                                                      for s in
                                                                                                      WebsiteTreeJsonFile.sort_ranks(
                                                                                                          v[x][z][o][f][
                                                                                                              g].keys()))}
                                                                                                 for g in
                                                                                                 WebsiteTreeJsonFile.sort_ranks(
                                                                                                     v[x][z][o][
                                                                                                         f].keys()))}
                                                                                            for f in
                                                                                            WebsiteTreeJsonFile.sort_ranks(
                                                                                                v[x][z][o].keys()))} for
                                                                          o in WebsiteTreeJsonFile.sort_ranks(
                                                             v[x][z].keys()))} for z in
                                                        WebsiteTreeJsonFile.sort_ranks(v[x].keys()))} for x in
                                      WebsiteTreeJsonFile.sort_ranks(v.keys()))}
        cleaned_dict = WebsiteTreeJsonFile.clean_empty_children(dict_rank)
        cleaned_dict = WebsiteTreeJsonFile.count_children(cleaned_dict)
        cleaned_dict = WebsiteTreeJsonFile.add_type_species(cleaned_dict, type_spe_list)
        return cleaned_dict

    @staticmethod
    def clean_empty_children(d):
        if "children" in d:
            v = d.get("children")
            if not v:
                del d["children"]
                return d
            else:
                tmp_list = []
                for it in v:
                    dict_tmp = WebsiteTreeJsonFile.clean_empty_children(it)
                    tmp_list.append(dict_tmp)
                d["children"] = tmp_list
                return d
        else:
            return d

    @staticmethod
    def count_children(d):
        if "extra_genomes" in d:
            d["countgen"] = d["extra_genomes"]
            return d
        if "children" in d and "value" not in d:
            v = d.get("children")
            countgen = 0
            for it in v:
                dict_tmp = WebsiteTreeJsonFile.count_children(it)
                countgen += dict_tmp["countgen"]
            d["countgen"] = countgen
            return d

        elif "children" not in d and "value" not in d:
            d["countgen"] = 1
            return d
        return d

    @staticmethod
    def add_type_species(d, type_spe_list):
        if "children" in d:
            v = d.get("children")
            for it in v:
                if it['name'] in type_spe_list:
                    it['type_material'] = 'type_species'
                else:
                    dict_tmp = WebsiteTreeJsonFile.add_type_species(it, type_spe_list)
        return d

    @staticmethod
    def sort_ranks(list_ranks):
        valid_names = sorted([x for x in list_ranks if WebsiteTreeJsonFile.is_valid_name(x)])
        placeholder_names = sorted(
            [x for x in list_ranks if x not in valid_names])
        results = valid_names + placeholder_names
        return results

    @staticmethod
    def is_valid_name(rank_name):
        rank_name = rank_name[3:].split('_')[0]
        if rank_name.replace(' ', '').isalpha() and rank_name[0].isupper() and rank_name[1:].islower():
            return True
        return False

    @staticmethod
    def process_genome(genome, meta_information):
        result_genome = {"name": genome, "type": 'genome'}
        if genome.startswith('+'):
            result_genome["extra_genomes"] = int(
                genome.split(' ')[1].replace(',', ''))
            return result_genome
        if meta_information.get(genome).get('type').startswith('type strain'):
            result_genome["type_material"] = meta_information.get(
                genome).get('type').replace(' ', '_')
        if meta_information.get(genome).get('rep') == 't':
            result_genome["rep"] = True
        return result_genome

    @staticmethod
    def sorted_genome(list_genomes, meta_information):
        list_representatives = [
            gen for gen in list_genomes if meta_information.get(gen).get('rep') == 't']
        list_type_species = [gen for gen in list_genomes if meta_information.get(gen).get(
            'type') == 'type strain of species' and gen not in list_representatives]
        list_type_subspecies = [gen for gen in list_genomes if meta_information.get(gen).get(
            'type') == 'type strain of subspecies' and gen not in list_representatives]

        dict_quality = {gen: meta_information.get(gen).get('comp') - 5 * meta_information.get(gen).get(
            'conta') for gen in list_genomes if
                        gen not in list_representatives + list_type_species + list_type_subspecies}
        list_quality = sorted(dict_quality, key=dict_quality.get)
        full_list = list_representatives + list_type_species + \
                    list_type_subspecies + list_quality
        if len(full_list) > 100:
            shorten_list = full_list[0:100]
            remaining_genomes = len(full_list) - len(shorten_list)
            shorten_list.append(f'+ {remaining_genomes:,} additional genomes.')
            return shorten_list
        return full_list
