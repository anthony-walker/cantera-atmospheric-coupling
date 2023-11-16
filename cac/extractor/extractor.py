import re
import os
import time
import click
import string
import random
import shutil
import datetime
import ruamel.yaml
from cac.constants import DATA_DIR
from cac.extractor.section_extractor import extract_from_fac
from cac.extractor.species_extractor import write_species_extraction
from cac.extractor.reaction_extractor import write_balanced_reaction_list

yaml = ruamel.yaml.YAML(typ='safe')
default_flow_style=False

def regenerate_files(facfile, dirname):
    prefix = facfile.split(".")[0]
    extract_from_fac(facfile, dirname)
    write_species_extraction(prefix, dirname)
    write_balanced_reaction_list(prefix, dirname)


@click.command()
@click.argument('facfile', nargs=1)
@click.option('--regenerate', default=True, help='Regenerate all needed files.')
def main(facfile, regenerate=True):
    """ This file is the main file used to construct the completed mechanism.

        To run "python species_extractor.py" and "python reaction_extractor.py".
        This will generate mcm-species.yaml and mcm-reactions.yaml. In between these
        steps a few species may need to be corrected, composition, etc. Specifically I
        know that NA, SA, and HO2NO2 will all need to be corrected.

        The run this file to construct the completed mechanism.

        Thermodynamic data is assumed to be ideal gas if it cannot be found. It
        searches through all the models in the `all-models` folder in an attempt to
        find it. Using InCHI and SMILES alias as well in this search.

        Some of the rates in these reactions are quite complex and custom which
        involves ExtensibleRateData in cantera.
    """
    prefix = facfile.split(".")[0]
    prefix_dir = os.path.join(DATA_DIR, prefix)
    if regenerate or not os.path.exists(prefix_dir):
        try:
            os.mkdir(prefix_dir)
        except Exception as e:
            pass
        regenerate_files(facfile, prefix_dir)
    # open template file
    tfile = os.path.join(DATA_DIR, "template.yaml")
    with open(tfile, "r") as f:
        aerosol_data = yaml.load(f)
        aerosol_data["date"] = datetime.datetime.now()
    # open mcm-species
    sfile = os.path.join(prefix_dir, f"{prefix}-species.yaml")
    with open(sfile, "r") as f:
        species_data = yaml.load(f)
    # open reactions
    rfile = os.path.join(prefix_dir, f"{prefix}-reactions.yaml")
    with open(rfile, "r") as f:
        reaction_data = yaml.load(f)
    # get all species names
    species_names, species_items = zip(*species_data.items())
    # get all elements
    complete_comp = {}
    for si in species_items:
        complete_comp.update(si["composition"])
    elements = list(complete_comp.keys())
    elements.sort()
    # atmosphere phase modification TODO: FIX elements
    aerosol_data["phases"][1]["species"] = species_names
    aerosol_data["phases"][1]["elements"] = elements
    # aerosol phase modification TODO: ALL
    aerosol_data["phases"][2]["species"] = ["NA", "SA"]
    aerosol_data["phases"][2]["elements"] = ["H", "N", "O", "S"]
    # transfer data
    aerosol_data["species"] = species_items
    aerosol_data["atmosphere-reactions"] = reaction_data["atmosphere-reactions"]
    aerosol_data["aerosol-reactions"] = reaction_data["aerosol-reactions"]
    with open(f"{prefix}.yaml", "w") as f:
        yaml.dump(aerosol_data, f, default_flow_style=False, sort_keys=False)
    # Format unruly lists like elements and species in phases to make more readable.
    with open(f"{prefix}.yaml", "r") as f:
        content = f.read()
        content = re.sub(r"[']NO[']", "NO", content)
        content = re.sub(r"\n  [-] (([A-Z]+[a-z]*[0-9]*)+)", r", \1", content)
        content = re.sub(r"[:][,]", ":", content)
        content = re.sub(r"[']([-]?\d+([.]\d*e?[+-]?\d*)?)[']", r"\1", content)
        # search for list
        random.seed(time.time())
        rep_str = "".join(random.choices(string.ascii_letters, k=30))
        res = re.search(r"(([A-Z]+[a-z]*[0-9]*)+[,][ ])+[^\n]*", content)
        replacements = []
        while res:
            replacements.append((rep_str, res.group(0)))
            content = content.replace(res.group(0), rep_str)
            rep_str = "".join(random.choices(string.ascii_letters, k=30))
            res = re.search(r"(([A-Z]+[a-z]*[0-9]*)+[,][ ])+[^\n]*", content)
        for rep, orig in replacements:
            orig = f"[{orig}]"
            content = content.replace(rep, orig)
    # rewrite content
    with open(f"{prefix}.yaml", "w") as f:
        f.write(content)
    pyprefix = prefix.replace("-", "_")
    relevant = list(filter(lambda x: prefix in x or pyprefix in x, os.listdir(prefix_dir)))
    try:
        os.mkdir(prefix)
    except Exception as e:
        pass
    for f in relevant:
        shutil.copyfile(os.path.join(prefix_dir, f), os.path.join(prefix, f))


if __name__ == "__main__":
    main()
