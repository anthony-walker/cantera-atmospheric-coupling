import sys
import ruamel.yaml
from cac.mechstats import *
from cac.constants import DATA_DIR

all_mechanisms()

# yaml = ruamel.yaml.YAML(typ='safe', pure=True)
# yaml.preserve_quotes = True
# yaml.sort_keys = False

# with open("scombustor.yaml", "r") as f:
#     ydata = yaml.load(f)
# photo = []
# nonphoto = []
# for r in ydata["reactions"]:
#     if r.get("type", "") == "zenith-angle-rate":
#         photo.append(r)
#     else:
#         nonphoto.append(r)

# for pr in photo:
#     for r in nonphoto:
#         if r["equation"] == pr["equation"]:
#             pr["duplicate"] = True
#             r.pop("duplicate", "")

# ydata["reactions"] = nonphoto
# ydata["photolysis-reactions"] = photo

# with open("scombustor.yaml", "w") as f:
#     yaml.dump(ydata, f)
