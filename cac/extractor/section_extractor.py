import re
import os
from cac.constants import DATA_DIR
from sympy.parsing.sympy_parser import parse_expr


def extract_from_fac(facfile, dirname, replace_generics=True, replace_complex=True):
    prefix = facfile.split(".")[0]
    # break down file into sections
    with open(facfile, "r") as f:
        content = f.read()
        delimiters = re.findall(r"[*][;]\n[*].*[;][\n][*][;]\n", content)
        sections = []
        for d in delimiters[::-1]:
            content, app = content.split(d)
            sections.append(app)
        sections = sections[::-1]
    # write species file
    species = re.findall(r"(?:[A-Z]+[a-z]*[0-9]*)+", sections.pop(0))
    species.remove("VARIABLE")
    sfile = os.path.join(dirname, f"{prefix}-species.txt")
    with open(sfile, "w") as f:
        for sp in species:
            f.write(f"{sp.strip()}\n")
    # substitute in complex rate coefficients where possible
    ccs = {}
    complex_coeffs = re.sub(r"\s*", "", sections.pop(0))
    complex_coeffs = re.sub(r"[*]+[;]", "", complex_coeffs)
    for cc in complex_coeffs.split(";")[:-1]:
        cc_key, cc_rate = cc.split("=")
        cc_rate = re.sub(r"EXP", r"exp", cc_rate)
        cc_rate = re.sub(r"(\d+)[DE]([+-]?)(\d+)", r"\1e\2\3", cc_rate)
        cc_rate = re.sub(r"[@]([-]?\d+([.]\d+)?([e][+-]?\d+)?)", r"**(\1)", cc_rate)
        ccs[cc_key] = cc_rate
    # get the simplified complex expressions
    ccs_keys = ccs.keys()
    change = True
    temp_data = {x: y for x, y in ccs.items()}
    while change:
        change = False
        crates = sorted(temp_data.items(), key=lambda x: len(x[0]), reverse=True)
        for rname, rexpr in crates:
                found = re.findall(r"(?:[A-Z]+[0-9]*)+", rexpr)
                found.sort(key=lambda x: len(x), reverse=True)
                for fexp in found:
                    if fexp in ccs_keys:
                        rexpr = rexpr.replace(fexp, ccs[fexp])
                        change = True
                temp_data[rname] = rexpr
    for rname, rate in temp_data.items():
        rate = re.sub(r"EXP", r"exp", rate)
        rate = re.sub(r"(\d+)[DE]([+-]?)(\d+)", r"\1e\2\3", rate)
        rate = re.sub(r"[@]", r"**", rate.strip())
        ccs[rname] = str(parse_expr(rate))
    # write ro2 species file
    species = re.findall(r"(?:[A-Z]+[a-z]*[0-9]*)+", sections.pop(0).split("=")[1])
    ro2_file = os.path.join(dirname, f"{prefix}-ro2-sum.txt")
    with open(ro2_file, "w") as f:
        for sp in species:
            f.write(f"{sp.strip()}\n")
    # separate reactions and rates
    equations = re.findall(r"[%].*[;]", sections.pop(0))
    rrs = []
    for eq in equations:
        ceq = re.sub(r"[%](.*)[;]", r"\1", eq)
        rate, reaction = ceq.split(":")
        rate = re.sub(r"EXP", r"exp", rate)
        rate = re.sub(r"(\d+)[DE]([+-]?)(\d+)", r"\1e\2\3", rate)
        rate = re.sub(r"[@]([-]?\d+([.]\d+)?([e][+-]?\d+)?)", r"**(\1)", rate.strip())
        # check for the generic rate
        rrs.append((reaction.strip(), rate))
    rrs.sort()
    reactions, rates = zip(*rrs)

    if replace_complex:
        arrhen_regex = r"\d+([.]\d*)?([e][+-]\d+)?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?([*][(]TEMP[/]\d+[)][*][*][-]?\d+([.]\d*)?)?([*]exp[(][-]?\d+[/]TEMP[)])?(([*]\d+([.]\d*)?([e][+-]\d+)?)+)?"
        rate_str = "\n".join(rates)
        sort_cc = list(ccs.items())
        sort_cc.sort(key=lambda x: len(x[0]), reverse=True)
        for cc, cc_exp in sort_cc:
            if re.fullmatch(arrhen_regex,cc_exp):
                rate_str = rate_str.replace(cc, cc_exp)
        rates = rate_str.split("\n")
    # write reactions file
    react_file = os.path.join(dirname, f"{prefix}-reactions.txt")
    with open(react_file, "w") as f:
        f.write("\n".join(reactions))
    # write rates file
    rate_file = os.path.join(dirname, f"{prefix}-rates.txt")
    with open(rate_file, "w") as f:
        f.write("\n".join(rates))


if __name__ == "__main__":
    extract_from_fac("full-mcm.fac", replace_generics=True, replace_complex=False)
