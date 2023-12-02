import os
import re
import numpy as np
import cantera as ct
from cac.constants import DATA_DIR
from sympy.parsing.sympy_parser import parse_expr


def get_photolysis_parameterization(photo_string):
    """ Get photolysis parameterization yaml

    Args:
        photo_string (str): J<n> where n is the number associated with
        photolysis.txt
    """
    res = re.search(r"[J][<]\d+[>]", photo_string)
    assert res
    J_string = res.group(0)
    # calculate scalar if there is one
    photo_string = re.sub(r"[J][<]\d+[>]", "", photo_string)
    scalars = re.findall(r"(\d+([.]\d*)?([e][+-]\d+)?)+", photo_string)
    scalars = [float(s[0]) for s in scalars]
    scalar = np.prod(scalars) if scalars else 1
    # reassign photo string as J string
    photo_string = J_string
    with open(os.path.join(DATA_DIR, "photolysis.txt"), "r") as f:
        content = f.read()
    content = re.sub(r"[ \t\r\f\v]+", " ", content)
    lines = [l.strip() for l in content.split("\n")]
    number = photo_string.strip()[2:-1]
    lines = list(filter(lambda x: x.split(" ")[0] == number, lines))
    assert len(lines) == 1
    data = lines[0].split(" ")
    data = [re.sub(r"D", "e", d) for d in data]
    # add conversion into scalar to convert rate to m^3 / kmol /s
    scalar *= (ct.avogadro / 1e6)
    J, l, m, n = data
    return {"type": "zenith-angle-rate", "l": l, "m": m, "n": n, "scalar":str(scalar)}


def get_list_of_rate_data(prefix, dirname):
    """ Get a complete list of rate data

    Returns:
        list: data corresponding to each reaction
    """
    rate_file = os.path.join(dirname, f"{prefix}-rates.txt")
    generic_file = os.path.join(dirname, f"{prefix}-generic.txt")
    with open(rate_file, "r") as f:
        rates = f.read().split("\n")
    # get all generic rates
    with open(generic_file, "r") as f:
        gen_rates = f.read().split("\n")[:-1]
    gen_dict = {}
    for gr in gen_rates:
        grn, gexp = gr.split(":")
        grn = grn.strip()
        gexp = gexp.strip()
        gen_dict[grn] = gexp
    modified = True
    while modified:
        modified = False
        for n0, r0 in gen_dict.items():
            for n1, r1 in gen_dict.items():
                if re.search(f"(^|[^A-Z0-9]){n1}([^A-Z0-9]|$)", r0):
                    nexp = re.sub(f"(^|[^A-Z0-9]){n1}([^A-Z0-9]|$)", r"\1_REPLACE_ME_\2", r0)
                    nexp = re.sub("_REPLACE_ME_", r1, nexp)
                    gen_dict[n0] = nexp
                    modified = True
    for gn, gr in gen_dict.items():
        expr = re.sub("TEMP", "T", gr)
        expr = re.sub("[@]", "**", expr)
        expr = re.sub("[D]", "E", expr)
        expr = re.sub("LOG10", "log10", expr)
        expr = re.sub("EXP", "exp", expr)
        expr = re.sub("M", "(O2 + N2)", expr)
        parsed = parse_expr(expr).simplify()
        gen_dict[gn] = str(parsed)
    # now add generic expressions to all rates and simplify
    srates = []
    for r in rates:
        cr = r
        for gn, gr in gen_dict.items():
            if re.search(f"(^|[^A-Z0-9]){gn}([^A-Z0-9]|$)", cr):
                cr = re.sub(f"(^|[^A-Z0-9]){gn}([^A-Z0-9]|$)", r"\1_REPLACE_ME_\2", cr)
                cr = re.sub("_REPLACE_ME_", gr, cr)
        # other matches
        cr = re.sub("TEMP", "T", cr)
        cr = re.sub("[@]", "**", cr)
        cr = re.sub("[D]", "E", cr)
        cr = re.sub("LOG10", "log10", cr)
        cr = re.sub("EXP", "exp", cr)
        cr = re.sub("M", "(O2 + N2)", cr)
        cr = re.sub(r"J[<](\d+)[>]", r"J\1", cr)
        cr = parse_expr(cr).simplify()
        srates.append(str(cr))
    # regex and templates
    arrhen_regex = r"((?P<A>[-]?\d+([.]\d+(e[+-]\d+)?)?)[*]?|T[*][*](?P<b>[-]?\d+([.]\d+(e[+-]\d+)?)?)[*]?|exp[\(](?P<EaR>[-]?\d+([.]\d+(e[+-]\d+)?)?)[/][T][\)][*]?)+"
    returns = []
    functions = []
    pyrate_file = os.path.join(os.path.dirname(__file__), "mcm_complex_rates.py")
    pyprefix = prefix.replace("-", "_")
    pypref_file = os.path.join(dirname, f"{pyprefix}_rates.py")
    fcn_template = "\ndef KUNKNOWN{}({}):\n    return {}\n"
    repl_exp = "KUNKNOWN{}"
    # add to returns
    for rn, rate in enumerate(srates):
        print(f"Evaluating rate {rate}")
        if re.search(r"[J]\d+", rate):
            prate = re.sub(r"J(\d+)", r"J<\1>", rate)
            rate_data =get_photolysis_parameterization(prate)
        elif re.fullmatch(arrhen_regex, rate):
            match = re.fullmatch(arrhen_regex, rate)
            rate_data = {"rate-constant": {"A": ct.avogadro / 1e6, "b": 0.0, "Ea": 0.0}}
            A = match.group("A")
            b = match.group("b")
            EaR = match.group("EaR")
            rate_data["rate-constant"]["A"] *= float(A) if A is not None else 1
            rate_data["rate-constant"]["b"] = float(b) if b is not None else 0
            rate_data["rate-constant"]["Ea"] =  float(EaR) * ct.gas_constant if EaR is not None else 0
        else:
            cplx_rate = parse_expr(rate) * (ct.avogadro / 1e6)
            cplx_rate = cplx_rate.simplify()
            # get and sort arguments
            variables = [str(sym) for sym in cplx_rate.free_symbols]
            if "T" not in variables:
                variables.append("T")
            variables.sort(reverse=True)
            # setup rate data
            rate_data = {"type": "function-rate"}
            rate_data["function"] = repl_exp.format(rn)
            rate_data["pyfile"] = f"{pyprefix}_rates.py"
            rate_data["ro2file"] = f"{prefix}-ro2-sum.txt"
            # format and add function
            var_str = ", ".join(variables)
            rate_str = str(cplx_rate)
            rate_str = re.sub("log10", "math.log10", rate_str)
            rate_str = re.sub("exp", "math.exp", rate_str)
            fcn = fcn_template.format(rn, var_str, rate_str)
            functions.append(fcn)
        # set units for rate
        rate_data["units"] = {"length":"m", "quantity":"kmol", "activation-energy":"J/kmol"}
        returns.append(rate_data)
    # read in generic rates
    with open(pyrate_file, "r") as f:
        mcm_complex_content = f.read()
    for fcn in functions:
        mcm_complex_content += fcn
    # write out all function rates
    with open(pypref_file, "w") as f:
        f.write(mcm_complex_content)
    return returns
