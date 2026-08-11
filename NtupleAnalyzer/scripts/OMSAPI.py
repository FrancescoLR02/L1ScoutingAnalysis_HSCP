import sys
import json
from omsapi import OMSAPI

my_app_id = 'cms-l1scout-analysis'
my_app_secret = 'XXX'

omsapi = OMSAPI("https://cmsoms.cern.ch/agg/api", "v1", cert_verify=False)
omsapi.auth_oidc(my_app_id, my_app_secret)

iruns = 394633
fruns = 398903

run_to_scheme = {}
run_to_lumi = {}
run_to_fill = {}
total_lumi = 0.0

debug_printed = False

for run in range(iruns, fruns + 1):
    # --- Step 1: get the run, including fill_number and delivered_lumi ---
    q = omsapi.query("runs")
    q.filter("run_number", run)
    q.attrs(["run_number", "fill_number", "delivered_lumi"])
    q.paginate(page=1, per_page=1)
    resp = q.data()

    data = resp.json().get("data", [])
    if not data:
        print(f"Run {run}: NOT FOUND in OMS")
        run_to_scheme[run] = None
        run_to_lumi[run] = 0.0
        continue

    attrs = data[0]["attributes"]
    fill_number = attrs.get("fill_number")
    lumi = attrs.get("delivered_lumi", 0.0) or 0.0

    if fill_number is None:
        print(f"Run {run}: no fill_number found")
        run_to_scheme[run] = None
        run_to_lumi[run] = lumi
        total_lumi += lumi
        continue

    # --- Step 2: get the fill, and look at its filling scheme ---
    qf = omsapi.query("fills")
    qf.filter("fill_number", fill_number)
    qf.paginate(page=1, per_page=1)
    fresp = qf.data()
    fdata = fresp.json().get("data", [])

    if not fdata:
        print(f"Run {run} (fill {fill_number}): fill NOT FOUND in OMS")
        run_to_scheme[run] = None
        run_to_lumi[run] = lumi
        total_lumi += lumi
        continue

    fill_attrs = fdata[0]["attributes"]

    # One-time debug dump so we can find the right key name
    if not debug_printed:
        print("=== Available attributes on 'fills' endpoint (first match) ===")
        print(json.dumps(fill_attrs, indent=2))
        print("================================================================\n")
        debug_printed = True

    # NOTE: adjust this key once you've identified the right one from the debug dump above
    scheme = fill_attrs.get("fill_scheme", fill_attrs.get("injection_scheme", "UNKNOWN"))

    run_to_scheme[run] = scheme
    run_to_lumi[run] = lumi
    run_to_fill[run] = fill_number
    total_lumi += lumi

    print(f"Run {run} (fill {fill_number}): scheme={scheme}, delivered_lumi={lumi}")

print(f"\nTotal delivered luminosity: {total_lumi}")

# --- Write C++ header file, grouped by fill number ---
with open("run_bxscheme_map.h", "w") as f:
    f.write("#pragma once\n\n")
    f.write("#include <map>\n#include <string>\n\n")
    f.write("std::map<int, std::string> runToScheme = {\n")

    prev_fill = None
    for run in sorted(run_to_scheme.keys(), reverse=True):
        scheme = run_to_scheme[run]
        if scheme is None:
            continue
        fill = run_to_fill.get(run, "UNKNOWN")

        if fill != prev_fill:
            if prev_fill is not None:
                f.write("\n")  # blank line between fill groups
            f.write(f"   //{fill}\n")
            prev_fill = fill

        f.write(f'   {{{run}, "{scheme}"}},\n')

    f.write("};\n")

print("\nSaved to run_bxscheme_map.h")

# --- Also dump JSON with scheme + lumi for reference ---
with open("run_bxscheme_map.json", "w") as f:
    json.dump(
        {
            "run_to_scheme": run_to_scheme,
            "run_to_lumi": run_to_lumi,
            "run_to_fill": run_to_fill,
            "total_lumi": total_lumi,
        },
        f,
        indent=2,
    )

print("Saved to run_bxscheme_map.json")