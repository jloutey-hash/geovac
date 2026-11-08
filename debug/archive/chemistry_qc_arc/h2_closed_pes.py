"""Thin wrapper: the closed-form H2 PES now lives in geovac.qfd_assemble
(h2_closed_form_E, h2_pes_certify); this driver just runs the certification
and writes debug/data/h2_closed_pes.json.  See CHANGELOG v4.106.0."""
import json

from geovac.qfd_assemble import h2_closed_form_E, h2_pes_certify

if __name__ == "__main__":
    E, out = h2_pes_certify()
    for k, v in out.items():
        print(f"{k} = {v}")
    json.dump(out, open("debug/data/h2_closed_pes.json", "w"), indent=1)
