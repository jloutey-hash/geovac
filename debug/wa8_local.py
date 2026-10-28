"""WA8 derivation: the zero-contact bond-angle terms for four fresh symmetric cases (record
docs/predictions/2026-10-27_water_zero_contact_fresh.md).  Uses debug/wa7_local.py unchanged.
Usage: python debug/wa8_local.py"""
from wa7_local import water

CASES = {'c1': ('both', '2.25', '1.24', ('1.24', '3'), '1.809'),
         'c2': ('one', '2.25', '1.24', ('1.24', '3'), '1.809'),
         'c3': ('both', '1', '1', ('1', '2'), '2.4'),
         'c4': ('one', '1', '1', ('1', '2'), '2.4')}

if __name__ == "__main__":
    for name, args in CASES.items():
        co = water(*args)
        print(name, args, {k: float(x) for k, x in co.items()})
