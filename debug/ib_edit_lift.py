from pathlib import Path
p = Path(__file__).parent / "ib_tetra_lift.py"
s = p.read_text(encoding="utf-8")
old_start = s.index("    # Q-form: x^order")
old_end = s.index("    out = {")
block = s[old_start:old_end]
s = s[:old_start] + "    Q, lo = d_to_Q(ip)\n" + s[old_end:]
func = '''

def d_to_Q(ip):
    """integer d/dx form on S = Y/t^3 -> Q-form on Y in u = theta_t; returns (list of sympy exprs, x shift)."""
    order = len(ip) - 1
''' + block.replace("    # Q-form", "    # Q-form", 1) + "    return Q, lo\n"
s = s.replace("\n\ndef main():", func + "\n\ndef main():", 1)
p.write_text(s, encoding="utf-8")
