from pathlib import Path
p = Path(__file__).parent / "ib_tetra_pf_apparent.py"
s = p.read_text(encoding="utf-8")
s = s.replace("from flint import nmod_poly, fq_default_ctx, fq_default_poly_ctx",
              "from flint import nmod_poly, fq_default_ctx, fq_default_poly_ctx, fmpz_mod_poly_ctx")
s = s.replace("    K = fq_default_ctx(p=P, modulus=f, check_prime=False)",
              "    K = fq_default_ctx(modulus=fmpz_mod_poly_ctx(P)([int(c) for c in f.coeffs()]))")
p.write_text(s, encoding="utf-8")
