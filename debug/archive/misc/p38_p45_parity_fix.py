"""Make the CG-sum integrality (unit-step) constraint explicit in
Paper 45 Appendix A and Paper 38 L2(b). Values were correct; the
displayed summation condition was underspecified (caught by the
Phase-B scalar prototype when sigma(1/2) > 1 appeared under the loose
triangle reading)."""
import io

fails = []


def rep(path, old, new, tag):
    s = io.open(path, encoding='utf-8').read()
    if s.count(old) == 1:
        io.open(path, 'w', encoding='utf-8').write(s.replace(old, new))
    else:
        fails.append((tag, s.count(old)))


P45 = r'papers/group1_operator_algebras/paper_45_lorentzian_propinquity.tex'
P38 = r'papers/group1_operator_algebras/paper_38_su2_propinquity_convergence.tex'

rep(P45, r'''   c_{j'} \;=\; \sum_{\substack{j_1, j_2 \le j_{\max} \\
                |j_1 - j_2| \le j' \le j_1 + j_2}}
                \sqrt{(2j_1+1)(2j_2+1)},''',
    r'''   c_{j'} \;=\; \sum_{\substack{j_1, j_2 \le j_{\max} \\
                |j_1 - j_2| \le j' \le j_1 + j_2 \\
                j_1 + j_2 + j' \in \mathbb{Z}}}
                \sqrt{(2j_1+1)(2j_2+1)},''', 'p45-csum')

rep(P45, r'''Clebsch--Gordan rule contributes each admissible $j'$ once per pair;''',
    r'''Clebsch--Gordan rule contributes each admissible $j'$ once per pair
--- $j'$ running in \emph{unit} steps from $|j_1 - j_2|$ to
$j_1 + j_2$, so only pairs with $j_1 + j_2 + j' \in \mathbb{Z}$
contribute;''', 'p45-proof')

rep(P38, '\n'.join([
    r"$c_{j'} = \sum_{j_1, j_2 \le j_{\max},\, \triangle(j_1, j_2, j')}",
    r"\sqrt{(2j_1+1)(2j_2+1)}$ (Clebsch--Gordan triangle rule):\ ",
]).rstrip() + '\n' + r"$\sigma_{\nmax}(0) = 1$",
    '\n'.join([
    r"$c_{j'} = \sum_{j_1, j_2 \le j_{\max},\, \triangle(j_1, j_2, j')}",
    r"\sqrt{(2j_1+1)(2j_2+1)}$, where $\triangle$ is the Clebsch--Gordan",
    r"condition \emph{with unit steps} ($|j_1 - j_2| \le j' \le j_1 + j_2$",
    r"and $j_1 + j_2 + j' \in \mathbb{Z}$):\ ",
]).rstrip() + '\n' + r"$\sigma_{\nmax}(0) = 1$", 'p38-csum')

rep(P38, '\n'.join([
    r"symbol follows by expanding $|h|^{2}$ with the SU(2) Clebsch--Gordan",
    r"rule $\chi_{j_1}\chi_{j_2} = \sum_{J = |j_1 - j_2|}^{j_1+j_2}\chi_{J}$",
]),
    '\n'.join([
    r"symbol follows by expanding $|h|^{2}$ with the SU(2) Clebsch--Gordan",
    r"rule $\chi_{j_1}\chi_{j_2} = \sum_{J = |j_1 - j_2|}^{j_1+j_2}\chi_{J}$",
    r"(unit steps in $J$, so $j_1 + j_2 + J \in \mathbb{Z}$)",
]), 'p38-proof')

print('FAILS:', fails if fails else 'none')
