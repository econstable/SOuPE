# Scientific sources and scope

The main scientific references are:

- Amelin et al., [Phys. Rev. B 102, 134428 (2020)](https://doi.org/10.1103/PhysRevB.102.134428).
- Alexanian et al., [Phys. Rev. B 107, 224404 (2023)](https://doi.org/10.1103/PhysRevB.107.224404), including its supplementary material.
- Supporting detail: [Alexanian's thesis](https://hal.science/tel-03668858v1).
Rafik Ballou's [Stevens operator table](Stevens_Operators-Ballou.pdf).

| Component | Scientific basis | Repository component |
| --- | --- | --- |
| Stevens polynomials | R. Ballou's Stevens operator table; selected definitions in Amelin et al., Appendix A | `_stevens.py`; source-table exceptions documented in CONVENTIONS.md. |
| Coefficient projection | Amelin et al., Appendix B, Eq. (B5); Alexanian et al., supplementary Eq. (16) | `conversions.py`, B=lambda*theta*L. |
| Published numerical factors | Amelin et al., Tables II-III; Alexanian et al., supplementary Tables II-III | Explicit factors in the configs; CONVERSIONS.md distinguishes published entries from the general table. Additional orders support lower-symmetry crystal fields; normalisation factors follow the McPhase manual, using the stated real combinations.|
| Full 27-component conversion | General tensor relation in supplementary Eq. (16), supported by thesis Eq. (III.18) | Remaining factors follow Racah tensor and polynomial normalisations; a documented generalisation for lower symmetry. |
| Euler orientation | Local/global calculations in Amelin et al., Eqs. (5-7); axes in Alexanian et al., supplementary Table VII | `rotations.py`; config angles are derived from the axes, with the derivation in the example guide. |
| Stevens input example | Amelin et al., Table I and section IV | `tb2ti2o7_Amelin2020.json`, g=1.4; static and sweep scripts. |
| Wybourne input example | Alexanian et al., supplementary Eq. (1), Tables I-III and VII | `tb2ti2o7_Alexanian2023.json`; unit scaling and explicit theta/lambda conversion. |
| Illustrative Alexanian2023 field setting | Supporting thesis, section V.3.2.1 | g=1.35 is sourced separately; the default static field is zero. |
| Zero-field state reporting | Basis selection within a degenerate eigenspace using projected local Jz | Example-only helper; no physical perturbation or change of energy/subspace. |
| Energy units | E=hc*wavenumber, E=k_B*T, exact SI constants | `units.py`; utility functionality rather than a physical model. |
| LaTeX | O_k^q notation in the publications and operator table | `latex.py`; expressions are rendered from implementation syntax. |
| QuTiP copies, storage and input validation | Software requirements of an operator utility | Public API implementation choices. |

The examples construct independent-ion CEF and Zeeman Hamiltonians. Their
config reader, reporting basis selection and TSV output are pedagogical
example code, not package-level modelling or configuration APIs.