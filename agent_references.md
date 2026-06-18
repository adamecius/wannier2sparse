# agent_references.md — reference manifest (committed; the PDFs are not)

This is the **committed manifest** of the reference material an agent may consult
to justify a convention, a Fourier-transform sign, a spin/orbital unit, or a
disentanglement detail. Per `AGENTS.md` §5, the reference **files** live in
`agent_references/` and are **gitignored, local only** — only this manifest is in
the repo. Never `git add` a PDF; never paste its content into the repo. To use
one: grep this file by topic, open the local file (or its arXiv page), and cite
it in the PR body by title.

The two primary references for this project are the **Quantum ESPRESSO** and
**Wannier90** user guides; the method papers below give the physics behind the
operators `wannier2sparse` builds. arXiv ids are given where one exists so an
agent can fetch the source without the local PDF. The maintainer will review and
prune this list.

Manifest entry format:

```
## <Author, Venue Year> — <Title>
- file:    agent_references/<filename>   (local only, not committed)
- arxiv:   <id or "none (user guide)">
- topics:  <semicolon-separated topics>
- use-for: <when an agent should reach for this>
```

---

## User guides (primary, no arXiv)

## Wannier90 — User Guide (v3.x)
- file:    agent_references/wannier90_user_guide.pdf   (local only, not committed)
- arxiv:   none (user guide; wannier90.org/support and the docs/ tree of the release)
- topics:  _hr.dat format; ndegen; use_ws_distance / _wsvec.dat; _u.mat / _u_dis.mat; .amn / .spn / .eig; projections block; real-harmonic (l, mr) order
- use-for: the authoritative file-format and projection-order reference; pair with docs/conventions.md, which cites the exact source lines (W90 3.1.0)

## Quantum ESPRESSO — User Guide + pw2wannier90 / pw.x input docs
- file:    agent_references/qe_user_guide.pdf   (local only, not committed)
- arxiv:   none (user guide; INPUT_PW.html, pw2wannier90 docs, data-file-schema.xml)
- topics:  SCF/NSCF workflow; noncollinear/SOC (lspinorb, noncolin); write_spn / write_amn / write_unk; data-file-schema.xml provenance (lattice, atoms, symmetry, k-mesh, XC, pseudos, ecutwfc)
- use-for: how the DFT inputs that feed Wannier90 are produced; the provenance fields parsed for the bundle manifest

---

## Wannier functions and the gauge (method papers)

## Marzari & Vanderbilt, Phys. Rev. B 56, 12847 (1997) — Maximally localized generalized Wannier functions for composite energy bands
- file:    agent_references/marzari_vanderbilt_1997.pdf   (local only, not committed)
- arxiv:   cond-mat/9707145
- topics:  maximally localized Wannier functions; the gauge U(k); real-space O(R)
- use-for: the definition of the Wannier gauge U(k) that _u.mat stores and that the exact-spin / orbital-L routes carry

## Souza, Marzari & Vanderbilt, Phys. Rev. B 65, 035109 (2001) — Maximally localized Wannier functions for entangled energy bands
- file:    agent_references/souza_marzari_vanderbilt_2001.pdf   (local only, not committed)
- arxiv:   cond-mat/0108084
- topics:  disentanglement; U_dis(k) (num_bands x num_wann); the optimal subspace
- use-for: the bookkeeping behind V(k) = U_dis(k)·U(k) and _u_dis.mat (present only for entangled runs); docs/conventions.md §2

## Marzari, Mostofi, Yates, Souza & Vanderbilt, Rev. Mod. Phys. 84, 1419 (2012) — Maximally localized Wannier functions: theory and applications
- file:    agent_references/marzari_rmp_2012.pdf   (local only, not committed)
- arxiv:   1112.5411
- topics:  Wannier interpolation; H(k) = Σ_R e^{ik·R} H(R); operators in the Wannier basis
- use-for: the umbrella review; the interpolation identity behind docs/operators.md §1

## Mostofi, Yates, Pizzi, Lee, Souza, Vanderbilt & Marzari, Comput. Phys. Commun. 185, 2309 (2014) — An updated version of wannier90
- file:    agent_references/mostofi_w90_2014.pdf   (local only, not committed)
- arxiv:   1309.1827
- topics:  wannier90 v2; _hr.dat; interpolation machinery
- use-for: the code that produced our _hr.dat fixtures; format provenance

## Pizzi et al., J. Phys. Condens. Matter 32, 165902 (2020) — Wannier90 as a community code: new features and applications
- file:    agent_references/pizzi_w90_2020.pdf   (local only, not committed)
- arxiv:   1907.09788
- topics:  wannier90 v3; pw2wannier90 interface; spin (.spn), orbital projections
- use-for: the current Wannier90 reference (v3.x, our verified version line); cite for any W90 feature

---

## Quantum ESPRESSO (method papers)

## Giannozzi et al., J. Phys. Condens. Matter 21, 395502 (2009) — QUANTUM ESPRESSO: a modular and open-source software project
- file:    agent_references/qe_2009.pdf   (local only, not committed)
- arxiv:   0906.2569
- topics:  plane-wave DFT; pw.x; the SCF engine upstream of pw2wannier90
- use-for: the canonical QE citation for the DFT half of the pipeline

## Giannozzi et al., J. Phys. Condens. Matter 29, 465901 (2017) — Advanced capabilities for materials modelling with Quantum ESPRESSO
- file:    agent_references/qe_2017.pdf   (local only, not committed)
- arxiv:   1709.10010
- topics:  QE 6.x/7.x capabilities; noncollinear/SOC; data-file-schema.xml era
- use-for: the modern QE citation (our verified version line is QE 7.2)

---

## Spin and orbital operators in the Wannier basis

## Lopez, Vanderbilt, Thonhauser & Souza, Phys. Rev. B 85, 014435 (2012) — Wannier-based calculation of the orbital magnetization in crystals
- file:    agent_references/lopez_orbital_2012.pdf   (local only, not committed)
- arxiv:   1112.1938   (verify id when fetching)
- topics:  orbital angular momentum / magnetization via Wannier interpolation; local L; projector route
- use-for: the physics behind the --orbital-L operator (units ħ; the L_local + C(k)=A†V route); docs/conventions.md §5

## Tsirkin, npj Comput. Mater. 7, 33 (2021) — High performance Wannier interpolation of Berry curvature and related quantities with WannierBerri code
- file:    agent_references/tsirkin_wannierberri_2021.pdf   (local only, not committed)
- arxiv:   2008.07992
- topics:  spin / Berry-curvature Wannier interpolation; V†O_B V gauge transform; use_ws_distance conventions
- use-for: the independent cross-check codebase used for the committed golden in docs/conventions.md §6–§7; WS-convention trap

---

## Transport methodology (downstream consumer)

## Weisse, Wellein, Alvermann & Fehske, Rev. Mod. Phys. 78, 275 (2006) — The kernel polynomial method
- file:    agent_references/weisse_kpm_2006.pdf   (local only, not committed)
- arxiv:   cond-mat/0504627
- topics:  KPM; Chebyshev moments; Jackson kernel; spectral functions from sparse H
- use-for: why the CSR output exists and how the examples reconstruct the DOS; cite in example tutorials that use kernel details

## Fan, Garcia, Cummings et al., Phys. Rep. 903, 1 (2021) — Linear scaling quantum transport methodologies
- file:    agent_references/fan_linqt_2021.pdf   (local only, not committed)
- arxiv:   1811.07387
- topics:  LinQT; Kubo-Greenwood / Bastin; real-space sparse operators; the package wannier2sparse feeds
- use-for: the methodology footer for tutorials and the downstream context for the operators built here
</content>
