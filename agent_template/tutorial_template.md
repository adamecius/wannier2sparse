<!--
  wannier2sparse tutorial template.
  Copy this file to examples/NN_short_name/README.md and fill the [[ ... ]] slots.
  The <!-- ... --> comments are guidance for the author. They are HTML comments,
  so they stay invisible when the page renders on GitHub. Delete them or leave
  them; readers never see them.

  Read tutorial_style.md once before writing. The short version:
    - lead with the physics, earn every command with a physical reason
    - one controlling idea per tutorial
    - all math in LaTeX, all runnable lines in code blocks
    - no em dashes, no setup or PATH or build paragraphs
    - embed every figure; its alt text states the takeaway
-->

# Tutorial [[N]]: [[the physical question, in plain words]]

<!--
  THE HOOK. Two short paragraphs, no headings.
  Paragraph 1: pose the phenomenon as a small puzzle a curious person would feel.
  Paragraph 2: name what we are after and the single idea this tutorial teaches.
  Do not mention the tool, the build, or the file formats here. Only the physics.
-->

[[Open with the phenomenon. Make the reader curious before they meet a command.]]

[[State plainly what we want to compute and the one lesson that will survive into
later tutorials. End with a sentence of the form "the lesson here is ...".]]

## The physics

<!--
  The minimum theory needed to read the steps, no more. Define any term specific
  to this pipeline (the _hr.dat gauge, ndegen, the spin-label convention). All
  equations in LaTeX. Put the ONE conceptual result the tutorial turns on in its
  own sentence, and say in words why it matters before any code appears.
-->

[[Define the system and its Hamiltonian.]]

$$ [[ H(\mathbf{k}) = \sum_{\mathbf{R}} e^{i\mathbf{k}\cdot\mathbf{R}}\, H(\mathbf{R}) ]] $$

[[Give the closed form or known limit we will compare against, in LaTeX.]]

$$ [[ \text{closed form, e.g. } \rho(E) = \dots ]] $$

[[State the key conceptual result of the tutorial and what it means physically.]]

![[[alt text that states the takeaway of this opening figure]]](fig_[[name]].png)

## Step 1: build [[the system]]

<!--
  Light. One command per size or case. One sentence on what it produces, framed
  as physical objects (a Hamiltonian, its band edges, an orbital set), never as a
  walk through file formats. The reader does not open these files.
-->

```bash
[[python gen_models.py [[system]]   # writes LABEL_hr.dat, LABEL.uc, LABEL.xyz]]
```

[[One sentence: what physical objects this writes, and that the next step reads them.]]

## Step 2: [[expand the operator]]

<!--
  The wannier2sparse call. State which operator and why, then the command.
  Show the one output the reader will refer to. Do not explain CMake, Eigen,
  or "this is the heavy part".
-->

[[One sentence naming the operator and the physical quantity it represents.]]

```bash
[[wannier2sparse [[LABEL]] [[N1 N2 N3]] [[OP ...]] -o out]]
```

[[One sentence: this expands and PBC-wraps the primitive operator, and writes:]]

```
[[out/LABEL.HAM.CSR  and one out/LABEL.<OP>.CSR per operator]]
```

## Step [[k]]: [[a physical statement, not a verb]]

<!--
  Each later step is titled with a physical claim the step demonstrates, e.g.
  "the spectrum is intensive in the supercell size", "the velocity carries the
  Wannier center". Run a command, embed the figure, then read the figure as
  physics.
-->

```bash
[[python w2s_dos.py out/LABEL.HAM.CSR   or   python validate.py]]
```

![[[alt text stating what this figure shows]]](fig_[[name]].png)

[[Read the figure: what the curves do, and the physical reason. If there is a
subtlety the reader would otherwise get wrong, name it plainly here.]]

<!-- Repeat the "Step k" block as needed. Keep the count small. -->

## What to take away

<!-- Three to five bullets. Each is a physics statement, not a software note. -->

- [[the central result, stated as physics]]
- [[the scaling or invariance lesson, e.g. the spectrum is intensive in N]]
- [[the resolution or accuracy lesson, e.g. KPM moments set the broadening]]
- [[the convention or operator lesson, if relevant]]

[[One sentence linking forward: what the next tutorial reuses from this one.]]

## References and links

- wannier2sparse source and documentation: https://github.com/adamecius/wannier2sparse
- Operator and gauge conventions: docs/conventions.md and docs/operators.md.
- Wannier functions: N. Marzari et al., Rev. Mod. Phys. 84, 1419 (2012),
  arXiv:1112.5411. Wannier90: G. Pizzi et al., J. Phys. Condens. Matter 32,
  165902 (2020), arXiv:1907.09788.
- Transport methodology: Z. Fan, J. H. Garcia, A. W. Cummings et al., Linear
  scaling quantum transport methodologies, Phys. Rep. 903, 1 (2021),
  arXiv:1811.07387.
- Installation: see the main README of the repository.

<!--
  Optional, only if the tutorial leans on a specific method or result:
  ## Further reading
  - [[author, title, venue, year]]
  e.g. I. Souza, N. Marzari, D. Vanderbilt, Phys. Rev. B 65, 035109 (2001),
  arXiv:cond-mat/0108084 (disentanglement), grepped from agent_references.md.
-->
</content>
