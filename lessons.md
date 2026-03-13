# Lessons

## 2026-03-10
- When the user narrows the allowed edit surface, update the implementation plan immediately and isolate new functionality into new files or example-local code instead of pushing changes into shared modules.
- When the user points to a project-local environment contract, wire that exact activation sequence into the runnable workflow and verify using that environment instead of assuming the ambient shell is representative.
- When a planned file split is based on memory or setup assumptions, verify the files actually exist in the repo before implementing around them; if the assumption is wrong, collapse the design back into the actual allowed file surface immediately.
- For Upside ligand work, assume vigorous protein rigid-body motion unless proven otherwise; ligand restraints, interaction terms, and accuracy metrics must be formulated in a protein-relative frame, not in the lab frame.
- For Upside ligand smoke tests, do not prove pose retention with backbone-only anchor geometry when the intended physics is sidechain-mediated; make the ligand depend on pocket sidechain placements directly in the test path.
- When upgrading a hybrid model toward a physical force field, do not stop at the interaction form if the atom parameters are still placeholders; finish the parameterization of the coupled subsystem before claiming that path is physically modeled.
