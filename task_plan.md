Project Goal
- Write the `ConDiv_symlay` training-potential and constraint design principles into a standalone Markdown document.

Architecture & Key Decisions
- Put the design note under `ConDiv_symlay/` so it stays next to the workflow code and README.
- Keep the note focused on the actual implemented design:
  - direct implicit membrane training objective
  - hard symmetric-layer projection
  - soft dry-MARTINI distance-resolved teacher
  - fixed DOPC thickness
  - preservation of the protein-vs-bilayer frame already encoded in training coordinates
- Add one short README link so the note is discoverable without duplicating the full content there.

Execution Phases
- [x] Choose the doc location and scope.
- [x] Add the standalone Markdown design note.
- [x] Link the note from `ConDiv_symlay/README.md`.
- [x] Update tracking files and verify the new documentation content.

Known Errors / Blockers
- None.

Review
- Added [TRAINING_POTENTIAL_DESIGN.md](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/TRAINING_POTENTIAL_DESIGN.md) under `ConDiv_symlay/`.
- The note captures the implemented design in one place:
  - implicit membrane runtime potential
  - total update rule
  - hard symmetric DOPC layer projection
  - soft distance-resolved dry-MARTINI teacher
  - fixed DOPC thickness
  - preservation of the protein-bilayer frame already encoded in training coordinates
- Added a pointer to the new note near the top of [README.md](/Users/yinhan/Documents/upside2-md/ConDiv_symlay/README.md).
- Verification was a direct read-back of the new Markdown file and README link. No code execution changes were involved in this task.
