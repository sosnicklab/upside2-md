# Adding NMR Structures To ConDiv Training (Future Plan)

This note describes a practical way to add NMR structures to training while preserving stability.

## Goal

Use NMR ensembles to improve local-structure realism without destabilizing current training.

## Key Principle

Treat NMR as an ensemble signal, not a single native target.

## Recommended Strategy

1. Curate NMR entries first
- Keep only high-quality NMR structures (reasonable restraints/completeness).
- Remove entries with severe missing regions or obvious geometry problems.
- Deduplicate against existing proteins (avoid high sequence identity overlap).

2. Represent each NMR target as multiple conformers
- Do not keep only model 1.
- For training, sample one conformer per round/minibatch, or average over a small sampled subset.

3. Add source-aware weighting
- Start with conservative weighting:
  - X-ray weight: `1.0`
  - NMR weight: `0.2` to `0.5`
- Increase NMR weight only if validation improves.

4. Use softer native restraints for NMR
- NMR conformers are usually local minima.
- Use lower restraint strength than X-ray targets to avoid over-constraining.

5. Normalize by protein, not by number of conformers
- If one NMR entry has many models, do not let it dominate.
- Keep total contribution per protein bounded.

## Suggested Data/Metadata Changes

Extend training list metadata to include:
- `source_type` (`xray` or `nmr`)
- `weight` (float)
- `ensemble_id` or list of model paths for NMR proteins

Example conceptual row:
- `3abc source=nmr weight=0.3 models=20`

## Suggested Training Code Changes

1. Dataset loading (`main_initialize` path)
- For NMR proteins, collect available model structures.
- Store model list in the target object/state.

2. Minibatch target selection
- At runtime, sample one model for each NMR protein for that round.
- Keep X-ray behavior unchanged.

3. Gradient aggregation
- Multiply each protein’s contribution by its configured `weight`.
- Keep existing solver and update flow unchanged.

4. Restraint strength routing
- Use source-specific restraint constants (e.g., separate NMR/X-ray defaults).

## Minimal First Experiment

1. Add a small NMR subset (10 to 20 proteins).
2. Use one sampled model per NMR protein per minibatch.
3. Set `nmr_weight=0.3`.
4. Set NMR restraint strength lower than X-ray.
5. Run 1 to 2 epochs and compare against X-ray-only baseline.

## Success Criteria

- No increase in crash/explosion rate.
- Gradient norms remain in the same rough scale as baseline.
- Validation metrics improve or stay neutral while local structure quality improves.

## Risks

- Overweighting NMR can bias toward local minima and hurt global fold behavior.
- Counting all NMR models as independent proteins can severely skew training.
- Poorly curated NMR sets can inject noise and destabilize updates.
