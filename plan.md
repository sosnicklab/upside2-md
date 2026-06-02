**Observations & Primary Objective:**
In the workflows `example/16.MARTINI/run_sim_1afo.sh` and `example/16.MARTINI/run_sim_1rkl.sh`, the coarse-grained lipid (CGL) orientation appears largely correct in the bulk. However, the critical issue is that **anomalous, unphysical CGL orientations are visible directly at the protein-bilayer interfaces** (as well as at the box boundaries).

**Action Required:**
The explicit goal of this task is to **fix the distorted CGL orientations around the protein**. Please review the simulation results (`example/16.MARTINI/outputs/martini_1afo_hybrid/1afo.stage_7.0.vtf` and `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf`) to determine whether these interfacial artifacts stem from the force field parameterization or the visualization software, and correct them accordingly.

**Potential Fitting Methodology:**
To systematically calculate and fit the potentials, evaluate the pairwise interactions between two entities (molecules or coarse-grained particles) in a simulation box using the dryMARTINI force field. The full spline table potential is constructed by exhaustively sampling the configurational space through the following sequence:

1. **Spatial and Orientational Sampling:** Place the first entity at the origin and systematically translate the second entity across a range of radial distances within a surrounding sphere. At each coordinate, rotate both entities through all possible relative orientations (if one entity is an isotropic particle, rotate only the molecule).
2. **Energy Minimization (Steric Relaxation):** Because exhaustive spatial sampling artificially forces molecules into close proximity, high-energy steric clashes are inevitable. A rigorous energy minimization step must be applied at each configuration to resolve these overlaps before recording the interaction. This prevents unphysical repulsive spikes from corrupting the spline table.
3. **Conformational Constraint (Plane/Axis Restriction):** During sampling, the movement of the molecules must be restricted along the axis or plane of their molecular vectors. Without this constraint, coarse-grained lipids are prone to unnatural lateral bending and distortion when subjected to artificial orientations. This restriction ensures that the calculated potentials reflect stable, physically realistic lipid conformations.
4. **Orientation Vector Definition (Relative vs. Absolute Reference Frames):** Unlike protein side chains, which possess a fixed absolute orientation anchored to a rigid backbone root, a coarse-grained lipid (CGL) exists in a rotationally invariant fluid environment. The absolute spatial orientation of a single CGL is physically meaningless; the interaction energy is dictated entirely by its *relative* orientation to its neighbors. Therefore, while absolute spatial coordinates are discarded, the orientational vectors ($\mathbf{n}$) must be strictly retained to preserve the anisotropic, cylindrical geometry of the molecules (preventing them from being mathematically treated as isotropic spheres):
* **CGL:** Defined by its macroscopic principal axis (e.g., primary head group to terminal tail bead).
* **Side Chain (SC):** Defined by its standard root-to-tip axis.



**Suggested Mathematical Formulations (Optional):**
*Note: The following mathematical potential forms and Singular Value Decomposition (SVD) strategies are provided as suggestions only. You are not strictly required to use these exact equations, provided you adhere to the hard constraints listed below.*

First, the purely distance-dependent baseline, $V_{\text{radial}}(r_{12})$, is extracted by averaging the interaction energy over all sampled orientations at each specific distance $r_{12}$. Once this isotropic baseline is subtracted, SVD is applied to the anisotropic residual.

* **Suggested Form 1: CGL-CGL & SC-CGL (Two-Vector Interactions)**
To ensure the potential relies strictly on relative alignment, the interaction incorporates three distinct relative angular dependencies:

$$V = \kappa \left[ V_{\text{radial}}(r_{12}) + \text{Ang}_1(-\mathbf{n}_1 \cdot \mathbf{n}_{12}) \cdot \text{Ang}_2(\mathbf{n}_2 \cdot \mathbf{n}_{12}) \cdot \text{Ang}_3(\mathbf{n}_1 \cdot \mathbf{n}_2) \cdot V_{\text{angular}}(r_{12}) \right]$$


* **Suggested Form 2: CGL-Particle & SC-Particle (Vector-Point Interactions)**
The interaction depends solely on the particle's spatial position relative to the molecule's axis, defined by a single relative angle $\theta$:

$$V = \kappa \left[ V_{\text{radial}}(r_{12}) + \text{Ang}(\cos \theta) \cdot V_{\text{angular}}(r_{12}) \right]$$



**Testing & Iterative Validation Protocol:**
To ensure stability and accurately resolve the interfacial orientations, the force field must be developed through a strict, iterative validation loop:

1. **Phase 1 (Bilayer Baseline):** Construct the CGL-CGL force field first. Test this potential exclusively on a bilayer-only model to verify the formation of a structurally stable bilayer with consistent, physically realistic CGL orientations.
2. **Phase 2 (Comprehensive Application):** Once the CGL-CGL baseline is validated, apply the validated fitting method to generate all remaining potential forms (SC-CGL, SC-particle, CGL-particle).
3. **Phase 3 (Hybrid System Validation):** Test the complete, integrated force field on the target hybrid protein-membrane systems (`1afo` and `1rkl`).
4. **Phase 4 (Iterative Refinement):** Evaluate the hybrid simulations against three strict success metrics: a stable bilayer, **consistent CGL orientations (especially eliminating the anomalous alignments at the protein interface)**, and stable protein structures. If any of these metrics fail, the underlying generation methodology must be adjusted, and the entire pipeline—starting over from Phase 1—must be repeated until the systems are fully stabilized.

**System Constraints:**

* **Hard Constraints (Strictly Mandatory):**
* **Physical Accuracy:** All interactions must remain strictly physical. Do not apply any parameter tweaking or twisting.
* **Force Field Sourcing:** All calculations must be derived directly from the `.itp` force field files.
* **No Artificial Potentials:** Do not apply any additional orientation potentials to the CGLs.
* **No Force or Energy Capping:** Do not apply any force or energy capping. These artificial limits are unphysical and strictly prohibited.


* **Soft Constraints (Flexible):**
* **Universal Methodology:** Aim to keep the computational method universally consistent across all interaction pairs (CGL-CGL, SC-CGL, CGL-particle, SC-particle). However, if enforcing a universal method proves technically prohibitive, you may break this rule and adapt the methodology per interaction type, provided all hard constraints are still met.



**Documentation Requirement:**

* **Manuscript Preparation:** Continuously update `example/16.MARTINI/cg_lipid_potentials.tex` concurrently with the debugging and validation process. This file must accurately document the final, tested procedures so that it is fully prepared to serve as the Methods section of a publication.