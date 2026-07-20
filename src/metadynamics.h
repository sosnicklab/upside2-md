#ifndef METADYNAMICS_H
#define METADYNAMICS_H

#include "deriv_engine.h"
#include <cstdint>

// Well-tempered metadynamics bias on a single collective variable (a CoordNode scalar).
// The bias node reads the CV, adds a sum-of-Gaussians bias energy, and pushes dV/ds back onto the
// CV's sensitivity (the engine propagates it to the atoms). Deposition of new Gaussians is driven
// once per MD round by maybe_deposit(), keyed off the true round counter (compute() runs several
// times per round, so deposition must not be tied to compute_value).
namespace metadynamics {
    // Deposit a Gaussian in every metadynamics node of the engine that is due at this round.
    void maybe_deposit(DerivEngine& engine, uint64_t round_num, float temperature);
}

#endif
