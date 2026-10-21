#include "martini.h"
#include "box.h"

void martini_unregister_engine(DerivEngine& engine) {
    martini_masses::unregister_masses_for_engine(&engine);
    martini_fix_rigid::unregister_fix_rigid_for_engine(engine);
    martini_stage_params::unregister_stage_params_for_engine(&engine);
    martini_hybrid::unregister_hybrid_for_engine(engine);
    martini_brownian::unregister_brownian_for_engine(&engine);
    simulation_box::npt::unregister_barostat_for_engine(engine);
}
