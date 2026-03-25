// MIT Shocks Simulation: every shock is a surprise (learnt only when it occurs).
// Uses perfect_foresight_with_expectation_errors_solver.
//
// Requires workspace variables:
//   shockssim_mit  (MIT_T_ACTIVE x n_sectors) - shock values for active periods
//
// Requires macro variables (set in model_config.mod):
//   MIT_T_ACTIVE   - number of active periods with shocks
//
// NOTE: periods= must be hardcoded (Dynare limitation). Currently set to 200.
// Update this if changing simul_T_mit in run_dynare_analysis.m.
//
// Shock structure:
//   - Shocks at periods 1 to MIT_T_ACTIVE: each shock is declared with
//     shocks(learnt_in=t), so agents only learn about it when it occurs.
//     After learning, they expect zero future shocks.
@#include "model_config.mod"
@#include "ProdNetRbc_base.mod"

initval;
@#for j in 1:n_sectors
    e_@{j}=0;
    a_@{j}=0;
    k_@{j}=k_ss(@{j});
    c_@{j}=policies_ss(@{j});
    l_@{j}=policies_ss(parn_sectors+@{j});
    pk_@{j}=policies_ss(2*parn_sectors+@{j});
    pm_@{j}=policies_ss(3*parn_sectors+@{j});
    m_@{j}=policies_ss(4*parn_sectors+@{j});
    mout_@{j}=policies_ss(5*parn_sectors+@{j});
    i_@{j}=policies_ss(6*parn_sectors+@{j});
    iout_@{j}=policies_ss(7*parn_sectors+@{j});
    p_@{j}=policies_ss(8*parn_sectors+@{j});
    q_@{j}=policies_ss(9*parn_sectors+@{j});
    y_@{j}=policies_ss(10*parn_sectors+@{j});
@#endfor
cagg = policies_ss(11*parn_sectors+1);
lagg = policies_ss(11*parn_sectors+2);
yagg = policies_ss(11*parn_sectors+3);
iagg = policies_ss(11*parn_sectors+4);
magg = policies_ss(11*parn_sectors+5);
        
end;

steady (solve_algo=3);

endval;
@#for j in 1:n_sectors
    e_@{j}=0;
    a_@{j}=0;
    k_@{j}=k_ss(@{j});
    c_@{j}=policies_ss(@{j});
    l_@{j}=policies_ss(parn_sectors+@{j});
    pk_@{j}=policies_ss(2*parn_sectors+@{j});
    pm_@{j}=policies_ss(3*parn_sectors+@{j});
    m_@{j}=policies_ss(4*parn_sectors+@{j});
    mout_@{j}=policies_ss(5*parn_sectors+@{j});
    i_@{j}=policies_ss(6*parn_sectors+@{j});
    iout_@{j}=policies_ss(7*parn_sectors+@{j});
    p_@{j}=policies_ss(8*parn_sectors+@{j});
    q_@{j}=policies_ss(9*parn_sectors+@{j});
    y_@{j}=policies_ss(10*parn_sectors+@{j});
@#endfor
cagg = policies_ss(11*parn_sectors+1);
lagg = policies_ss(11*parn_sectors+2);
yagg = policies_ss(11*parn_sectors+3);
iagg = policies_ss(11*parn_sectors+4);
magg = policies_ss(11*parn_sectors+5);
        
end;

// periods = MIT_T_ACTIVE = 200 (no burn-in/burn-out for MIT shocks)
// IMPORTANT: Must hardcode this value (Dynare doesn't allow macro vars in periods=).
// If changing simul_T_mit in run_dynare_analysis.m, update this!
perfect_foresight_with_expectation_errors_setup(periods=200);

// Each period's shock is a surprise: learnt_in = the period itself.
// shockssim_mit(t,j) is the shock value for period t, sector j.
@#for t in 1:MIT_T_ACTIVE
shocks(learnt_in=@{t});
  @#for j in 1:n_sectors
    var e_@{j}; periods @{t}; values shockssim_mit(@{t},@{j});
  @#endfor
end;
@#endfor

perfect_foresight_with_expectation_errors_solver;
