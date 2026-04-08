// First-Order (Linear) Approximation Solution
//
// This file avoids the Dynare command/file-name collision from stoch_simul.mod
// when running under Octave.

@#include "model_config.mod"
@#include "base_model_include.mod"

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
c_util = policies_ss(11*parn_sectors+1);
l_util = policies_ss(11*parn_sectors+2);
c_agg = policies_ss(11*parn_sectors+3);
l_agg = policies_ss(11*parn_sectors+4);
gdp_agg = policies_ss(11*parn_sectors+5);
i_agg = policies_ss(11*parn_sectors+6);
k_agg = policies_ss(11*parn_sectors+7);
utility_intratemp = policies_ss(11*parn_sectors+8);

end;

shocks;
@#for j in 1:n_sectors
   var e_@{j} = parSigma_A(@{j},@{j});
@#endfor
@#for j in 1:n_sectors-1
    @#for i in j+1:n_sectors
        var e_@{j}, e_@{i} = parSigma_A(@{i},@{j});
    @#endfor
@#endfor
end;

// IMPORTANT:
// - keep the full endogenous list (no trailing var list) so oo_.var and
//   oo_.autocorr include the full state/policy system needed by TheoStats.
// - omit the periods option so Dynare reports theoretical moments in oo_.var.
stoch_simul(order=1, irf=0, ar=1, nocorr, nograph, nofunctions);
