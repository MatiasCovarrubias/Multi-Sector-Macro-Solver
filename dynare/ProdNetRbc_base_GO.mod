load ModStruct_temp.mat
var 
    @#for j in 1:n_sectors
        k_@{j}
    @#endfor
    @#for j in 1:n_sectors
        a_@{j}
    @#endfor
    @#for j in 1:n_sectors
        c_@{j}
    @#endfor
    @#for j in 1:n_sectors
        l_@{j}
    @#endfor
    @#for j in 1:n_sectors
        pk_@{j}
    @#endfor
    @#for j in 1:n_sectors
        pm_@{j}
    @#endfor
    @#for j in 1:n_sectors
        m_@{j}
    @#endfor
    @#for j in 1:n_sectors
        mout_@{j}
    @#endfor
    @#for j in 1:n_sectors
        i_@{j}
    @#endfor
    @#for j in 1:n_sectors
        iout_@{j}
    @#endfor
    @#for j in 1:n_sectors
        p_@{j}
    @#endfor
    @#for j in 1:n_sectors
        q_@{j}
    @#endfor
    @#for j in 1:n_sectors
        y_@{j}
    @#endfor
    c_util l_util c_agg l_agg gdp_agg i_agg k_agg utility_intratemp
    ;

predetermined_variables 
    @#for j in 1:n_sectors
        k_@{j}
    @#endfor
    @#for j in 1:n_sectors
        a_@{j}
    @#endfor
    ;

varexo 
    @#for j in 1:n_sectors
        e_@{j}
    @#endfor
    ;

parameters eps_l eps_c beta phi theta sigma_c sigma_m sigma_q sigma_y sigma_I sigma_l
    @#for j in 1:n_sectors
        xi_@{j} alpha_@{j} mu_@{j} rho_@{j} delta_@{j} pss_@{j} pkss_@{j} pmss_@{j} Css_@{j}
        @#for i in 1:n_sectors
            Gamma_M_@{i}_@{j} Gamma_I_@{i}_@{j}
        @#endfor
    @#endfor
;


eps_l = pareps_l;
eps_c = pareps_c;
beta=parbeta;
phi=parphi;
theta = partheta;
sigma_c = parsigma_c;
sigma_m = parsigma_m;
sigma_q = parsigma_q;
sigma_y = parsigma_y;
sigma_I = parsigma_I;
sigma_l = parsigma_l;
@#for j in 1:n_sectors
    pss_@{j}=p_ss_log(@{j});
    pkss_@{j}=pk_ss_log(@{j});
    pmss_@{j}=policies_ss(3*parn_sectors+@{j});
    Css_@{j}=exp(policies_ss(@{j}));
    xi_@{j}=parxi(@{j});
    alpha_@{j}=paralpha(@{j});
    mu_@{j}=parmu(@{j});
    rho_@{j}=parrho(@{j});
    delta_@{j}=pardelta(@{j});
    @#for i in 1:n_sectors
        Gamma_M_@{i}_@{j}=parGamma_M(@{i},@{j});
        Gamma_I_@{i}_@{j}=parGamma_I(@{i},@{j});
    @#endfor
@#endfor

model;
# MU = (exp(c_util) - theta*(1/(1+eps_l^(-1))) * exp(l_util)^(1+eps_l^(-1)))^(-eps_c^(-1));

@#for j in 1:n_sectors
    #   Pktil_@{j} = (0
        @#for i in 1:n_sectors
            +Gamma_I_@{i}_@{j}*exp(p_@{i})^(1-sigma_I)
        @#endfor
        )^(1/(1-sigma_I));  
 
    exp(k_@{j}(+1))=(1-delta_@{j})*exp(k_@{j})+exp(i_@{j})-phi/2 * exp(k_@{j}) * ( exp(i_@{j})/exp(k_@{j}) - delta_@{j} )^2;
    a_@{j}(+1)=rho_@{j}*a_@{j}+e_@{j};
    
    exp(p_@{j}) = (exp(c_util) * xi_@{j} / exp(c_@{j}))^(1/sigma_c);

    theta * exp(l_util)^(eps_l^(-1))*(exp(l_@{j})/exp(l_util))^(1/sigma_l) = (exp(p_@{j})) * (exp(a_@{j}))^((sigma_q-1)/sigma_q) * ( (mu_@{j}) * exp(q_@{j})/exp(y_@{j}) )^(1/sigma_q) * ( (1-alpha_@{j}) * exp(y_@{j})/exp(l_@{j}) )^(1/sigma_y);
    
    exp(pk_@{j}) = beta * (((exp(c_util(+1)) - theta*(1/(1+eps_l^(-1))) * exp(l_util(+1))^(1+eps_l^(-1)))^(-eps_c^(-1))) / MU) * (exp(p_@{j}(+1)) * exp(a_@{j}(+1))^((sigma_q-1)/sigma_q) * ((mu_@{j}) * exp(q_@{j}(+1))/exp(y_@{j}(+1)))^(1/sigma_q) * (alpha_@{j} * exp(y_@{j}(+1))/exp(k_@{j}(+1)))^(1/sigma_y) + exp(pk_@{j}(+1)) * ((1-delta_@{j}) + phi/2 * ((exp(i_@{j}(+1))/exp(k_@{j}(+1)))^2 - delta_@{j}^2)));

    exp(pm_@{j}) = (0
        @#for i in 1:n_sectors
            +Gamma_M_@{i}_@{j}*exp(p_@{i})^(1-sigma_m)
        @#endfor
        )^(1/(1-sigma_m));

    exp(m_@{j}) = (1-mu_@{j}) * (exp(pm_@{j})/((exp(a_@{j}))^((sigma_q-1)/sigma_q)*exp(p_@{j})))^(-sigma_q) * exp(q_@{j});

    exp(mout_@{j}) = (0
            @#for i in 1:n_sectors
                +Gamma_M_@{j}_@{i}*(exp(p_@{j})/exp(pm_@{i}))^(-sigma_m) * exp(m_@{i})
            @#endfor
            ) ;  
    

    exp(pk_@{j}) = Pktil_@{j} * ( 1 - phi*( exp(i_@{j})/exp(k_@{j}) - delta_@{j} ) )^(-1);
    
    exp(iout_@{j}) = (0
            @#for i in 1:n_sectors
                +Gamma_I_@{j}_@{i}*(exp(p_@{j})/exp(pk_@{i}))^(-sigma_I) * exp(i_@{i}) * ( 1 - phi*( exp(i_@{i})/exp(k_@{i}) - delta_@{i} ) )^(sigma_I)
            @#endfor
            ) ;  
    exp(q_@{j}) = exp(c_@{j}) + exp(mout_@{j}) + exp(iout_@{j});

    exp(q_@{j}) = exp(a_@{j}) * ((mu_@{j})^(1/sigma_q) * (exp(y_@{j}))^((sigma_q-1)/sigma_q) + (1-mu_@{j})^(1/sigma_q) * (exp(m_@{j}))^((sigma_q-1)/sigma_q) )^(sigma_q/(sigma_q-1));
    exp(y_@{j}) = (alpha_@{j}^(1/sigma_y) * exp(k_@{j})^((sigma_y-1)/sigma_y) + (1-alpha_@{j})^(1/sigma_y) * exp(l_@{j})^((sigma_y-1)/sigma_y) )^(sigma_y/(sigma_y-1));
    
@#endfor

exp(c_util)=(0
        @#for j in 1:n_sectors
            +xi_@{j}^(1/sigma_c)*(exp(c_@{j}))^((sigma_c-1)/sigma_c)
        @#endfor
        )^(sigma_c/(sigma_c-1));

exp(l_util)=(
        @#for j in 1:n_sectors
            +exp(l_@{j})^((sigma_l+1)/sigma_l)
        @#endfor
        )^(sigma_l/(sigma_l+1));

exp(c_agg) = (
        @#for j in 1:n_sectors
            +exp(p_@{j})*exp(c_@{j})
        @#endfor
        );

exp(l_agg) = (
        @#for j in 1:n_sectors
            +exp(l_@{j})
        @#endfor
        );

exp(gdp_agg) = (
        @#for j in 1:n_sectors
            +exp(p_@{j})*(exp(q_@{j}) - exp(mout_@{j}))
        @#endfor
        );

exp(i_agg) = (
        @#for j in 1:n_sectors
            +exp(p_@{j})*exp(iout_@{j})
        @#endfor
        );

exp(k_agg) = (
        @#for j in 1:n_sectors
            +exp(pkss_@{j})*exp(k_@{j})
        @#endfor
        );

utility_intratemp = exp(c_util) - theta*(1/(1+eps_l^(-1))) * exp(l_util)^(1+eps_l^(-1));

end;


