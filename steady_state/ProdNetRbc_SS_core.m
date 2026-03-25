function [base_losses, StStval, econ] = ProdNetRbc_SS_core(sol_base, params, alpha, mu, xi, Gamma_M, Gamma_I, print_flag)

beta = params.beta;
delta = params.delta;
rho = params.rho;
eps_l = params.eps_l;
eps_c = params.eps_c;
phi = params.phi;
sigma_c = params.sigma_c;
sigma_m = params.sigma_m;
sigma_q = params.sigma_q;
sigma_y = params.sigma_y;
sigma_I = params.sigma_I;
sigma_l = params.sigma_l;
n_sectors = params.n_sectors;
Sigma_A = params.Sigma_A;
theta = params.theta;
ss_idx = get_steady_state_indices(n_sectors);

%% -------------------- Endogenous variables -------------------- %%

C          = exp(sol_base(ss_idx.c(1):ss_idx.c(2)));
L          = exp(sol_base(ss_idx.l(1):ss_idx.l(2)));
Pk         = exp(sol_base(ss_idx.pk(1):ss_idx.pk(2)));
Pm         = exp(sol_base(ss_idx.pm(1):ss_idx.pm(2)));
M          = exp(sol_base(ss_idx.m(1):ss_idx.m(2)));
Mout       = exp(sol_base(ss_idx.mout(1):ss_idx.mout(2)));
I          = exp(sol_base(ss_idx.i(1):ss_idx.i(2)));
Iout       = exp(sol_base(ss_idx.iout(1):ss_idx.iout(2)));
P          = exp(sol_base(ss_idx.p(1):ss_idx.p(2)));
Q          = exp(sol_base(ss_idx.q(1):ss_idx.q(2)));
Y          = exp(sol_base(ss_idx.y(1):ss_idx.y(2)));
C_util     = exp(sol_base(ss_idx.c_util));
L_util     = exp(sol_base(ss_idx.l_util));

%% -------------------- Model equations -------------------- %%

K = I./delta;
utility_intratemp = C_util-theta*(1/(1+eps_l^(-1)))*L_util^(1+eps_l^(-1));
Pdef = (C_util * xi./C).^(1/sigma_c);
Lsup = theta*L_util^(eps_l^(-1)) * (L/L_util).^(1/sigma_l);

MPLmod = P .* ((mu).*Q./Y).^(1/sigma_q) .* ((1-alpha).*Y./L).^(1/sigma_y);
MPKmod = (beta./(1-beta*(1-delta))).*(P.*((mu).*Q./Y).^(1/sigma_q) .* (alpha.*Y./K).^(1/sigma_y));
Pmdef = (Gamma_M.'* (P.^(1-sigma_m))).^(1/(1-sigma_m)) ;
Mmod = (1-mu) .* (Pm./P).^(-sigma_q) .* Q  ;
Moutmod = P.^(-sigma_m) .* (Gamma_M*((Pm.^(sigma_m)).*M)); 
Pkdef =(Gamma_I.'* (P.^(1-sigma_I))).^(1/(1-sigma_I));
Ioutmod = P.^(-sigma_I) .* (Gamma_I* ((Pk.^(sigma_I)).*I)) ;
Qrc = C+Mout+Iout;
Qdef = (((mu).^(1/sigma_q).*Y.^((sigma_q-1)/sigma_q) + (1-mu).^(1/sigma_q).*M.^((sigma_q-1)/sigma_q) ).^(sigma_q/(sigma_q-1)));
Ydef = ((alpha.^(1/sigma_y).*K.^((sigma_y-1)/sigma_y) + (1-alpha).^(1/sigma_y).*L.^((sigma_y-1)/sigma_y) ).^(sigma_y/(sigma_y-1)));
C_util_def = ( (xi.^(1/sigma_c))' * (C.^((sigma_c-1)/sigma_c)) )^(sigma_c/(sigma_c-1));
L_util_def = (sum(L.^((sigma_l+1)/sigma_l)))^(sigma_l/(sigma_l+1));

Cagg_exp = P'*C;
Lagg_sum = sum(L);
GDPagg_exp = P'*(Q-Mout);
Iagg_exp = P'*Iout;
Kagg_price = Pk'*K;

C_loss = P./Pdef - 1;
L_loss = Lsup ./MPLmod - 1;
K_loss = Pk./MPKmod - 1;
Pm_loss = Pm./Pmdef - 1;
M_loss = M./Mmod - 1;
Mout_loss = Mout./Moutmod - 1;
Pk_loss = Pk./Pkdef - 1;
Iout_loss = Iout./Ioutmod - 1;
Qrc_loss = Q./Qrc - 1;
Qdef_loss = Q./ Qdef - 1;
Ydef_loss = Y./ Ydef - 1;
C_util_loss = C_util/C_util_def - 1;
L_util_loss = L_util/L_util_def - 1;

base_losses = zeros(ss_idx.n_base, 1);
base_losses(ss_idx.c(1):ss_idx.c(2)) = C_loss;
base_losses(ss_idx.l(1):ss_idx.l(2)) = L_loss;
base_losses(ss_idx.pk(1):ss_idx.pk(2)) = K_loss;
base_losses(ss_idx.pm(1):ss_idx.pm(2)) = Pm_loss;
base_losses(ss_idx.m(1):ss_idx.m(2)) = M_loss;
base_losses(ss_idx.mout(1):ss_idx.mout(2)) = Mout_loss;
base_losses(ss_idx.i(1):ss_idx.i(2)) = Pk_loss;
base_losses(ss_idx.iout(1):ss_idx.iout(2)) = Iout_loss;
base_losses(ss_idx.p(1):ss_idx.p(2)) = Qrc_loss;
base_losses(ss_idx.q(1):ss_idx.q(2)) = Qdef_loss;
base_losses(ss_idx.y(1):ss_idx.y(2)) = Ydef_loss;
base_losses(ss_idx.c_util) = C_util_loss;
base_losses(ss_idx.l_util) = L_util_loss;

%% Print Section %%

if print_flag
    disp(' ');

    disp('Analytic steady state');

    disp(' ');

    disp ('--- Prices ---')
    disp('Sectoral Price Indices (P_j): ');
    disp(P)
    disp('Capital Price Indices (P.^k_j): ');
    disp(Pk)
    disp('Intermediate Price Indices (P.^m_j): ');
    disp(Pm)

    disp('--- Quantities ---')
    disp('Sectoral Gross Output (Q_j): ')
    disp(Q)
    disp('Sectoral Value Added (Y_j): ')
    disp(Y)
    disp('Sectoral Intermediates (M_j): ')
    disp(M)
    disp('Sectoral Investment (I_j): ')
    disp(I)
    disp('Sectoral Labor (L_j): ')
    disp(L)
    disp('Sectoral Capital (K_j): ')
    disp(K)
    
    disp('--- Shares over nominal consumption  ---')
    disp('Sectoral Gross Output (P_jQ_j/PC): ')
    disp(P.*Q/(P'*C))
    disp('Sectoral Value Added (P_jY_j/PC): ')
    disp(P.*Y/(P'*C))
    disp('Sectoral Intermediates (Pm_jM_j/PC): ')
    disp(Pm.*M/(P'*C))
    disp('Sectoral Investment (Pk_jI_j/PC): ')
    disp(Pk.*I/(P'*C))
    disp('Sectoral Capital (Pk_jK_j/PC): ')
    disp(Pk.*K/(P'*C))
    
    disp('--- Aggregate Quantities ---')
    disp(['Aggregate Consumption Expenditure (C): ',num2str(Cagg_exp)])
    disp(['Aggregate Labor Headcount (L): ',num2str(Lagg_sum)])
    disp(['Aggregate GDP Expenditure (GDP): ',num2str(GDPagg_exp)])
    disp(['Aggregate Investment Expenditure (I): ',num2str(Iagg_exp)])
    disp(['Aggregate Capital (K): ',num2str(Kagg_price)])
    
    disp(['Intratemporal Utility: ',num2str(utility_intratemp)])
    
    disp('--- Aggregate Shares over nominal consumption  ---')
    disp('Sectoral Gross Output (PQ/PC): ')
    disp(P'*Q/(P'*C))
    disp('Aggregate Share Value Added (PY/PC): ')
    disp(P'*Y/(P'*C))
    disp('Aggregate Share of Intermediates (Pm M/PC): ')
    disp(Pm'*M/(P'*C))
    disp('Aggregate Share of Investment (Pk I/PC): ')
    disp(Pk'*I/(P'*C))
    disp('capital over Consumption (Pk K/PC): ')
    disp(Pk'*K/(P'*C))
    
    
end

%% Return Output

StStval.parameters    = struct('parn_sectors', n_sectors, 'parbeta', beta, 'pareps_c', eps_c, 'pareps_l', eps_l, 'parphi', phi, 'partheta', theta, ...
    'parsigma_c', sigma_c, 'parsigma_m', sigma_m, 'parsigma_q', sigma_q, 'parsigma_y', sigma_y,'parsigma_I', sigma_I, 'parsigma_l', sigma_l, ...
    'paralpha', alpha, 'pardelta', delta, 'parmu', mu, 'parrho', rho, 'parxi', xi, ...
    'parGamma_I', Gamma_I, 'parGamma_M', Gamma_M, 'parSigma_A', Sigma_A);
StStval.policies_ss = [ ...
    sol_base(ss_idx.c(1):ss_idx.y(2)); ...
    log(C_util); ...
    log(L_util); ...
    log(Cagg_exp); ...
    log(Lagg_sum); ...
    log(GDPagg_exp); ...
    log(Iagg_exp); ...
    log(Kagg_price); ...
    log(utility_intratemp)];
StStval.endostates_ss = log(K);
StStval.C_ss = Cagg_exp;
StStval.L_ss = Lagg_sum;
StStval.GDP_ss = GDPagg_exp;
StStval.I_ss = Iagg_exp;
StStval.K_ss = Kagg_price;
StStval.utility_intratemp_ss = utility_intratemp;

%% Intermediate variables for moment matching

econ.C = C;
econ.L = L;
econ.K = K;
econ.Iout = Iout;
econ.Mout = Mout;
econ.Q = Q;
econ.P = P;
econ.Pm = Pm;
econ.Pk = Pk;
econ.Ydef = Ydef;
econ.Qdef = Qdef;

end
