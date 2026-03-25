function [oo, M, options] = ensure_dynare_solution(order, dynare_folder, verbose)
    cache = get_dynare_solution_cache(order);
    if verbose, fprintf('Solving order=%d (%s)...\n', order, cache.mod_name); end
    tic;
    run_dynare_mod(dynare_folder, cache.mod_name);
    if verbose, fprintf('Order=%d solved (%.1fs)\n', order, toc); end

    oo = read_base_workspace_var('oo_');
    M = read_base_workspace_var('M_');
    options = read_base_workspace_var('options_');
end

function cache = get_dynare_solution_cache(order)
cache = struct();

if order == 1
    cache.oo_var = 'oo_1st_';
    cache.M_var = 'M_1st_';
    cache.options_var = 'options_1st_';
    cache.mod_name = 'model_1st_order';
    cache.required_dr_field = 'ghx';
else
    cache.oo_var = 'oo_2nd_';
    cache.M_var = 'M_2nd_';
    cache.options_var = 'options_2nd_';
    cache.mod_name = 'stoch_simul_2ndOrder';
    cache.required_dr_field = 'ghxx';
end
end
