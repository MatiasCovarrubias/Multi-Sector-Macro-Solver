function fsolve_options = build_quiet_fsolve_options()
% BUILD_QUIET_FSOLVE_OPTIONS Standard quiet fsolve options for tests/fixtures

fsolve_options = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10, ...
    'MaxFunEvals', 10000000, 'MaxIter', 10000);
end
