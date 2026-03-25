function run_dynare_mod(dynare_folder, mod_name)
    current_dir = pwd;
    cd(dynare_folder);
    try
        evalin('base', ['dynare ' mod_name]);
    catch ME
        cd(current_dir);
        rethrow(ME);
    end
    cd(current_dir);
end
