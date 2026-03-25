function write_val_block(fid, block_type, n_sectors_macro, a_expr, k_expr)
% write_val_block  Write a Dynare initval/endval block for the production network model.
%
%   write_val_block(fid, block_type, n_sectors_macro)
%   write_val_block(fid, block_type, n_sectors_macro, a_expr)
%   write_val_block(fid, block_type, n_sectors_macro, a_expr, k_expr)
%
%   Inputs:
%     fid              - file identifier from fopen
%     block_type       - 'initval' or 'endval'
%     n_sectors_macro  - Dynare macro variable name for number of sectors
%                        (e.g. 'n_sectors')
%     a_expr           - (optional, default '0') expression for the TFP
%                        state variable a_@{j}. Use 'shocksim_0(@{j},1)'
%                        when the initval block must reflect a non-zero
%                        initial TFP level.
%     k_expr           - (optional, default 'k_ss(@{j})') expression for
%                        the capital state variable k_@{j}.

if nargin < 4
    a_expr = '0';
end

if nargin < 5
    k_expr = 'k_ss(@{j})';
end

fprintf(fid, '%s;\n', block_type);
fprintf(fid, '@#for j in 1:%s\n', n_sectors_macro);
fprintf(fid, '    e_@{j}=0;\n');
fprintf(fid, '    a_@{j}=%s;\n', a_expr);
fprintf(fid, '    k_@{j}=%s;\n', k_expr);
fprintf(fid, '    c_@{j}=policies_ss(@{j});\n');
fprintf(fid, '    l_@{j}=policies_ss(parn_sectors+@{j});\n');
fprintf(fid, '    pk_@{j}=policies_ss(2*parn_sectors+@{j});\n');
fprintf(fid, '    pm_@{j}=policies_ss(3*parn_sectors+@{j});\n');
fprintf(fid, '    m_@{j}=policies_ss(4*parn_sectors+@{j});\n');
fprintf(fid, '    mout_@{j}=policies_ss(5*parn_sectors+@{j});\n');
fprintf(fid, '    i_@{j}=policies_ss(6*parn_sectors+@{j});\n');
fprintf(fid, '    iout_@{j}=policies_ss(7*parn_sectors+@{j});\n');
fprintf(fid, '    p_@{j}=policies_ss(8*parn_sectors+@{j});\n');
fprintf(fid, '    q_@{j}=policies_ss(9*parn_sectors+@{j});\n');
fprintf(fid, '    y_@{j}=policies_ss(10*parn_sectors+@{j});\n');
fprintf(fid, '@#endfor\n');
fprintf(fid, 'c_util = policies_ss(11*parn_sectors+1);\n');
fprintf(fid, 'l_util = policies_ss(11*parn_sectors+2);\n');
fprintf(fid, 'c_agg = policies_ss(11*parn_sectors+3);\n');
fprintf(fid, 'l_agg = policies_ss(11*parn_sectors+4);\n');
fprintf(fid, 'gdp_agg = policies_ss(11*parn_sectors+5);\n');
fprintf(fid, 'i_agg = policies_ss(11*parn_sectors+6);\n');
fprintf(fid, 'k_agg = policies_ss(11*parn_sectors+7);\n');
fprintf(fid, 'utility_intratemp = policies_ss(11*parn_sectors+8);\n');
fprintf(fid, 'end;\n\n');

end
