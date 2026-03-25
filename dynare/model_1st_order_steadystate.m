function [ys, params, check] = model_1st_order_steadystate(ys, exo, M_, options_)
check = 0;
params = M_.params;
load('ModStruct_temp.mat', 'policies_ss', 'k_ss', 'parn_sectors');
policies_ss = real(policies_ss);
k_ss = real(k_ss);
n = parn_sectors;
ys(1:n) = k_ss;
ys(n+1:2*n) = 0;
ys(2*n+1:3*n) = policies_ss(1:n);
ys(3*n+1:4*n) = policies_ss(n+1:2*n);
ys(4*n+1:5*n) = policies_ss(2*n+1:3*n);
ys(5*n+1:6*n) = policies_ss(3*n+1:4*n);
ys(6*n+1:7*n) = policies_ss(4*n+1:5*n);
ys(7*n+1:8*n) = policies_ss(5*n+1:6*n);
ys(8*n+1:9*n) = policies_ss(6*n+1:7*n);
ys(9*n+1:10*n) = policies_ss(7*n+1:8*n);
ys(10*n+1:11*n) = policies_ss(8*n+1:9*n);
ys(11*n+1:12*n) = policies_ss(9*n+1:10*n);
ys(12*n+1:13*n) = policies_ss(10*n+1:11*n);
ys(13*n+1) = policies_ss(11*n+1);
ys(13*n+2) = policies_ss(11*n+2);
ys(13*n+3) = policies_ss(11*n+3);
ys(13*n+4) = policies_ss(11*n+4);
ys(13*n+5) = policies_ss(11*n+5);
end
