function Z = mod_sir_sfo_sim(x,xdata)
%% Numeric implementation  -- JULY 2020
%  Paper: On an Alternative Susceptible-Infected-Removed Epidemic Model in Discrete-time
%  Authors: Jorge A. Costa Jr.; Amanda C. Martinez; Jose C. Geromel
%  Code developed on MATLAB Version R14, see [9]

z0 = xdata(1:3)';
M = xdata(4);
nu = xdata(5);
Nh = xdata(6);

gamma = x(1);
alpha = x(2);

H = [1 1 1];
Z = [];
z = z0;

A = diag([1 1-alpha 1]);
G = gamma*[-1; 1; 1];
for k = 1:Nh
    w = z(2)*((z(1)/M)^nu);
    z = A*z + G*w;
    Z = [Z z];
end
