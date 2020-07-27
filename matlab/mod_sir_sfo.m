function y = mod_sir_sfo(x,xdata)
%% Numeric implementation  -- JULY 2020
%  Paper: On an Alternative Susceptible-Infected-Removed Epidemic Model in Discrete-time
%  Authors: Jorge A. Costa Jr.; Amanda C. Martinez; Jose C. Geromel
%  Code developed on MATLAB Version R14, see [9]
%  This program is a free software: you can redistribute it and/or modify it under the terms of the MIT license.

z0 = xdata(1:3)';
M = xdata(4);
nu = xdata(5);
Nh = xdata(6);

gamma = x(1);
alpha = x(2);

H = [0 0 1];
y = [];
z = z0;

A = diag([1 1-alpha 1]);
G = gamma*[-1; 1; 1];
for k = 1:Nh
    y = [y H*z];
    w = z(2)*((z(1)/M)^nu);
    z = A*z + G*w;
end
