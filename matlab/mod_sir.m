function y = mod_sir(x,xdata)
%% Numeric implementation  -- JULY 2020
%  Paper: On an Alternative Susceptible-Infected-Removed Epidemic Model in Discrete-time
%  Authors: Jorge A. Costa Jr.; Amanda C. Martinez; Jose C. Geromel
%  Code developed on MATLAB Version R14, see [9]

M = xdata(1);
nu = xdata(2);
Nh = xdata(3:end);

N = length(Nh)-1;

z0 = x(1:3)';
gamma = x(4:3+N);
alpha = x(4+N:3+2*N);

H = [0 0 1];
y = [];
z = z0;

for j = 1:N
    A = diag([1 1-alpha(j) 1]);
    G = gamma(j)*[-1; 1; 1];
    for k = Nh(j)+1:Nh(j+1)
        y = [y H*z];
        w = z(2)*((z(1)/M)^nu);
        z = A*z + G*w;
    end
end

