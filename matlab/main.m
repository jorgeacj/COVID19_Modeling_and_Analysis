%% MAIN ROUTINE
%% Numeric implementation  -- JULY 2020
%  Paper: On an Alternative Susceptible-Infected-Removed Epidemic Model in Discrete-time
%  Authors: Jorge A. Costa Jr.; Amanda C. Martinez; Jose C. Geromel
%  Code developed on MATLAB Version R14, see [9]
%  This program is a free software: you can redistribute it and/or modify it under the terms of the MIT license.

%% Initialization
clear all
close all

%% Data loading
% The following dataset was obtained from World Health Organization [13]
% for the pandemics in Italy and Brazil, respectively. 
%
%%%%%%%%%%%%%%%%%%%%%%  Italy
M = 60e+6;
am = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 17 79 132 229 322 400 ...
      650 888 1128 1689 2036 2502 3089 3858 4636 5883 7375 9172 10149 ...
      12462 15113 17660 21157 23980 27980 31506 35713 41035 47021 53578 ...
      59138 63927 69176 74386 80539 86498 92472 97689 101739 105792 ...
      110574 115242 119827 124632 128948 132547 135586 139422 143626 ...
      147577 152271 156363 159516 162488 165155 168941 172434 175925 ...
      178972 181228 183957 187327 189973 192994 195351 197675 199414 ...
      201505 203591 205463 207428 209328 210717 211938 213013 214457 ...
      215858 217185 218268 219070 219814 221216 222104 223096 223885 ...
      224760 225435 225886 226699 227364 228006 228658 229327 229858 ...
      230158 230555 231139 231732 232248 232664 233019 233197];
%
%%%%%%%%%%%%%%%%%%%%%%%  Brazil
% M = 210e+6;
% am = [1 1 1 1 2 2 2 2 3 8 13 13 25 25 34 52 77 98 121 200 234 291 428 ...
%       621 904 1128 1546 1891 2201 2433 2915 3417 3904 4256 4579 5717 ...
%       6836 7910 9056 10278 11130 12056 13717 15927 17857 19638 20727 ...
%       22169 23430 25262 28320 30425 33682 36599 38654 40581 43079 45757 ...
%       49492 52995 58509 61888 66501 71886 78162 85380 91589 96559 101147 ...
%       107780 114715 125218 135106 145328 155939 162699 168331 177589 ...
%       188974 202918 218223 233142 241080 254220 271628 291579 310087 ...
%       330890 347398 363211 374898 391222 411821 438238 465166 498440 ...
%       514849 526447 555383 584016 614932 645762 672837 691758 707412 ...
%       739503 772416 802828 828810 850514 867624 888271 923189 955377 ...
%       978142 1032913 1067579 1085038 1106470 1145906 1188631 1228114 ...
%       1274974 1313667 1344143 1368195 1402041 1448753 1496858 1539081 ...
%       1577004 1603055 1623284 1668589 1713160 1755779 1800827 1839850 ...
%       1864681 1884967 1926824 1966748 2012151 2046328 2074860 2098389]; 
%
%% Algorithm and model selection
%
% Algorithm selection:
% alg = 1 - Time-invariant parameter optimization
% alg = 2 - Time-varying parameter optimization
% alg = 3 - Sequential forward optimization
%
alg = 3;

% Model selection:
% nu = 1 - Classical SIR model
% nu = 2 - Alternative SIR model
%
nu = 2;

%% Data selectioning for parameter identification
if alg == 1
    Nh = [0 length(am)];
elseif alg == 2 || alg == 3
    Nh = [0 42 49 77 98 length(am)];
else
    disp('Error in algorithm selection. Please select a valid algorithm.')
    return
end

%% Number of time subintervals and simulation horizon
N = length(Nh)-1;
Nsim = 124;

%% Optimization using function LSQCURVEFIT
options = optimset('Display','final',...
    'MaxFunEvals',20000,...
    'MaxIter',2000,...
    'TolFun',1e-8,...
    'TolX',1e-8);

if alg == 1 || alg == 2
    xdata = [M nu Nh];
    p0 = [M-am(1) am(1) am(1) (nu/2)*ones(1,N) (1/2)*ones(1,N)];
    lb = [0 0 0 0*ones(1,N) (1/10)*ones(1,N)];
    ub = [M M M nu*ones(1,N) ones(1,N)];
    [p_min, e_rr] = lsqcurvefit(@mod_sir,p0,xdata,am,lb,ub,options);
elseif alg == 3
    z0 = [M-am(1) am(1) am(1)];
    p_min = [];
    e_min = [];
    p0 = [nu/2 1/2];
    lb = [0 1/10];
    ub = [nu 1];
    for j = 1:N
        xdata = [z0 M nu Nh(j+1)-Nh(j)];
        ydata = am(Nh(j)+1:Nh(j+1));
        [p,e_rr] = lsqcurvefit(@mod_sir_sfo,p0,xdata,ydata,lb,ub,options);
        p_min = [p_min p'];
        e_min = [e_min e_rr];
        p0 = p;
        Z = mod_sir_sfo_sim(p0,xdata);
        z0 = Z(:,end)';
    end
end

%% Printing metrics
if alg == 1 || alg == 2
    z0 = p_min(1:3)
    gamma = p_min(4:3+N)
    alpha = p_min(4+N:3+2*N)
    R0 = gamma./alpha
    e_min = (1/2)*log10(e_rr/Nh(end))
%% The constraint s0 + i0 =< M is rarely infeasible. For this reason, it is not imposed but tested
    if (z0(1) + z0(2) > M)
        disp('INVALID SOLUTION!!!')
    end
%% If this condition is verified the calculated solution is useless!!!
elseif alg == 3
    z0 = [M-am(1) am(1) am(1)]
    gamma = p_min(1,:)
    alpha = p_min(2,:)
    R0 = gamma./alpha
    e_min = (1/2)*log10(sum(e_min)/Nh(end))
end

%% Testing if the solution adheres the official data
% Simulating the model with the optimal parameters for the simulation horizon
p = [z0 gamma alpha];
xdata = [M nu Nh(1:end-1) Nsim];
a = mod_sir(p,xdata);

%% Calculating R0 as function of time
r0 = [];
for j = 1:length(Nh)-1
    r0 = [r0 R0(j)*ones(1,(Nh(j+1) - Nh(j)))];
end

%% Plotting results
% Defining time basis
Time = 0:(Nh(end) - 1);
Timesim = 0:(Nsim - 1);

%% Daily confirmed cases
figure(1)
subplot(2,1,1), plot(Time,[1 diff(a(1:Nh(end)))],'-b',Time,[1 diff(am)],':ro')
xlabel('x'), ylabel('y')
subplot(2,1,2), plot(Timesim,[1 diff(a(1:Nsim))],'-b')
xlabel('x'), ylabel('y')

%% Total of confirmed cases
figure(2)
subplot(2,1,1), plot(Time,a(1:Nh(end)),'-b',Time,am,':ro')
xlabel('x'), ylabel('y')
subplot(2,1,2), plot(Timesim,a(1:Nsim),'-b')
xlabel('x'), ylabel('y')

%% Total of confirmed cases (log10 scale)
figure(3)
subplot(2,1,1), plot(Time,log10(a(1:Nh(end))),'-b',Time,log10(am),':ro'), grid
xlabel('x'), ylabel('y')
subplot(2,1,2), plot(Timesim,log10(a(1:Nsim)),'-b'), grid
xlabel('x'), ylabel('y')

%% R0 as a function of time
figure(4)
plot(Time,ones(size(Time)),'--r',0,0,'b')
hold
stairs(Time,r0,'-b')
xlabel('x'), ylabel('r')

%% Daily confirmed cases and total confirmed cases
figure(5)
subplot(2,1,1), plot(Timesim,[1 diff(a(1:Nsim))],'-b',Time,[1 diff(am)],':r.')
xlabel('x'), ylabel('n')
subplot(2,1,2), plot(Timesim,a(1:Nsim),'-b',Time,am,':r.')
xlabel('x'), ylabel('a')
