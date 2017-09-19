function [ z_evol_pause ] = z_evol_pause3( t, z_i, r1a, r1b, kba, kab )
%z_evol_pause Evolution in time of "a" and "b
% pools.  See lab notebook 26, p.117
%   Detailed explanation goes here

% treat t like ode45 does.  First value is t0 (with value z_i).  Solution
% goes to t(end)

t0 = t(1);
za_i = z_i(1);
zb_i = z_i(2);

% using first order approx:
z_evol_pause(:,1) = kab/kba*(za_i - zb_i)*exp(-kba*(t-t0)) + (za_i - kab/kba*(za_i - zb_i) - 1)*exp(-r1a*(t-t0)) + 1; % za
z_evol_pause(:,2) = -1*(za_i - zb_i)*exp(-kba*(t-t0)) + (za_i - 1)*exp(-r1a*(t-t0)) + 1; % zb

end

