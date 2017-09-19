function [ za_evol_pulse_M  zb_evol_pulse_M] = z_evol_pulse_M( t, t0, za_i, zb_i, r1a, r2a, r1b , r2b, kba, kab, dwa, dwb, w1, c1, c2 )
%z_evol_pulse_M Evolution in time of "a" pool using Meissner et al., NMR
%in Biomed, 2015, p.1196
%   Detailed explanation goes here

r1p_shaped_eff = r1a + (r2a-r1a).*c1.*w1.^2./(w1.^2 + c2^2*dwa.^2);
r1p_shaped_ex = kab.*c1.*w1.^2 ./ (w1.^2 + (kba.*(kba + r2b) + kba.*dwb.^2 ./(kba+r2b))*c2^2);
r1p_shaped = r1p_shaped_eff + r1p_shaped_ex;
za_cw_ss = r1a./r1p_shaped;     % following Meissner eqn 11, dropping cos

za_evol_pulse_M = (za_i - za_cw_ss)*exp(-r1p_shaped*(t-t0)) + za_cw_ss;
zb_evol_pulse_M = za_evol_pulse_M;      % I think this is the key assumption leading to single exponential recovery during pause
end

