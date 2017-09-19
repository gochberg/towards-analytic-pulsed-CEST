function [ z_returned ] = z_evol_pulse_rot_full_3_derivedTorrey( t, z_i, r1a, r2a, r1b , r2b, kba, kab, dwa, dwb, w1)
% z_evol_pulse_rot_full_3_derivedTorrey
% use derived short time and long time solns as in
% z_evol_pulse_rot_full_3_derived, but with Torrey approximation for the
% eigenvalues, which are then used both in the rates and in the amplitude
% calculations.  Also, use large-shift approximation for R1p rate and
% corresponding steady state value R1a/R1p.


%z_evol_pulse_rot_full_3_derived: based on z_evol_pulse_rot_full_3.m
% key change from previous version: base evolution on derived version,
% using short time coupling to za_ic and long time coupling to za_ss.  See
% Bloch_coupled_eigenvalues_z3.m with short_time_approach = 1,
% long_time_approach = 1,
%  make_zb_ic_match = 1;

% z_evol_pulse_rot_full_3:
% Evolution in time of "a" and "b"
% pools.  See lab notebook 26, p.107 
% For conection to shaped pulse, see constant avg w1 and w1^2 calcs on notebook 26, p.127
%   Gives function gives analytic approximation of evolution during hard pulse with small
%   b pool


% a pool: just R1p, like Meissner, but since using pause_hard_pause approach, use c1=1, c2=1 during this pulse segment.  This should
% be very close to Meissner.  He does R1p = R1a + c1* ( exch effects with
% w1avg replaced w1avg/c2.  We'll do time distinct regions, and we will
% preserve avg w1 over tp duration.  

% In addition, I'm adding rotation and
% rapid decay to the b pool during the pulse


za_i = z_i(:,1);
zb_i = z_i(:,2);
t0 = t(1);
ab_offset = dwb-dwa;
pbpa = kab/kba;

% prealocate, matching code in Bloch_coupled_eigenvalues_z4.m, but with no
% j,k,l
z_an_ss = zeros([1 2]); 
z_an_comp_red2_amp = zeros([2 6]);  % [a or b,  R1p or R1pfast or Rabi_a cos or Rabi_a sin or Rabi_b cos or Rabi_b sin components]
z_an_comp_red2= zeros([length(t) 2 6]);
z_an = zeros([length(t) 2]);

% rates used for a pool (from Bloch_coupled_eigenvalues_z3.m)
    rabi_rate_a = sqrt( w1^2 + dwa^2);
    rabi_rate_b = sqrt( w1^2 + dwb^2);
    Reff = r1a + (r2a - r1a) * w1^2 / (w1^2 + dwa^2);
    gamma = 2*sqrt( (kba+r2b)*w1^2/kba + (kba+r2b)^2 );
    Rex = kab * ab_offset^2 * w1^2 / ((w1^2 + dwa^2)*(1/4 * gamma^2 + dwb^2)) + ...
        pbpa * r2b * w1^2 / (1/4 * gamma^2 + dwb^2) + ...
        kab * w1^2 / (w1^2 + dwa^2) * r2b*(r2b + kba) / (1/4 * gamma^2 + dwb^2); 
    Rex_LS = kab*w1^2/(w1^2 + kba*(r2b+kba) + kba*dwb^2/(r2b+kba)); % Meissner LS calcs email
    R1p_LS = Reff + Rex_LS;
    R1p = Reff + Rex;
    
% additional rates and amps used for b pool (from Bloch_coupled_eigenvalues_z3.m)

    % rates use Torrey with r1b->r1b+kba and r2b-> r2b+kba. with the exception of R1p, calculated above.  amps (below)
    % use full solution, derived at short time (ST).
    R2p_b_Torrey = r2b - 0.5*(r2b-r1b)/(1+(dwb/w1)^2);     % Torrey, eqn 59b        
    R2p_b = kba + R2p_b_Torrey;  % my guess for b pool (relation to R2rho_a when exch?? -> see Garwood papers)
    R1p_fast_b_Torrey = (r2b + r1b*(dwb/w1)^2)/(1+(dwb/w1)^2);  % Torrey eqn 59a, with same assumptions listed above, 
    R1p_fast = kba + R1p_fast_b_Torrey;
    
        % short time approach.  I'll use the Morris/Torrey notation that
        % normailzes variables by w1.  It doesn't make any difference in
        % A,B,C amplitudes, and it matches my notebook and Mathematica
        % work.  I'll divide the amps by Mb0, to match zb.
        
        kST = kba/w1;
        alphaST = r1b/w1;
        betaST = r2b/w1;
        deltaST = dwb/w1;
        deltaAST = dwa/w1;
        
        DST_offset = (alphaST*(1 + deltaAST^2)*(betaST^2 + deltaST^2 + 2*betaST*kST + kST^2) )/ ...
            ((1 + deltaAST^2)*(betaST + kST + (alphaST + kST)*(deltaST^2 + (betaST + kST)^2)));
        DST_za_coef = deltaAST*kST*(deltaST + deltaST^2*deltaAST + deltaAST*(betaST + kST)^2)/ ...
            ((1 + deltaAST^2)*(betaST + kST + (alphaST + kST)*(deltaST^2 + (betaST + kST)^2)));
        
        % put in Torrey type approx for rates when calculating amps for
        % 'derivedTorrey' case.  aST->R1p_fast, bST->R2p_b, 
        % sST->rabi_rate_b, gammaST->gammaST_Torrey
        aST_Torrey = R1p_fast/w1;
        bST_Torrey = R2p_b/w1;
        sST_Torrey = rabi_rate_b/w1;
        gammaST_Torrey = aST_Torrey*( (bST_Torrey - aST_Torrey)^2 + sST_Torrey^2 );
        DST_offset_Torrey = DST_offset;
        DST_za_coef_Torrey = DST_za_coef;
        AST_offset_Torrey = -((alphaST*(aST_Torrey^2 + betaST^2 + deltaST^2 + 2*betaST*kST + kST^2 - 2*aST_Torrey*(betaST + kST)))/gammaST_Torrey); 
        AST_za_coef_Torrey =  - ((deltaST*deltaAST*kST + deltaST^2*deltaAST^2*kST + deltaAST^2*kST*(-aST_Torrey + betaST + kST)^2))/ ...
            ((1 + deltaAST^2)*gammaST_Torrey);
        AST_zb_coef_Torrey =  - ((-(aST_Torrey*deltaST^2) - aST_Torrey*deltaST^2*deltaAST^2 - aST_Torrey*(-aST_Torrey + betaST + kST)^2 - aST_Torrey*deltaAST^2*(-aST_Torrey + betaST + kST)^2))/ ...
            ((1 + deltaAST^2)*gammaST_Torrey);
        BST_offset_Torrey =  -AST_offset_Torrey-DST_offset_Torrey;
        BST_za_coef_Torrey = -AST_za_coef_Torrey-DST_za_coef_Torrey;
        BST_zb_coef_Torrey = 1 -AST_zb_coef_Torrey; 
        CST_offset_Torrey = alphaST + aST_Torrey*AST_offset_Torrey + bST_Torrey*BST_offset_Torrey;
        CST_za_coef_Torrey = (kST*deltaAST^2/(1 + deltaAST^2) ) + aST_Torrey*AST_za_coef_Torrey + bST_Torrey*BST_za_coef_Torrey;
        CST_zb_coef_Torrey = -(alphaST + kST) + aST_Torrey*AST_zb_coef_Torrey + bST_Torrey*BST_zb_coef_Torrey;

% za (just R1p + ss)    

        z_an_ss(1) = (dwa^2/(dwa^2+w1^2))*r1a./R1p; % Used cos2theta, even though likely adiabatic transition.  Correct?
        z_an_comp_red2_amp(1,1)  = za_i - z_an_ss(1); % z_an_comp_red2_amp(a or b, component #, j,k,l)
                            % z_an_ss has only two components, z_exact_i has 6.
                            % Amps for components 2-6 are 0, i.e. just R1p + ss 
                            
        % put in time dependence
        z_an_comp_red2(:,1,1) = exp(-R1p*(t-t0)) * z_an_comp_red2_amp(1,1);  
        
            % ss + sum of components
        z_an(:,1) = ones(size(t))*z_an_ss(1) + sum(squeeze(z_an_comp_red2(:,1,:)), 2); 
         
% zb (R1p + R1p_fast + Rabi_b cos + Rabi_b sin + ss, ignoring Rabi_a cos, and Rabi_a sin)

        % see analytic_CEST_zb_ss.nb and lab notebook 27, p.29.  Solve Bloch when za=za_ss by Moritz.  Use R1p for 'derivedTorrey' case   
        z_an_ss(2) = (dwa^2*(kba*r1a + r1b*R1p)*(dwb^2 + (kba + r2b)^2) + dwa*dwb*kba*r1a*w1^2 + ...
            r1b*R1p*(dwb^2 + (kba + r2b)^2)*w1^2)/(R1p*(dwa^2 + w1^2)* ...
            (dwb^2*(kba + r1b) + (kba + r2b)*(kba^2 + r1b*r2b + kba*(r1b + r2b) + w1^2)));
        
        % for amp calcs below, see analytic_CEST_short_times2.nb
        % (Mathematica notebook).  Use Torrey type approx for rates in amp
        % calc for 'derivedTorrey case'.
        z_an_comp_red2_amp(2,2) = zb_i*AST_zb_coef_Torrey + za_i*AST_za_coef_Torrey + AST_offset_Torrey;  % R1p_fast
        z_an_comp_red2_amp(2,5) = zb_i*BST_zb_coef_Torrey + za_i*BST_za_coef_Torrey + BST_offset_Torrey;  % Rabi_b cos
        z_an_comp_red2_amp(2,6) = zb_i*CST_zb_coef_Torrey/sST_Torrey + za_i*CST_za_coef_Torrey/sST_Torrey + CST_offset_Torrey/sST_Torrey ; % Rabi b sin
   
        % make ic match, but still ignoring rabi_a cos component (component 3):
        z_an_comp_red2_amp(2,1) = zb_i - z_an_comp_red2_amp(2,2) - z_an_comp_red2_amp(2,5)-z_an_ss(2);
  
        % put in time dependence (haven't modified all the rates to
            % match Morris.)
        z_an_comp_red2(:,2,1) = exp(-R1p*(t-t0)) * z_an_comp_red2_amp(2,1);  
        z_an_comp_red2(:,2,2) = exp(-R1p_fast*(t-t0)) * z_an_comp_red2_amp(2,2);
        z_an_comp_red2(:,2,5) = exp(-R2p_b*(t-t0)) .* cos( rabi_rate_b* (t-t0) ) * ...
            z_an_comp_red2_amp(2,5); % Rabi_b cos component
        z_an_comp_red2(:,2,6) = exp(-R2p_b*(t-t0)) .* sin( rabi_rate_b* (t-t0) ) * ...
            z_an_comp_red2_amp(2,6); % Rabi_b sin component

                % ss + sum of components
        z_an(:,2) = ones(size(t))*z_an_ss(2) + sum(squeeze(z_an_comp_red2(:,2,:)), 2); 

z_returned = z_an;

end

