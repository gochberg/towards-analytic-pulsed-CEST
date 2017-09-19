% figures4to9.m
%
% produced by D.F. Gochberg for paper:
% D.F. Gochberg, M.D. Does, Z. Zu, C.L. Lankford.  Towards an analytic
% solution for pulsed CEST.  NMR in Biomedicine 2017
%
% You are free to use this code for non-commercial purposes, but please
% cite the above manuscript if you use the code, or parts thereof, to help
% produce a manuscript or presentation figure, as appropriate.  Thanks.
%
% Based on: Bloch_coupled_eigenvalues_z4_R2016b.m (Dan's not to himself)


clear all;

% cd to whichever directory contains this script, avoiding confusion 
% between called functions with same name but different directories


if(~isdeployed)
  	cd(fileparts(mfilename('fullpath')));
end

% PARAMETERS   ************************************************************

% sample parameteres
pbpa=0.004;
z_exact_i=[0;0;.9;0;0;.7]; % initial value;
r1a=1;
r2a=10;
r1b=1;
r2b=100;

kbas = [.0001 50 1000];
B0 = 9.4;
ab_offsets = [3.5 3.5 2.0]*(267.5*B0);     %ppm to rad/s.  Change ab_offsets in parallel with kbas
kba_disp = 2;

    % calculated
    n_exch = length(kbas);     

% sequence parameters
w1s = [1 3 6]*267.5;  % uT to rad/s
w1_disp = 1;

    % calculated
    n_powers = length(w1s);

% plot parameters
plot_n = 10;
n_offsets = 161; % use number such that ab_offsets land on point %141;
w1_offsets_ppm_start = -4.5;
w1_offsets_ppm_end = 0.5;
delta_from_b_for_plot_ppm = 0.4;
n_times = 600;

    % calculated
    w1_offsets_ppm = w1_offsets_ppm_start:(w1_offsets_ppm_end-w1_offsets_ppm_start)/(n_offsets-1):w1_offsets_ppm_end; 
    w1_offsets= w1_offsets_ppm*267.5*B0;    %ppm to rad/s
    for j=1:length(ab_offsets)
        [dummy, offset_res_disp(j)] = min(abs(w1_offsets + ab_offsets(j))); % use +, since different signs
        [dummy, offset_delta_disp(j)] = min(abs(w1_offsets + ab_offsets(j) + delta_from_b_for_plot_ppm*267.5*B0));
    end 
    offset_disp = offset_res_disp(1);
    

    t=(0:1.5*pi/min(w1s)/(n_times-1):1.5*pi/min(w1s))';  % need to start at t=0 for initial condition code (below) to work correctly




% VARIABLES ***************************************************************
 
% dimensionality
    % Used to lable dimensions when preallocating (so meaning
    % of each dimension is clear), but is not used to dynamically change any
    % actual calculations.
   
n_vector_elem = 6;  % [Max/Ma0 May/Ma0 za=Maz/Ma0 Mbx/Mb0 Mby/Mb0 zb=Mbz/Mb0]
n_eig_values = 6;
n_comp = 6; % # of components is dictated by n_eig_values
n_red2_comp = 6;  % 6 'reduced' components: R1p, R1p_fast, Rabi_a cos, Rabi_a sin, Rabi_b cos, Rabi_b sin 
n_an_vector_elem = 2;   % za or zb
n_an_comp = n_red2_comp; % n_red2_comp breaks down numeric results. n_an_comp are analytic components that use the same component breakdown.
                         % I use separate names since you could breakdown both the numeric and analytic into many different component combinations.

% preallocating 
    % analytic guesses for rates
rabi_rate_a = zeros([n_offsets,n_powers,n_exch]);
rabi_rate_b = zeros([n_offsets,n_powers,n_exch]);
Reff = zeros([n_offsets,n_powers,n_exch]);
Rex = zeros([n_offsets,n_powers,n_exch]);
Rex_LS = zeros([n_offsets,n_powers,n_exch]);
R1p = zeros([n_offsets,n_powers,n_exch]);
R1p_LS = zeros([n_offsets,n_powers,n_exch]);
R2p_b = zeros([n_offsets,n_powers,n_exch]);
R2p_a = zeros([n_offsets,n_powers,n_exch]);
R2p_a_Moran = zeros([n_offsets,n_powers,n_exch]);
R2p_a_Garwood = zeros([n_offsets,n_powers,n_exch]);
R1p_fast = zeros([n_offsets,n_powers,n_exch]);

    % analytic solution
z_an = zeros([n_times, n_an_vector_elem, n_offsets,n_powers,n_exch] );
z_an_ss = zeros([n_an_vector_elem, n_offsets,n_powers,n_exch] );   
z_an_comp_red2_amp = zeros([n_an_vector_elem, n_an_comp, n_offsets,n_powers,n_exch] ); % use red2 in name, since components match numeric z_comp_red2_amp
z_an_comp_red2 = zeros([n_times, n_an_vector_elem, n_an_comp, n_offsets,n_powers,n_exch] );  

    % numeric solution
lamda = zeros([n_eig_values,n_offsets,n_powers,n_exch]);
V = zeros([n_vector_elem,n_eig_values,n_offsets,n_powers,n_exch]);
D = zeros([n_eig_values,n_eig_values,n_offsets,n_powers,n_exch]);
z_exact = zeros( [n_times, n_vector_elem, n_offsets,n_powers,n_exch] ); % the 'z' here refers to normalized magnetization (like z-spectrum) of Mx, My, and Mz and not only Mz
z_exact_ss = zeros( [n_vector_elem, n_offsets,n_powers,n_exch] );
z_comp_weightings = zeros( [n_comp, n_offsets,n_powers,n_exch] );  
z_comp = zeros( [n_times, n_vector_elem, n_comp, n_offsets,n_powers,n_exch] );
z_comp_amp = zeros( [n_vector_elem, n_comp, n_offsets,n_powers,n_exch] );
z_comp_red2 = zeros( [n_times, n_vector_elem, n_red2_comp, n_offsets,n_powers,n_exch] );
z_comp_red2_amp = zeros( [n_vector_elem, n_red2_comp, n_offsets,n_powers,n_exch] );


% CALCULATIONS ************************************************************

% looping through kba/ab_offset, w1, and w1_offset

for l=1:n_exch
    kba = kbas(l);
    kab= kba*pbpa;
    ab_offset = ab_offsets(l);  % change in parallel with kba to match amides and amines
for k=1:n_powers
    w1 = w1s(k);
for j=1:n_offsets
    w1_offset = w1_offsets(j);
    dwa=w1_offset;
    dwb=w1_offset + ab_offset;
    
% analytic solution

    % analytic guesses for rates:  
    
    rabi_rate_a(j,k,l) = sqrt( w1^2 + dwa^2);
    rabi_rate_b(j,k,l) = sqrt( w1^2 + dwb^2);
        % see table 2 of Zaiss et al review 2013 for following:
    Reff(j,k,l) = r1a + (r2a - r1a) * w1^2 / (w1^2 + dwa^2);
    Reff_LS = r1a;
    gammaM = 2*sqrt( (kba+r2b)*w1^2/kba + (kba+r2b)^2 );
    Rex(j,k,l) = kab * ab_offset^2 * w1^2 / ((w1^2 + dwa^2)*(1/4 * gammaM^2 + dwb^2)) + ...
        pbpa * r2b * w1^2 / (1/4 * gammaM^2 + dwb^2) + ...
        kab * w1^2 / (w1^2 + dwa^2) * r2b*(r2b + kba) / (1/4 * gammaM^2 + dwb^2);
    Rex_LS(j,k,l) = kab*w1^2/(w1^2 + kba*(r2b+kba) + kba*dwb^2/(r2b+kba)); % Meissner LS calcs email
    R1p_LS(j,k,l) = Reff(j,k,l) + Rex_LS(j,k,l);    
    R1p(j,k,l) = Reff(j,k,l) + Rex(j,k,l);    
    R2p_b_Torrey = r2b - 0.5*(r2b-r1b)/(1+(dwb/w1)^2);     % Torrey, eqn 59b 
        % assumes, w1 >> R2 (and single pool)      
    R2p_b(j,k,l) = kba + R2p_b_Torrey;  
    R2p_a_Torrey = r2a - 0.5*(r2a-r1a)/(1+(dwa/w1)^2);     % Torrey, eqn 59b
    R2p_a(j,k,l) = R2p_a_Torrey + kab;  
    R1p_fast_b_Torrey = (r2b + r1b*(dwb/w1)^2)/(1+(dwb/w1)^2);  % Torrey eqn 59a, with same assumptions listed above, 
    R1p_fast(j,k,l) = kba + R1p_fast_b_Torrey;
    

        % short time approach.  I'll use the Morris/Torrey notation that
        % normailzes variables by w1.  It doesn't make any difference in
        % A,B,C amplitudes, and it matches my notebook and Mathematica
        % work.  I'll divide the amps by Mb0, to match zb.
        
        kST(j,k,l) = kba/w1;
        alphaST(j,k,l) = r1b/w1;
        betaST(j,k,l) = r2b/w1;
        deltaST(j,k,l) = dwb/w1;
        deltaAST(j,k,l) = dwa/w1;
        thetaST(j,k,l) = alphaST(j,k,l) - betaST(j,k,l);
        etaST(j,k,l) = 3*(1 + deltaST(j,k,l)^2) - thetaST(j,k,l)^2;
        zetaST(j,k,l) = thetaST(j,k,l)*(9 - 2*thetaST(j,k,l)^2 - 18*deltaST(j,k,l)^2);
        epsilonP_ST(j,k,l) =((zetaST(j,k,l) + sqrt(zetaST(j,k,l)^2 + 4*etaST(j,k,l)^3))/54)^(1/3);
        epsilonM_ST(j,k,l) =nthroot((zetaST(j,k,l) - sqrt(zetaST(j,k,l)^2 + 4*etaST(j,k,l)^3))/54,3);
        
        aST(j,k,l) = 1/3 * (alphaST(j,k,l) + kST(j,k,l) + 2*(betaST(j,k,l)+kST(j,k,l))) - (epsilonP_ST(j,k,l) + epsilonM_ST(j,k,l)); % note alpha->alpha+kba, beta->beta+kba, as dictated by differential eqns
        bST(j,k,l) = 1/3 * (alphaST(j,k,l) + kST(j,k,l) + 2*(betaST(j,k,l)+kST(j,k,l))) + 1/2*(epsilonP_ST(j,k,l) + epsilonM_ST(j,k,l)); % note alpha->alpha+kba, beta->beta+kba, as dictated by differential eqns
        sST(j,k,l) = sqrt(3)/2 * (epsilonP_ST(j,k,l) - epsilonM_ST(j,k,l));
        gammaST(j,k,l) = aST(j,k,l)*( (bST(j,k,l) - aST(j,k,l))^2 + sST(j,k,l)^2 );
        DST_offset(j,k,l) = (alphaST(j,k,l)*(1 + deltaAST(j,k,l)^2)*(betaST(j,k,l)^2 + deltaST(j,k,l)^2 + 2*betaST(j,k,l)*kST(j,k,l) + kST(j,k,l)^2) )/ ...
             ((1 + deltaAST(j,k,l)^2)*(betaST(j,k,l) + kST(j,k,l) + (alphaST(j,k,l) + kST(j,k,l))*(deltaST(j,k,l)^2 + (betaST(j,k,l) + kST(j,k,l))^2)));
        DST_za_coef(j,k,l) = deltaAST(j,k,l)*kST(j,k,l)*(deltaST(j,k,l) + deltaST(j,k,l)^2*deltaAST(j,k,l) + deltaAST(j,k,l)*(betaST(j,k,l) + kST(j,k,l))^2)/ ...
            ((1 + deltaAST(j,k,l)^2)*(betaST(j,k,l) + kST(j,k,l) + (alphaST(j,k,l) + kST(j,k,l))*(deltaST(j,k,l)^2 + (betaST(j,k,l) + kST(j,k,l))^2)));
        AST_offset(j,k,l) = -((alphaST(j,k,l)*(aST(j,k,l)^2 + betaST(j,k,l)^2 + deltaST(j,k,l)^2 + 2*betaST(j,k,l)*kST(j,k,l) + kST(j,k,l)^2 - 2*aST(j,k,l)*(betaST(j,k,l) + kST(j,k,l))))/gammaST(j,k,l)); 
        AST_za_coef(j,k,l) =  - ((deltaST(j,k,l)*deltaAST(j,k,l)*kST(j,k,l) + deltaST(j,k,l)^2*deltaAST(j,k,l)^2*kST(j,k,l) + deltaAST(j,k,l)^2*kST(j,k,l)*(-aST(j,k,l) + betaST(j,k,l) + kST(j,k,l))^2))/ ...
            ((1 + deltaAST(j,k,l)^2)*gammaST(j,k,l));
        AST_zb_coef(j,k,l) =  - ((-(aST(j,k,l)*deltaST(j,k,l)^2) - aST(j,k,l)*deltaST(j,k,l)^2*deltaAST(j,k,l)^2 - aST(j,k,l)*(-aST(j,k,l) + betaST(j,k,l) + kST(j,k,l))^2 - aST(j,k,l)*deltaAST(j,k,l)^2*(-aST(j,k,l) + betaST(j,k,l) + kST(j,k,l))^2))/ ...
            ((1 + deltaAST(j,k,l)^2)*gammaST(j,k,l));
        BST_offset(j,k,l) =  -AST_offset(j,k,l)-DST_offset(j,k,l);
        BST_za_coef(j,k,l) = -AST_za_coef(j,k,l)-DST_za_coef(j,k,l);
        BST_zb_coef(j,k,l) = 1 -AST_zb_coef(j,k,l); 
        CST_offset(j,k,l) = alphaST(j,k,l) + aST(j,k,l)*AST_offset(j,k,l) + bST(j,k,l)*BST_offset(j,k,l);
        CST_za_coef(j,k,l) = (kST(j,k,l)*deltaAST(j,k,l)^2/(1 + deltaAST(j,k,l)^2) ) + aST(j,k,l)*AST_za_coef(j,k,l) + bST(j,k,l)*BST_za_coef(j,k,l);
        CST_zb_coef(j,k,l) = -(alphaST(j,k,l) + kST(j,k,l)) + aST(j,k,l)*AST_zb_coef(j,k,l) + bST(j,k,l)*BST_zb_coef(j,k,l);

        % put in Torrey type approx for rates when calculating amps for
        % 'derivedTorrey' case.  aST->R1p_fast, bST->R2p_b, 
        % sST->rabi_rate_b, gammaST->gammaST_Torrey
        aST_Torrey(j,k,l) = R1p_fast(j,k,l)/w1;
        bST_Torrey(j,k,l) = R2p_b(j,k,l)/w1;
        sST_Torrey(j,k,l) = rabi_rate_b(j,k,l)/w1;
        gammaST_Torrey(j,k,l) = aST_Torrey(j,k,l)*( (bST_Torrey(j,k,l) - aST_Torrey(j,k,l))^2 + sST_Torrey(j,k,l)^2 );
        DST_offset_Torrey(j,k,l) = DST_offset(j,k,l);
        DST_za_coef_Torrey(j,k,l) = DST_za_coef(j,k,l);
        AST_offset_Torrey(j,k,l) = -((alphaST(j,k,l)*(aST_Torrey(j,k,l)^2 + betaST(j,k,l)^2 + deltaST(j,k,l)^2 + 2*betaST(j,k,l)*kST(j,k,l) + kST(j,k,l)^2 - 2*aST_Torrey(j,k,l)*(betaST(j,k,l) + kST(j,k,l))))/gammaST_Torrey(j,k,l)); 
        AST_za_coef_Torrey(j,k,l) =  - ((deltaST(j,k,l)*deltaAST(j,k,l)*kST(j,k,l) + deltaST(j,k,l)^2*deltaAST(j,k,l)^2*kST(j,k,l) + deltaAST(j,k,l)^2*kST(j,k,l)*(-aST_Torrey(j,k,l) + betaST(j,k,l) + kST(j,k,l))^2))/ ...
            ((1 + deltaAST(j,k,l)^2)*gammaST_Torrey(j,k,l));
        AST_zb_coef_Torrey(j,k,l) =  - ((-(aST_Torrey(j,k,l)*deltaST(j,k,l)^2) - aST_Torrey(j,k,l)*deltaST(j,k,l)^2*deltaAST(j,k,l)^2 - aST_Torrey(j,k,l)*(-aST_Torrey(j,k,l) + betaST(j,k,l) + kST(j,k,l))^2 - aST_Torrey(j,k,l)*deltaAST(j,k,l)^2*(-aST_Torrey(j,k,l) + betaST(j,k,l) + kST(j,k,l))^2))/ ...
            ((1 + deltaAST(j,k,l)^2)*gammaST_Torrey(j,k,l));
        BST_offset_Torrey(j,k,l) =  -AST_offset_Torrey(j,k,l)-DST_offset_Torrey(j,k,l);
        BST_za_coef_Torrey(j,k,l) = -AST_za_coef_Torrey(j,k,l)-DST_za_coef_Torrey(j,k,l);
        BST_zb_coef_Torrey(j,k,l) = 1 -AST_zb_coef_Torrey(j,k,l); 
        CST_offset_Torrey(j,k,l) = alphaST(j,k,l) + aST_Torrey(j,k,l)*AST_offset_Torrey(j,k,l) + bST_Torrey(j,k,l)*BST_offset_Torrey(j,k,l);
        CST_za_coef_Torrey(j,k,l) = (kST(j,k,l)*deltaAST(j,k,l)^2/(1 + deltaAST(j,k,l)^2) ) + aST_Torrey(j,k,l)*AST_za_coef_Torrey(j,k,l) + bST_Torrey(j,k,l)*BST_za_coef_Torrey(j,k,l);
        CST_zb_coef_Torrey(j,k,l) = -(alphaST(j,k,l) + kST(j,k,l)) + aST_Torrey(j,k,l)*AST_zb_coef_Torrey(j,k,l) + bST_Torrey(j,k,l)*BST_zb_coef_Torrey(j,k,l);

    
    
    
% 'derived' with  aST->aST_Torrey, etc.
       
      % za (just R1p + ss)    

        z_an_ss(1,j,k,l) = (dwa^2/(dwa^2+w1^2))*r1a./R1p(j,k,l); % Used cos2theta, even though likely adiabatic transition.  Correct?
        z_an_comp_red2_amp(1,1,j,k,l)  = z_exact_i(3) - z_an_ss(1,j,k,l); % z_an_comp_red2_amp(a or b, component #, j,k,l)
                            % z_an_ss has only two components, z_exact_i has 6.
                            % Amps for components 2-6 are 0, i.e. just R1p + ss 
                            
        % put in time dependence
        z_an_comp_red2(:,1,1,j,k,l) = exp(-R1p(j,k,l)*t) * z_an_comp_red2_amp(1,1,j,k,l);   
        
         % ss + sum of components
        z_an(:,1,j,k,l) = ones(size(t))*z_an_ss(1,j,k,l) + sum(squeeze(z_an_comp_red2(:,1,:,j,k,l)), 2); 
        
      % zb (R1p + R1p_fast + Rabi_b cos + Rabi_b sin + ss, ignoring Rabi_a cos, and Rabi_a sin)

        % Solve Bloch when za=za_ss by Moritz.  Use R1p for 'derivedTorrey' case 
%   I used the commented out code below for z_an_ss(2,... in previous calculations, but changed to
%   a version that matches the paper's eqn 17 for this code in order to
%   increase clarity.
%         z_an_ss(2,j,k,l) = (dwa^2*(kba*r1a + r1b*R1p(j,k,l))*(dwb^2 + (kba + r2b)^2) + dwa*dwb*kba*r1a*w1^2 + ...
%             r1b*R1p(j,k,l)*(dwb^2 + (kba + r2b)^2)*w1^2)/(R1p(j,k,l)*(dwa^2 + w1^2)* ...
%             (dwb^2*(kba + r1b) + (kba + r2b)*(kba^2 + r1b*r2b + kba*(r1b + r2b) + w1^2)));

        if dwa == 0     % want to avoid 1/dwa calcuation issue, so solve for dwa==0 as special case
            z_an_ss(2,j,k,l) = ...
            ( r1b*(dwb^2 + (kba + r2b)^2) ) / ...
                ( (kba + r1b)*(dwb^2 + (kba + r2b)^2) + (kba + r2b)*w1^2 );
        else
            z_an_ss(2,j,k,l) = ...
            ( r1b*(dwb^2 + (kba + r2b)^2) + kba*(dwb^2 + (kba + r2b)^2 + (dwb*w1^2)/dwa)*z_an_ss(1,j,k,l) ) / ...
                ( (kba + r1b)*(dwb^2 + (kba + r2b)^2) + (kba + r2b)*w1^2 );
        end
        
        % Use Torrey type approx for rates in amp
        % calc for 'derivedTorrey case'.
        z_an_comp_red2_amp(2,2,j,k,l) = z_exact_i(6)*AST_zb_coef_Torrey(j,k,l) + z_exact_i(3)*AST_za_coef_Torrey(j,k,l) + AST_offset_Torrey(j,k,l);  % R1p_fast(j,k,l)
        z_an_comp_red2_amp(2,5,j,k,l) = z_exact_i(6)*BST_zb_coef_Torrey(j,k,l) + z_exact_i(3)*BST_za_coef_Torrey(j,k,l) + BST_offset_Torrey(j,k,l);  % Rabi_b cos
        z_an_comp_red2_amp(2,6,j,k,l) = z_exact_i(6)*CST_zb_coef_Torrey(j,k,l)/sST_Torrey(j,k,l) + z_exact_i(3)*CST_za_coef_Torrey(j,k,l)/sST_Torrey(j,k,l) + CST_offset_Torrey(j,k,l)/sST_Torrey(j,k,l) ; % Rabi b sin
   
        % make ic match, but still ignoring rabi_a cos component (component 3):
        z_an_comp_red2_amp(2,1,j,k,l) = z_exact_i(6) - z_an_comp_red2_amp(2,2,j,k,l) - z_an_comp_red2_amp(2,5,j,k,l)-z_an_ss(2,j,k,l);
  
        % put in time dependence 
        z_an_comp_red2(:,2,1,j,k,l) = exp(-R1p(j,k,l)*t) * z_an_comp_red2_amp(2,1,j,k,l);  
        z_an_comp_red2(:,2,2,j,k,l) = exp(-R1p_fast(j,k,l)*t) * z_an_comp_red2_amp(2,2,j,k,l);
        z_an_comp_red2(:,2,5,j,k,l) = exp(-R2p_b(j,k,l)*t) .* cos( rabi_rate_b(j,k,l)* t ) * ...
            z_an_comp_red2_amp(2,5,j,k,l); % Rabi_b cos component
        z_an_comp_red2(:,2,6,j,k,l) = exp(-R2p_b(j,k,l)*t) .* sin( rabi_rate_b(j,k,l)* t ) * ...
            z_an_comp_red2_amp(2,6,j,k,l); % Rabi_b sin component

                % ss + sum of components
        z_an(:,2,j,k,l) = ones(size(t))*z_an_ss(2,j,k,l) + sum(squeeze(z_an_comp_red2(:,2,:,j,k,l)), 2); 
        
    
% numeric solution  

    % calculate eigenvalues and eigenvectors
    
    A = [ ...   
        -r2a-kab,   dwa,       -imag(w1),  kab,        0,          0; ...
        -dwa,        -r2a-kab,   real(w1),   0,          kab,        0; ...
        imag(w1),   -real(w1),  -r1a-kab,   0,          0,          kab; ...
        kba,        0,          0,          -r2b-kba,   dwb,       -imag(w1); ...
        0,          kba,        0,          -dwb,        -r2b-kba,   real(w1); ...
        0,          0,          kba,        imag(w1),   -real(w1),  -r1b-kba ...
    ];
 
     B = [0 0 r1a 0 0 r1b]'; 


    [V_returned, D_returned] = eig(A);
    lamda_returned = diag(D_returned);
    
   % sort eigenvalues and eigenvectors into order that matches analytic
   % rate guesses
   
    [~, R1p_eig_index] = min(abs(lamda_returned - (-R1p(j,k,l))));
    [~, R1p_fast_eig_index] = min(abs(lamda_returned - (-R1p_fast(j,k,l))));
    [~, R_osc_a_1_eig_index] = min(abs(lamda_returned - (-R2p_a(j,k,l) + 1i*rabi_rate_a(j,k,l)) ));
    [~, R_osc_a_2_eig_index] = min(abs(lamda_returned - (-R2p_a(j,k,l) - 1i*rabi_rate_a(j,k,l)) ));
    [~, R_osc_b_1_eig_index] = min(abs(lamda_returned - (-R2p_b(j,k,l) + 1i*rabi_rate_b(j,k,l)) ));
    [~, R_osc_b_2_eig_index] = min(abs(lamda_returned - (-R2p_b(j,k,l) - 1i*rabi_rate_b(j,k,l)) ));
    
    eig_reorder = [R1p_eig_index R1p_fast_eig_index R_osc_a_1_eig_index R_osc_a_2_eig_index R_osc_b_1_eig_index R_osc_b_2_eig_index]; 
   
    for n=1:length(eig_reorder);
        if sum(eig_reorder(n) == eig_reorder) > 1  % return if repeated values
            j
            k
            l
            return;
        end
    end
    
    V(:,:,j,k,l) = V_returned(:,eig_reorder);
    lamda(:,j,k,l) = lamda_returned(eig_reorder);
    D(:,:,j,k,l) = diag(lamda(:,j,k,l));
 
    % calculate exact numeric solution and components
    
        % dzdt = A*z + B; -> z= expm(A*t)*(z_i-z_ss) + z_ss, where z_ss = -Inv(A)*B;
        % Rewrite expm(A*t) as V*expm(D*t)*inv(V), since then only need to calc eig
        % values once (via eig command above) and (more importantly) can
        % combine components (next section) to compare to analytic components
    
    z_exact_ss(:,j,k,l) = -inv(A)*B;

    invV = inv(squeeze(V(:,:,j,k,l)));
    for m=1:n_times     
        z_exact(m,:,j,k,l) = V(:,:,j,k,l)*diag(exp(diag(D(:,:,j,k,l)*t(m))))*invV*(z_exact_i-z_exact_ss(:,j,k,l)) + z_exact_ss(:,j,k,l);
    end
    z_comp_weightings(:,j,k,l) = invV*(z_exact_i-z_exact_ss(:,j,k,l));
    for m=1:6 % m = component number 
        z_comp_amp(:,m,j,k,l) = transpose(z_comp_weightings(m,j,k,l) * V(:,m,j,k,l));
        z_comp(:,:,m,j,k,l) =  exp(lamda(m,j,k,l)*t) * transpose(z_comp_amp(:,m,j,k,l));
    end

    % calculate numeric 'reduced' components = amplitude * function(t). Sum
    % of reduced components = sum of components = z_exact 
    
        % previous version had z_comp_red, which just added eig vectors to get reduced components.  Problem with
        % this approach:  lose phase info in Rabi components (3 & 4), i.e. I don't know how much sin vs cos.
        % z_comp_red2 breaks down into sin and cos components.  
    
    z_comp_red2(:,:,1:2,j,k,l) = z_comp(:,:,1:2,j,k,l); % R1p, R1p_fast(j,k,l) components
    z_comp_red2_amp(:,1:2,j,k,l) = squeeze(z_comp_red2(1,:,1:2,j,k,l)); % at time = 0, (t(1)=0)
    z_comp_red2_amp(:,3,j,k,l) = 2*real(z_comp_amp(:,3,j,k,l));  % Rabi_a cos component amp
    z_comp_red2_amp(:,4,j,k,l) = -2*imag(z_comp_amp(:,3,j,k,l)); % Rabi_a sin component amp
    z_comp_red2_amp(:,5,j,k,l) = 2*real(z_comp_amp(:,5,j,k,l));  % Rabi_b cos component amp
    z_comp_red2_amp(:,6,j,k,l) = -2*imag(z_comp_amp(:,5,j,k,l)); % Rabi_b sin component amp
    
    z_comp_red2(:,:,3,j,k,l) = exp(real(lamda(3,j,k,l))*t) .* cos( imag(lamda(3,j,k,l)) * t ) * ...
        transpose(z_comp_red2_amp(:,3,j,k,l)); % Rabi_a cos component
    z_comp_red2(:,:,4,j,k,l) = exp(real(lamda(3,j,k,l))*t) .* sin( imag(lamda(3,j,k,l)) * t ) * ...
        transpose(z_comp_red2_amp(:,4,j,k,l)); % Rabi_a sin component
    z_comp_red2(:,:,5,j,k,l) = exp(real(lamda(5,j,k,l))*t) .* cos( imag(lamda(5,j,k,l)) * t ) * ...
        transpose(z_comp_red2_amp(:,5,j,k,l)); % Rabi_b cos component
    z_comp_red2(:,:,6,j,k,l) = exp(real(lamda(5,j,k,l))*t) .* sin( imag(lamda(5,j,k,l)) * t ) * ...
        transpose(z_comp_red2_amp(:,6,j,k,l)); % Rabi_b sin component
   
end
end
end


% FIGURES  ****************************************************************

ps = 4;     % point skip in plot of analytic solns, freq dim
ps_t = 28;   % point skip in plot of analytic solns, time dim
fs =12;     % xlabel, ylabel font size
markers = 7;
mlw = 0.2;
LineWidth = 1.0;
axis_fontsize = 12;
legend_fontsize = 13;
xlabel_fontsize = 13;
sp_edge = 0.08;      % edge size of subplots
sp_sep = 0.02;     % separation of subplots

% ************** paper version of  ***********************
figure(5);
clf;
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    w1 = w1s(k);
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend

p_num=semilogy( ...
    -w1_offsets_ppm, [-real(lamda([1 2 3 5],:,k,l)); imag(lamda([3 5],:,k,l))], ...     % FLIP SIGN OF X OFFSET +3.5ppm amides
    'LineWidth', LineWidth ...
    );

hold on;
ax = gca;
ax.ColorOrderIndex=1;
ax.FontSize = axis_fontsize;
ax.XDir = 'reverse';                                    % flip

p_an1=plot(...
    -w1_offsets_ppm, R1p(:,k,l), '.:', ...              % flip
    -w1_offsets_ppm, R1p_fast(:,k,l), '.:', ...
    -w1_offsets_ppm, R2p_a(:,k,l), '.:', ...
    -w1_offsets_ppm, R2p_b(:,k,l), '.:', ...
    -w1_offsets_ppm, rabi_rate_a(:,k,l), '.:', ...
    -w1_offsets_ppm, rabi_rate_b(:,k,l), '.:', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps:length(w1_offsets_ppm) ...
    );

p_handle = [p_num;p_an1];

%plot range
xlim( -flip([w1_offsets_ppm_start w1_offsets_ppm_end]));    % flip sign and direction
set(gca,'ytick',[1 100 10000]);


% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '$R_{{\rm{1}}\rho }$'; ...
    '${R_{{\rm{1}}\rho {\rm{,}}\;{\rm{fast}}}}$';
    '${R_{{\rm{2}}\rho {\rm{,}}\;{\rm{a}}}}$'; ...
    '${R_{{\rm{2}}\rho {\rm{,}}\;{\rm{b}}}}$'; ...
    '${\omega _{{\rm{eff,}}\;{\rm{a}}}}$'; ... 
    '${\omega _{{\rm{eff,}}\;{\rm{b}}}}$'; ...
    ' '; ... % no labels for corresponding analytic labels
    ' '; ...
    ' '; ...
    ' '; ...
    ' '; ...
    ' ' ...
    };
legend_order = [1 7 2 8 3 9 4 10 5 11 6 12]; % remove analytic labels 

if l==1 && k==1
    ylabel('no exchange','FontSize',fs,'FontWeight','bold')
    title(['B_1 = ' num2str(w1s(1)/267.5) ' \rm{\muT}'],'FontSize',fs)
    legend_h = legend( ...  
        p_handle(legend_order), ...
        legend_text(legend_order), ... 
        'interpreter','latex','fontsize',legend_fontsize ...
        );
    set(legend_h, 'Location','west');
    axis_w_legend = gca;
elseif l==2 && k==1
    ylabel('amide','FontSize',fs,'FontWeight','bold')
elseif l==3 && k==1
    ylabel('amine','FontSize',fs,'FontWeight','bold')
    xlabel('$\Delta {\omega _{\rm{a}}}\;\left( {{\rm{ppm}}} \right)$','FontSize',xlabel_fontsize,'interpreter','latex')
elseif l==1 && k==2
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==3 && k>1
    xlabel('$\Delta {\omega _{\rm{a}}}\;\left( {{\rm{ppm}}} \right)$','FontSize',xlabel_fontsize,'interpreter','latex')
end
 

% remove tick marks
if k==1 && l<3
    set(gca,'xtick',[]);
elseif l==3 && k>1
    set(gca,'ytick',[]);
elseif not(k==1 && l==3)
    set(gca,'xtick',[],'ytick',[]);
end


% set size of subplots
sp_w = (1-2*sp_edge-3*sp_sep)/4;
sp_h = (1-2*sp_edge-2*sp_sep)/3;
set(sp_handle,'position',[sp_edge+(k-1)*(sp_w+sp_sep), sp_edge+(3-l)*(sp_h+sp_sep), sp_w, sp_h])

end
end

% set legend position
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep), sp_w, sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp)
 



% ************* paper version of  ***************************

figure(6);
clf;
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    w1 = w1s(k);
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend

p_num=plot( ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(3, 1,:,k,l)),  ... %FLIP SIGN TO MATCH amides at 3.5ppm
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(3, 2,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(3, 3,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(3, 4,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(3, 5,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(3, 6,:,k,l)),  ...
    -w1_offsets_ppm, z_exact_ss(3,:,k,l),  ...
    'LineWidth',LineWidth ...
    );
hold on;
ax = gca;
ax.ColorOrderIndex=1;
ax.FontSize = axis_fontsize;
ax.XDir = 'reverse';                                % flip


p_an1=plot(...
    -w1_offsets_ppm, squeeze(z_an_comp_red2_amp(1, 1,:,k,l)), '.:', ... % flip
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps:length(w1_offsets_ppm) ...
    );
ax.ColorOrderIndex=7;
p_an2=plot( ...
    -w1_offsets_ppm, z_an_ss(1,:,k,l),'.:',  ...    % flip
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps:length(w1_offsets_ppm) ...
    );
p_zi=plot( -w1_offsets_ppm, squeeze(z_exact(1, 3, :,k,l)), '--k','linewidth',LineWidth); % 1st value is t=0, which should be initial value
            % flip
p_skip = plot(nan,nan,'LineStyle','none');

p_handle = [p_num;p_an1;p_an2;p_zi;p_skip];

%plot range
xlim( -flip([w1_offsets_ppm_start w1_offsets_ppm_end]));    % flip sign and direction
ylim([-0.2 1.0]);

% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '$c_{\rm{a}}^{{R_{{\rm{1}}\rho }}}$'; ...
    '$c_{\rm{a}}^{{R_{{\rm{1}}\rho {\rm{,}}\;{\rm{fast}}}}}$'; ...
    '$c_{\rm{a}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}$'; ...
    '$c_{\rm{a}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}$'; ... 
    '$c_{\rm{a}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}$'; ...
    '$c_{\rm{a}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}$'; ... 
    '$Z_{\rm{a}}^{{\rm{cw,}}\;{\rm{ss}}}$'; ...
    ' '; ... % no labels for corresponding analytic labels
    ' '; ... 
    '${Z_{\rm{a}}^{{\rm{i}}}}$'; ...  %'${Z_{{\rm{ia}}}}$'; ...
    ' ' ... % added white space for line skipping
    };
legend_order = [1 8 2 11 3 11 4 11 5 11 6 11 7 9 10]; % remove analytic labels and add line skips between c values

if l==1 && k==1
    ylabel('no exchange','FontSize',fs,'FontWeight','bold')
    title(['B_1 = ' num2str(w1s(1)/267.5) ' \rm{\muT}'],'FontSize',fs)
    legend_h = legend( ...  
        p_handle(legend_order), ...
        legend_text(legend_order), ... 
        'interpreter','latex','fontsize',legend_fontsize ...
        );
    set(legend_h, 'Location','west');
    axis_w_legend = gca;
elseif l==2 && k==1
    ylabel('amide','FontSize',fs,'FontWeight','bold')
elseif l==3 && k==1
    ylabel('amine','FontSize',fs,'FontWeight','bold')
    xlabel('$\Delta {\omega _{\rm{a}}}\;\left( {{\rm{ppm}}} \right)$','FontSize',xlabel_fontsize,'interpreter','latex')
elseif l==1 && k==2
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==3 && k>1
    xlabel('$\Delta {\omega _{\rm{a}}}\;\left( {{\rm{ppm}}} \right)$','FontSize',xlabel_fontsize,'interpreter','latex')
end
 
% remove tick marks
if k==1 && l<3
    set(gca,'xtick',[]);
elseif l==3 && k>1
    set(gca,'ytick',[]);
elseif not(k==1 && l==3)
    set(gca,'xtick',[],'ytick',[]);
end


% set size of subplots
sp_w = (1-2*sp_edge-3*sp_sep)/4;
sp_h = (1-2*sp_edge-2*sp_sep)/3;
set(sp_handle,'position',[sp_edge+(k-1)*(sp_w+sp_sep), sp_edge+(3-l)*(sp_h+sp_sep), sp_w, sp_h])

end
end

% set legend position
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep), sp_w, sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp)

  
  
  % *********** paper display version of    ***********************
  
figure(7);
clf;
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    w1 = w1s(k);
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend

p_num=plot( ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(6, 1,:,k,l)),  ... %FLIP SIGN TO MATCH amide at 3.5ppm
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(6, 2,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(6, 3,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(6, 4,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(6, 5,:,k,l)),  ...
    -w1_offsets_ppm, squeeze(z_comp_red2_amp(6, 6,:,k,l)),  ...
    -w1_offsets_ppm, z_exact_ss(6,:,k,l),  ...
    'LineWidth',LineWidth ...
    );
hold on;
ax = gca;
ax.ColorOrderIndex=1;
ax.FontSize = axis_fontsize;
ax.XDir = 'reverse';                                        % flip


p_an1=plot(...
    -w1_offsets_ppm, squeeze(z_an_comp_red2_amp(2, 1,:,k,l)), '.:',  ...  % flip
    -w1_offsets_ppm, squeeze(z_an_comp_red2_amp(2, 2,:,k,l)), '.:',  ...  % flip
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps:length(w1_offsets_ppm) ...
    );
ax.ColorOrderIndex=5;
p_an2=plot(...
    -w1_offsets_ppm, squeeze(z_an_comp_red2_amp(2, 5,:,k,l)), '.:',  ...  % flip
    -w1_offsets_ppm, squeeze(z_an_comp_red2_amp(2, 6,:,k,l)), '.:',  ...  % flip
    -w1_offsets_ppm, z_an_ss(2,:,k,l),'.:', 'MarkerIndices', 1:ps:length(w1_offsets_ppm), ...  % flip
    'markers', markers, 'linewidth', mlw ...
    );
p_zi=plot( -w1_offsets_ppm, squeeze(z_exact(1, 6, :,k,l)), '--k','linewidth',LineWidth); % 1st value is t=0, which should be initial value
            % flip
p_skip = plot(nan,nan,'LineStyle','none');
p_handle = [p_num;p_an1;p_an2;p_zi;p_skip];

%plot range
xlim( -flip([w1_offsets_ppm_start w1_offsets_ppm_end]));    % flip sign and direction
ylim([-0.2 1.0]);

% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '$c_{\rm{b}}^{{R_{{\rm{1}}\rho }}}$'; ...
    '$c_{\rm{b}}^{{R_{{\rm{1}}\rho {\rm{,}}\;{\rm{fast}}}}}$'; ...
    '$c_{\rm{b}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}$'; ...
    '$c_{\rm{b}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}$'; ... 
    '$c_{\rm{b}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}$'; ...
    '$c_{\rm{b}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}$'; ... 
    '$Z_{\rm{b}}^{{\rm{cw,}}\;{\rm{ss}}}$'; ...
    ' '; ... % no labels for corresponding analytic labels
    ' '; ... 
    ' '; ...
    ' '; ...
    ' '; ...
    '${Z_{\rm{b}}^{{\rm{i}}}}$'; ...
    ' ' ... % added white space for line skipping
    };
legend_order = [1 8 2 9 3 14 4 14 5 10 6 11 7 12 13]; % remove analytic labels and add line skips between c values

if l==1 && k==1
    ylabel('no exchange','FontSize',fs,'FontWeight','bold')
    title(['B_1 = ' num2str(w1s(1)/267.5) ' \rm{\muT}'],'FontSize',fs)
    legend_h = legend( ...  
        p_handle(legend_order), ...
        legend_text(legend_order), ... 
        'interpreter','latex','fontsize',legend_fontsize ...
        );
    set(legend_h, 'Location','west');
    axis_w_legend = gca;
elseif l==2 && k==1
    ylabel('amide','FontSize',fs,'FontWeight','bold')
elseif l==3 && k==1
    ylabel('amine','FontSize',fs,'FontWeight','bold')
    xlabel('$\Delta {\omega _{\rm{a}}}\;\left( {{\rm{ppm}}} \right)$','FontSize',xlabel_fontsize,'interpreter','latex')
elseif l==1 && k==2
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==3 && k>1
    xlabel('$\Delta {\omega _{\rm{a}}}\;\left( {{\rm{ppm}}} \right)$','FontSize',xlabel_fontsize,'interpreter','latex')
end
 
% remove tick marks
if k==1 && l<3
    set(gca,'xtick',[]);
elseif l==3 && k>1
    set(gca,'ytick',[]);
elseif not(k==1 && l==3)
    set(gca,'xtick',[],'ytick',[]);
end


% set size of subplots
sp_w = (1-2*sp_edge-3*sp_sep)/4;
sp_h = (1-2*sp_edge-2*sp_sep)/3;
set(sp_handle,'position',[sp_edge+(k-1)*(sp_w+sp_sep), sp_edge+(3-l)*(sp_h+sp_sep), sp_w, sp_h])

end
end

% set legend position
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep), sp_w, sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp)
  



% ************* paper version of : Za vs t ***************************

figure(8);
clf;
tms = 1000*t;   % time in ms
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    w1 = w1s(k);
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend

p_num=plot( ...
    tms, squeeze(z_comp_red2(:, 3, 1,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 3, 2,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 3, 3,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 3, 4,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 3, 5,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 3, 6,offset_res_disp(l),k,l)), ...
    tms, z_exact_ss(3,offset_res_disp(l),k,l)*ones(size(t)),  ...
    tms, z_exact(:, 3, offset_res_disp(l),k,l), 'k', ...
    'LineWidth',LineWidth ...
    );
hold on;
ax = gca;
ax.ColorOrderIndex=1;
ax.FontSize = axis_fontsize;

p_an1=plot( ...
    tms, squeeze(z_an_comp_red2(:, 1, 1,offset_res_disp(l),k,l)), '.:', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
ax.ColorOrderIndex=7;
p_an2=plot( ...
    tms, z_an_ss(1,offset_res_disp(l),k,l)*ones(size(t)),'.:',  ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
p_an3=plot(...
    tms, z_an(:, 1, offset_res_disp(l),k,l), '.:k', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
p_skip = plot(nan,nan,'LineStyle','none');

p_handle = [p_num;p_an1;p_an2;p_an3;p_skip];

%plot range
xlim( [tms(1) tms(end)]);
ylim([-0.2 1.0]);

% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '$c_{\rm{a}}^{{R_{{\rm{1}}\rho }}}{e^{ - {R_{{\rm{1}}\rho }}t}}$'; ...
    '$c_{\rm{a}}^{{R_{{\rm{1}}\rho {\rm{,}}\;{\rm{fast}}}}}{e^{ - {R_{{\rm{1}}\rho {\rm{,}}\;{\rm{fast}}}}t}}$'; ...
    '$\begin{array}{l}c_{\rm{a}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{a}}}}t}}\\\; \times \cos \left( {{\omega _{{\rm{eff,}}\;{\rm{a}}}}t} \right)\end{array}$'; ...
    '$\begin{array}{l}c_{\rm{a}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{a}}}}t}}\\\; \times \sin \left( {{\omega _{{\rm{eff,}}\;{\rm{a}}}}t} \right)\end{array}$'; ... 
    '$\begin{array}{l}c_{\rm{a}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{b}}}}t}}\\\; \times \cos \left( {{\omega _{{\rm{eff,}}\;{\rm{b}}}}t} \right)\end{array}$'; ...
    '$\begin{array}{l}c_{\rm{a}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{b}}}}t}}\\\; \times \sin \left( {{\omega _{{\rm{eff,}}\;{\rm{b}}}}t} \right)\end{array}$'; ... 
    '$Z_{\rm{a}}^{{\rm{cw,}}\;{\rm{ss}}}$'; ...
    '${Z_{\rm{a}}}(t)$'; ...
    ' '; ... % no labels for corresponding analytic labels
    ' '; ...
    ' '; ...
    ' ' ... % added white space for line skipping
    };
legend_order = [1 9 2 12 3 12 4 12 5 12 6 12 7 10 8 11]; 

if l==1 && k==1
    ylabel('no exchange','FontSize',fs,'FontWeight','bold')
    title(['B_1 = ' num2str(w1s(1)/267.5) ' \rm{\muT}'],'FontSize',fs)
    legend_h = legend( ...  
        p_handle(legend_order), ...
        legend_text(legend_order), ... 
        'interpreter','latex','fontsize',legend_fontsize ...
        );
    set(legend_h, 'Location','west');
    axis_w_legend = gca;
elseif l==2 && k==1
    ylabel('amide','FontSize',fs,'FontWeight','bold')
elseif l==3 && k==1
    ylabel('amine','FontSize',fs,'FontWeight','bold')
    xlabel('t \rm{(ms)}','FontSize',fs,'FontWeight','bold')
elseif l==1 && k==2
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==3 && k>1
    xlabel('t \rm{(ms)}','FontSize',fs,'FontWeight','bold')
end
 
% remove tick marks
if k==1 && l<3
    set(gca,'xtick',[]);
elseif l==3 && k>1
    set(gca,'ytick',[]);
elseif not(k==1 && l==3)
    set(gca,'xtick',[],'ytick',[]);
end

% set size of subplots
sp_w = (1-2*sp_edge-3*sp_sep)/4;
sp_h = (1-2*sp_edge-2*sp_sep)/3;
set(sp_handle,'position',[sp_edge+(k-1)*(sp_w+sp_sep), sp_edge+(3-l)*(sp_h+sp_sep), sp_w, sp_h])

end
end

% set legend position
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep), sp_w, sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp)
  
  
  
  
% ************* paper version of : Zb vs t ***************************

figure(9);
clf;
tms = 1000*t;   % time in ms
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    w1 = w1s(k);
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend

p_num=plot( ...
    tms, squeeze(z_comp_red2(:, 6, 1,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 6, 2,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 6, 3,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 6, 4,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 6, 5,offset_res_disp(l),k,l)), ...
    tms, squeeze(z_comp_red2(:, 6, 6,offset_res_disp(l),k,l)), ...
    tms, z_exact_ss(6,offset_res_disp(l),k,l)*ones(size(t)),  ...
    tms, z_exact(:, 6, offset_res_disp(l),k,l), 'k', ...
    'LineWidth',LineWidth ...
    );
hold on;
ax = gca;
ax.ColorOrderIndex=1;
ax.FontSize = axis_fontsize;

p_an1=plot( ...
    tms, squeeze(z_an_comp_red2(:, 2, 1,offset_res_disp(l),k,l)), '.:', ...
    tms, squeeze(z_an_comp_red2(:, 2, 2,offset_res_disp(l),k,l)), '.:', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
ax.ColorOrderIndex=5;
p_an2=plot( ...
    tms, squeeze(z_an_comp_red2(:, 2, 5,offset_res_disp(l),k,l)), '.:', ...
    tms, squeeze(z_an_comp_red2(:, 2, 6,offset_res_disp(l),k,l)), '.:', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
ax.ColorOrderIndex=7;
p_an3=plot( ...
    tms, z_an_ss(2,offset_res_disp(l),k,l)*ones(size(t)),'.:',  ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
p_an4=plot(...
    tms, z_an(:, 2, offset_res_disp(l),k,l), '.:k', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
p_skip = plot(nan,nan,'LineStyle','none');

p_handle = [p_num;p_an1;p_an2;p_an3;p_an4;p_skip];

%plot range
xlim( [tms(1) tms(end)]);
ylim([-0.5 1.0]);

% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '$c_{\rm{b}}^{{R_{{\rm{1}}\rho }}}{e^{ - {R_{{\rm{1}}\rho }}t}}$'; ...
    '$c_{\rm{b}}^{{R_{{\rm{1}}\rho {\rm{,}}\;{\rm{fast}}}}}{e^{ - {R_{{\rm{1}}\rho {\rm{,}}\;{\rm{fast}}}}t}}$'; ...
    '$\begin{array}{l}c_{\rm{b}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{a}}}}t}}\\\; \times \cos \left( {{\omega _{{\rm{eff,}}\;{\rm{a}}}}t} \right)\end{array}$'; ...
    '$\begin{array}{l}c_{\rm{b}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{a}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{a}}}}t}}\\\; \times \sin \left( {{\omega _{{\rm{eff,}}\;{\rm{a}}}}t} \right)\end{array}$'; ... 
    '$\begin{array}{l}c_{\rm{b}}^{\cos ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{b}}}}t}}\\\; \times \cos \left( {{\omega _{{\rm{eff,}}\;{\rm{b}}}}t} \right)\end{array}$'; ...
    '$\begin{array}{l}c_{\rm{b}}^{\sin ,\;{\omega _{{\rm{eff,}}\;{\rm{b}}}}}{e^{ - {R_{{\rm{2}}\rho {\rm{,}}\;{\rm{b}}}}t}}\\\; \times \sin \left( {{\omega _{{\rm{eff,}}\;{\rm{b}}}}t} \right)\end{array}$'; ... 
    '$Z_{\rm{b}}^{{\rm{cw,}}\;{\rm{ss}}}$'; ...
    '${Z_{\rm{b}}}(t)$'; ...
    ' '; ... % no labels for corresponding analytic labels
    ' '; ...
    ' '; ...
    ' '; ...
    ' '; ...
    ' '; ...
    ' ' ... % added white space for line skipping
    };
legend_order = [1 9 2 10 3 15 4 15 5 11 6 12 7 13 8 14]; 

if l==1 && k==1
    ylabel('no exchange','FontSize',fs,'FontWeight','bold')
    title(['B_1 = ' num2str(w1s(1)/267.5) ' \rm{\muT}'],'FontSize',fs)
    legend_h = legend( ...  
        p_handle(legend_order), ...
        legend_text(legend_order), ... 
        'interpreter','latex','fontsize',legend_fontsize ...
        );
    set(legend_h, 'Location','west');
    axis_w_legend = gca;
elseif l==2 && k==1
    ylabel('amide','FontSize',fs,'FontWeight','bold')
elseif l==3 && k==1
    ylabel('amine','FontSize',fs,'FontWeight','bold')
    xlabel('t \rm{(ms)}','FontSize',fs,'FontWeight','bold')
elseif l==1 && k==2
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==3 && k>1
    xlabel('t \rm{(ms)}','FontSize',fs,'FontWeight','bold')
end
 
% remove tick marks
if k==1 && l<3
    set(gca,'xtick',[]);
elseif l==3 && k>1
    set(gca,'ytick',[]);
elseif not(k==1 && l==3)
    set(gca,'xtick',[],'ytick',[]);
end


% set size of subplots
sp_w = (1-2*sp_edge-3*sp_sep)/4;
sp_h = (1-2*sp_edge-2*sp_sep)/3;
set(sp_handle,'position',[sp_edge+(k-1)*(sp_w+sp_sep), sp_edge+(3-l)*(sp_h+sp_sep), sp_w, sp_h])

end
end

% set legend position
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep), sp_w, sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp)


  
% ******** paper version: Za, Zai, Zb vs t to justify short time approx****

figure(4);
clf;
tms = 1000*t;   % time in ms
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    w1 = w1s(k);
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend

zai_proj = z_exact_i(3) * ...
    w1_offsets(offset_res_disp(l))^2/(w1^2+w1_offsets(offset_res_disp(l))^2);
    % z component of projection of Ma_initial onto Beff 
    
p_num_a=plot( ...
    tms, z_exact(:, 3, offset_res_disp(l),k,l), 'b', ...
    'LineWidth',LineWidth ...
    );
hold on;
p_an_a=plot( ...
    tms, zai_proj*ones(size(t)), '.:b', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(tms) ...
    );
p_num_b=plot( ...
    tms, z_exact(:, 6, offset_res_disp(l),k,l), 'r', ...
    'LineWidth',LineWidth ...
    );

p_handle = [p_num_a;p_an_a;p_num_b];

%plot range
xlim( [tms(1) tms(end)]);
ylim([-0.5 1.0]);

% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '${Z_{\rm{a}}}(t)$'; ...
    '${Z_{\rm{a}}^{\rm{i}}\;{\rm{proj}}}$'; ... 
    '${Z_{\rm{b}}}(t)$' ...
    };
legend_order = [1 2 3]; 

if l==1 && k==1
    ylabel('no exchange','FontSize',fs,'FontWeight','bold')
    title(['B_1 = ' num2str(w1s(1)/267.5) ' \rm{\muT}'],'FontSize',fs)
    legend_h = legend( ...  
        p_handle(legend_order), ...
        legend_text(legend_order), ... 
        'interpreter','latex','fontsize',legend_fontsize ...
        );
    set(legend_h, 'Location','west');
    axis_w_legend = gca;
elseif l==2 && k==1
    ylabel('amide','FontSize',fs,'FontWeight','bold')
elseif l==3 && k==1
    ylabel('amine','FontSize',fs,'FontWeight','bold')
    xlabel('t \rm{(ms)}','FontSize',fs,'FontWeight','bold')
elseif l==1 && k==2
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_1 = ' num2str(w1s(k)/267.5) ' \rm{\muT}'],'FontSize',fs)
elseif l==3 && k>1
    xlabel('t \rm{(ms)}','FontSize',fs,'FontWeight','bold')
end
 
% remove tick marks
if k==1 && l<3
    set(gca,'xtick',[]);
elseif l==3 && k>1
    set(gca,'ytick',[]);
elseif not(k==1 && l==3)
    set(gca,'xtick',[],'ytick',[]);
end

% set size of subplots
sp_w = (1-2*sp_edge-3*sp_sep)/4;
sp_h = (1-2*sp_edge-2*sp_sep)/3;
set(sp_handle,'position',[sp_edge+(k-1)*(sp_w+sp_sep), sp_edge+(3-l)*(sp_h+sp_sep), sp_w, sp_h])

end
end

% set legend position
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep), sp_w, sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp)
  
  
  