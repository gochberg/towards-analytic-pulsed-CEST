% figure2.m
%
% produced by D.F. Gochberg for paper:
% D.F. Gochberg, M.D. Does, Z. Zu, C.L. Lankford.  Towards an analytic
% solution for pulsed CEST.  NMR in Biomedicine 2017
%
% You are free to use this code for non-commercial purposes, but please
% cite the above manuscript if you use the code, or parts thereof, to help
% produce a manuscript or presentation figure, as appropriate.  Thanks.
%
% Based on: cest_pulsed_ss_an_eqns_vs_num2_R2016b.m (Dan's note to himself)



clear all;

% cd to whichever directory contains this script, avoiding confusion 
% between called functions with same name but different directories
if(~isdeployed)
  	cd(fileparts(mfilename('fullpath')));
end

%pulse sequence

g = 267.5; % uT to rad/s
rf_shape_file = 'gauss.RF';
p1 = rf_p1(rf_shape_file);
p2 = rf_p2(rf_shape_file);


hard_pulse_w1s = [1 3 6]*g; % uT to rad/s
w1_theta = pi;
Bavgp = [0.6 1.8 3.6]; % uT

w1_avg_disp = 1;

    % calculated
    w1_avgs = hard_pulse_w1s*p1^2/p2;
    tps = w1_theta./w1_avgs;
    tds = tps.*( (w1_avgs./Bavgp).^2 * p2/(g^2 * p1^2) - 1);
    n_powers = length(w1_avgs);


%sample variables
pbpa=0.004;
r1a=1;
r2a=10;
r1b=1;
r2b=100;
za_initial=1;
zb_initial=1;

kbas = [.0001 50 1000];
B0 = 9.4;
ab_offsets = [3.5 3.5 2.0]*(267.5*B0);     %ppm to rad/s.  Change ab_offsets in parallel with kbas
kba_disp = 2;

    % calculated
    n_exch = length(kbas); 

%output
max_n_repetitions = 2000;
end_condition = .00001;  % fraction of equilibrium za and zb;

            
n_offsets = 161; 
w1_offsets_ppm_start = -4.5;
w1_offsets_ppm_end = 0.5; 
delta_from_b_for_plot_ppm = -0.5;

    % calculated
    w1_offsets_ppm = w1_offsets_ppm_start:(w1_offsets_ppm_end-w1_offsets_ppm_start)/(n_offsets-1):w1_offsets_ppm_end; 
    w1_offsets= w1_offsets_ppm*267.5*B0;    %ppm to rad/s
    for j=1:length(ab_offsets)
        [dummy, offset_res_disp(j)] = min(abs(w1_offsets + ab_offsets(j))); % use +, since different signs
        [dummy, offset_delta_disp(j)] = min(abs(w1_offsets + ab_offsets(j) + delta_from_b_for_plot_ppm*267.5*B0));
    end 
    offset_disp = offset_res_disp(1);
    
% set rf shape and amp variables
rf_shape = read_phased_rf_shape(rf_shape_file);



 % predefine variables

 ss_an_dwb0_dwainf = zeros([n_offsets, n_powers, n_exch]);
 z_after_pulse = zeros([max_n_repetitions, 2, n_offsets, n_powers, n_exch]);
    z_after_pulse(1,1,:,:,:) = za_initial * ones([n_offsets, n_powers, n_exch]);
    z_after_pulse(1,2,:,:,:) = zb_initial * ones([n_offsets, n_powers, n_exch]);
 pulses_applied = zeros([ n_offsets, n_powers, n_exch]);
 ss_an_dwb0_dwainf = zeros([ n_offsets, n_powers, n_exch]);
 ss_an_dwainf = zeros([ n_offsets, n_powers, n_exch]);
 ss_an = zeros([ n_offsets, n_powers, n_exch]);
 

 
 
for l=1:n_exch
    kba = kbas(l);
    kab= kba*pbpa;
    ab_offset = ab_offsets(l);  % change in parallel with kba to match amides and amines
for k=1:n_powers
    w1_avg = w1_avgs(k);
    tp = tps(k);
    td = tds(k);
parfor j=1:n_offsets
    w1_offset = w1_offsets(j);
    dwa=w1_offset;
    dwb=w1_offset + ab_offset;
    
 [num2str( (l-1)*n_powers*n_offsets+(k-1)*n_offsets+j ) ' of ' num2str(n_exch*n_powers*n_offsets)]
 
 
% analytic solutions:

% anayltic solution: za_an
% Replace shaped with hard pulse with same integrated w1^2 and w1, but
% shorter duration. Dan's self note: See lab note book 26, p.127

% Divide shaped pulse of duration tp into three times tp1 (pause), tp2 (pulse), and tp3 (pause).


tp2 = tp * p1^2/p2;
tp1 = (tp - tp2)/2;
tp3 = (tp - tp2)/2;
w1_hard = w1_avg*p2/(p1^2);

% calculate analytic steady state results using ss results and Torrey rate approximations 

rabi_rate_b = sqrt( w1_hard^2 + dwb^2);
Reff = r1a + (r2a - r1a) * w1_hard^2 / (w1_hard^2 + dwa^2);
gammaM = 2*sqrt( (kba+r2b)*w1_hard^2/kba + (kba+r2b)^2 );
Rex = kab * ab_offset^2 * w1_hard^2 / ((w1_hard^2 + dwa^2)*(1/4 * gammaM^2 + dwb^2)) + ...
        pbpa * r2b * w1_hard^2 / (1/4 * gammaM^2 + dwb^2) + ...
        kab * w1_hard^2 / (w1_hard^2 + dwa^2) * r2b*(r2b + kba) / (1/4 * gammaM^2 + dwb^2);
R1p = Reff+Rex;
R2p_b_Torrey = r2b - 0.5*(r2b-r1b)/(1+(dwb/w1_hard)^2);           
R2p_b = kba + R2p_b_Torrey;
R1p_fast_b_Torrey = (r2b + r1b*(dwb/w1_hard)^2)/(1+(dwb/w1_hard)^2);  % Torrey eqn 59a, with same assumptions listed above, 
R1p_fast = kba + R1p_fast_b_Torrey;


ProjectionFactor = ((dwa^2*tp2)/(w1_hard^2 + dwa^2) + tp1+tp3+td)/(tp + td);

ss_an_dwb0_dwainf(j,k,l) = ...
    r1a / ...
        (...
        r1a*(1-tp2/(tp+td)) + ...
        R1p*tp2/(tp+td) + ...
        pbpa*(1-exp(-kba*(td+tp1+tp3)))*(1 - exp(-R2p_b*tp2)*cos(rabi_rate_b*tp2) - exp(-R2p_b*tp2)*(kba/w1_hard)*sin(rabi_rate_b*tp2)) / ...  % 10/26/16: CORRECTION. (kba/w1_hard) factor in sin term from analytic_pulsed_CEST4.nb
            ( (1+(kba/w1_hard)^2)*(tp+td)*(1-exp(-kba*(td+tp1+tp3) - R2p_b*tp2)*cos(rabi_rate_b*tp2))) ...
        );

ss_an(j,k,l) = ...
    ProjectionFactor*r1a / ...
        (...
        r1a*(1-tp2/(tp+td)) + ...
        R1p*tp2/(tp+td) + ...
        (pbpa*(1-exp(-kba*(td+tp1+tp3)))*w1_hard* ...
        ( ...
            -(sqrt(1 + dwb^2/w1_hard^2)*w1_hard*(dwa^2*(dwb^2 + w1_hard^2) - (dwb^2*(-1+exp(-R1p_fast*tp2)) - w1_hard^2)*(dwb^2 + kba^2 + w1_hard^2) + ...
            dwa*dwb*(dwb^2*(-1+exp(-R1p_fast*tp2)) - w1_hard^2 + exp(-R1p_fast*tp2)*(kba^2 + w1_hard^2)))) + ...
            exp(-R2p_b*tp2)*sqrt(1 + dwb^2/w1_hard^2)*w1_hard*(dwa*dwb*kba^2 + dwa^2*(dwb^2 + w1_hard^2) + ...
            w1_hard^2*(dwb^2 + kba^2 + w1_hard^2))*cos(rabi_rate_b*tp2) + ...
            dwa*(dwa - dwb)*exp(-R2p_b*tp2)*kba*(dwb^2 + w1_hard^2)*sin(rabi_rate_b*tp2) ...
        )) / ...
         ( ...
            (tp+td)*sqrt(1 + dwb^2/w1_hard^2)*(dwa^2 + w1_hard^2)*(dwb^2 + kba^2 + w1_hard^2)*(dwb^2*(-1 + exp(-kba*(td+tp1+tp3) -R1p_fast*tp2)) - w1_hard^2 + exp(-kba*(td+tp1+tp3) - R2p_b*tp2)*w1_hard^2*cos(rabi_rate_b*tp2)) ...
         ) ...
       );
   
ss_an_dwainf(j,k,l) = ...
    r1a / ...
        (...
        r1a*(1-tp2/(tp+td)) + ...
        R1p*tp2/(tp+td) - ...
        ((-1 + exp(-kba*(td+tp1+tp3)))*pbpa*sqrt(1 + dwb^2/w1_hard^2)*w1_hard^3*(-(sqrt(1 + dwb^2/w1_hard^2)*w1_hard) + exp(-R2p_b*tp2)*sqrt(1 + dwb^2/w1_hard^2)*w1_hard*cos(rabi_rate_b*tp2) + kba*exp(-R2p_b*tp2)*sin(rabi_rate_b*tp2))) / ...
  ((tp+td)*(dwb^2 + kba^2 + w1_hard^2)*(dwb^2*(-1 + exp(-kba*(td+tp1+tp3) -R1p_fast*tp2)) - w1_hard^2 + w1_hard^2*exp(-kba*(td+tp1+tp3) - R2p_b*tp2)*cos(rabi_rate_b*tp2))));
    


% numeric solution

w1_vector = w1_avg*length(rf_shape)*rf_shape/sum(rf_shape);


za_change = 1;  %ensures runs at least once
za_after_pulse_previous = za_initial;  % kludge codeing to satisfy parfor requirement for z_after_pulse to be a 'sliced' variable.
z_i = [0 0 za_initial 0 0 zb_initial];


    
for n_pulses_plus_1=2:max_n_repetitions  % use index n_pulses_plus_1 so that z_after_pulse uses the for iteration variable as an index, as required by parfor
    if abs(za_change) <= end_condition
        break
    end
    
    % pulse (duration tp)
  
t0 = (tp+td)*n_pulses_plus_1;
    t_end = t0 + tp;
    tspan = [t0; t_end];  

      
    
    [dummy_t, z_returned_pulse ] = ode45(@(t_cont,z_cont) bloch_coupled_shaped_w1_zaiss(t_cont, z_cont, r1a,r2a,dwa,r1b,r2b,dwb,kab,kba,w1_vector, t0, tp), ...
        tspan([1 length(tspan)]), z_i); % t_cont and z_cont are just (continuous) function call names.  No other use.


    z_after_pulse(n_pulses_plus_1,:,j,k,l) = [z_returned_pulse(size(z_returned_pulse,1),3), z_returned_pulse(size(z_returned_pulse,1),6)];  % m pulses has index m+1 so first value of z_after_pulse corresponds to zero pulses=i.c.
    za_change = z_returned_pulse(size(z_returned_pulse,1),3) - za_after_pulse_previous;  % can't use the more clear z_after_pulse(m+1,1,j,k,l) or z_after_pulse(m,1,j,k,l), since it will not satisfy the requirements of a sliced variable in parfor
    za_after_pulse_previous = z_returned_pulse(size(z_returned_pulse,1),3);
    pulses_applied(j,k,l) = n_pulses_plus_1 - 1;
    z_i = z_returned_pulse(size(z_returned_pulse,1),:); % last value is initial value of next segment

    
    % spoiler is inherent, since only the z-components of the z_i will be
    % used during the pause below
    
    % pause
    
    t0 = (tp+td)*n_pulses_plus_1 + tp;
    t_end = t0 + td;
    tspan = [t0; t_end];  

    z_returned_pause = z_evol_pause3(tspan, [z_i(3) z_i(6)], r1a, r1b , kba, kab );
    z_i = [0 0 z_returned_pause(size(z_returned_pause,1),1) 0 0 z_returned_pause(size(z_returned_pause,1),2)]; % last value is initial value of next segment
    
end

z_after_pulse_ss(:,j,k,l) = [z_returned_pulse(size(z_returned_pulse,1),3), z_returned_pulse(size(z_returned_pulse,1),6)]; 
        %setting to last returned value.  In order to keep z_after_pulse a sliced variable, I can't use the more sensible form z_after_pulse(pulses_applied(j,k,l)+1,:,j,k,l);



end
end
end

%%
% plot


% **************** paper version plot parameters ********************

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



 % *********** paper display version of fig 2  **
  
figure(2);
clf;
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend




p_num=plot( ...
    -w1_offsets_ppm, squeeze(z_after_pulse_ss(1,:,k,l)), 'k', ... % flip freqs to match standard amide at +3.5ppm
    'LineWidth',LineWidth ...
    );
hold on;
ax = gca;
ax.ColorOrderIndex=1;
ax.FontSize = axis_fontsize;
ax.XDir = 'reverse';                  

p_an1=plot( ...
    -w1_offsets_ppm, squeeze(ss_an_dwainf(:,k,l)),'.:k', ...       
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps:length(w1_offsets_ppm) ...
    );




p_handle = [p_num;p_an1];

%plot range
xlim( -[w1_offsets_ppm_end w1_offsets_ppm_start]);   % flip order and sign
ylim([0 1.0]);


% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '$Z_{\rm{a}}^{{\rm{ss,}}\;{\rm{pul}}}\;{\rm{num}}$'; ...
    '${\left. {Z_{\rm{a}}^{{\rm{ss,}}\;{\rm{pul}}}} \right|_{\Delta {\omega _{\rm{a}}} = \infty }}$'; ...
    };

legend_order = [1 2]; 

if l==1 && k==1
    ylabel('no exchange','FontSize',fs,'FontWeight','bold')
    title(['B_{avgpower} = ' num2str(Bavgp(1)) ' \rm{\muT}'],'FontSize',fs)
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
    title(['B_{avgpower} = ' num2str(Bavgp(k)) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_{avgpower} = ' num2str(Bavgp(k)) ' \rm{\muT}'],'FontSize',fs)
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
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep)+0.25*sp_h, sp_w, 0.5*sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp) 
  

  
  
  

