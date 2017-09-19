% figure10.m
%
% produced by D.F. Gochberg for paper:
% D.F. Gochberg, M.D. Does, Z. Zu, C.L. Lankford.  Towards an analytic
% solution for pulsed CEST.  NMR in Biomedicine 2017
%
% You are free to use this code for non-commercial purposes, but please
% cite the above manuscript if you use the code, or parts thereof, to help
% produce a manuscript or presentation figure, as appropriate.  Thanks.
%
% Based on: cest_pulsed_analytic_vs_num5_R2016b.m (Dan's note to self)



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


% specifying hard_pulse_w1s, w1_theta, and Bavgp.

hard_pulse_w1s = [1 3 6]*g; % uT to rad/s
w1_theta = 3*pi;
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
kba_disp = 1;

    % calculated
    n_exch = length(kbas); 

%output
n_times_repetition = 200;    % for all segments, there is the i.c. followed by n_times_pulse points
n_repetitions = 2;
n_times = n_repetitions*n_times_repetition + 1; % +1, so can start at time = 0 and end at end of last pause
for j=1:length(tps)
    ts(:,j)=(0:(tps(j)+tds(j))/n_times_repetition:n_repetitions*(tps(j)+tds(j)))';
end
            
n_offsets = 5; % use number such that ab_offsets and delta_from_b_for_plot_ppm land on points
w1_offsets_ppm_start = -3.5;
w1_offsets_ppm_end = -1.5; 
delta_from_b_for_plot_ppm = -0.5;

    % calculated
    w1_offsets_ppm = w1_offsets_ppm_start:(w1_offsets_ppm_end-w1_offsets_ppm_start)/(n_offsets-1):w1_offsets_ppm_end; 
    w1_offsets= w1_offsets_ppm*267.5*B0;    %ppm to rad/s
    for j=1:length(ab_offsets)
        [dummy, offset_res_disp(j)] = min(abs(w1_offsets + ab_offsets(j))); % use +, since different signs
        [dummy, offset_delta_disp(j)] = min(abs(w1_offsets + ab_offsets(j) + delta_from_b_for_plot_ppm*267.5*B0));
    end 
    offset_disp = offset_res_disp(1);
    



 % predefine variables
 z_an = zeros([n_times 2 n_exch n_powers n_offsets]); % Order:  za_an=Maz/Ma0 zb_an=Mbz/Mb0
 z_num = zeros([n_times 6 n_exch n_powers n_offsets]); % Order:  Max/Ma0 May/Ma0 za_num Mbx/Mb0 Mby/Mb0 zb_num
 w1_num = zeros([n_times 1 n_exch n_powers n_offsets]);  % poorly defined function at pulse-pause interface, since value depends on direction of approach.  I'll make arbitrary assignment.
 w1_an = zeros([n_times 1 n_exch n_powers n_offsets]);
 ss_an_dwb0_dwainf = zeros([n_exch n_powers n_offsets]);

 
 
for l=1:n_exch
    kba = kbas(l);
    kab= kba*pbpa;
    ab_offset = ab_offsets(l);  % change in parallel with kba to match amides and amines
for k=1:n_powers
    w1_avg = w1_avgs(k);
    tp = tps(k);
    td = tds(k);
    t = squeeze(ts(:,k));
for j=1:n_offsets
    w1_offset = w1_offsets(j);
    dwa=w1_offset;
    dwb=w1_offset + ab_offset;
   
 
 
% analytic solutions:

% anayltic solution: za_an
% Replace shaped with hard pulse with same integrated w1^2 and w1, but
% shorter duration.

% Divide shaped pulse of duration tp into three times tp1 (pause), tp2 (pulse), and tp3 (pause).


tp2 = tp * p1^2/p2;
tp1 = (tp - tp2)/2;
tp3 = (tp - tp2)/2;
w1_hard = w1_avg*p2/(p1^2);

% calculate analytic steady state results using ss results and Torrey rate
% approximations


for m=1:n_repetitions
    

    % tp1 pause  
    t0 = (tp+td)*(m-1);
    t_end = t0 + tp1;
    pause_time_seg_ind = logical( (t0 <= t) .* (t <= t_end) ); % using .* as logical &&
    tspan = [t0; t(pause_time_seg_ind); t_end];
    if m==1   
        z_i = [za_initial zb_initial]; 
    end 
    
    z_returned = z_evol_pause3(tspan, z_i, r1a, r1b , kba, kab );        
    z_an(pause_time_seg_ind,:,j,k,l) = z_returned(2:(size(z_returned,1) - 1),:); % ignore first (t0) and last (tp1) results
    z_i = z_returned(size(z_returned,1),:); % last value is initial value of next segment
    w1_an(pause_time_seg_ind,j,k,l) = 0;  % just use 0, since using logical indexing
    
    
    % tp2 pulse
    t0 = (tp+td)*(m-1) + tp1;
    t_end = t0 + tp2;
    pulse_time_seg_ind = logical( (t0 <= t) .* (t <= t_end) ); % using .* as logical &&
    tspan = [t0; t(pulse_time_seg_ind); t_end];    
    
    z_returned = z_evol_pulse_rot_full_3_derivedTorrey( tspan, z_i, r1a, r2a, r1b , r2b, kba, kab, dwa, dwb, w1_hard);
   
    z_an(pulse_time_seg_ind,:,j,k,l) = z_returned(2:(size(z_returned,1) - 1),:); % ignore first (t0) and last (tp1) results
    z_i = z_returned(size(z_returned,1),:); % last value is initial value of next segment
    w1_an(pulse_time_seg_ind,j,k,l) = w1_hard;
    
    
    % tp3 + td pause
    t0 = (tp+td)*(m-1) + tp1 + tp2;
    t_end = t0 + tp3 + td;
    pause_time_seg_ind = logical( (t0 <= t) .* (t <= t_end + .000000001) ); % using .* as logical &&.  Added .000000001 since there seems to be a numerical issue at the last time point. Why?
    tspan = [t0; t(pause_time_seg_ind); t_end];

    z_returned = z_evol_pause3(tspan, z_i, r1a, r1b , kba, kab );
    z_an(pause_time_seg_ind,:,j,k,l) = z_returned(2:(size(z_returned,1) - 1),:); % ignore first (t0) and last (tp1) results
    z_i = z_returned(size(z_returned,1),:); % last value is initial value of next segment
    w1_an(pause_time_seg_ind,j,k,l) = 0;
end


% numeric solution

% clear looping variables (since z goes from 2 to 6 components)
clear z_returned z_i;

% set rf shape and amp variables
rf_shape = read_phased_rf_shape(rf_shape_file);
w1_vector = w1_avg*length(rf_shape)*rf_shape/sum(rf_shape);

for m=1:n_repetitions
    
    % pulse (duration tp)
    
    t0 = (tp+td)*(m-1);
    t_end = t0 + tp;
    pulse_time_seg_ind = logical( (t0 <= t) .* (t <= t_end) ); % using .* as logical &&
    tspan = [t0; t(pulse_time_seg_ind); t_end];
    if m==1                                     % otherwise, z_i set by previous segment
        z_i = [0 0 za_initial 0 0 zb_initial]; 
    end 
      
    z_returned = zeros( [length(tspan) 6] );
    [dummy_t, z_temp ] = ode45(@(t_cont,z_cont) bloch_coupled_shaped_w1_zaiss(t_cont, z_cont, r1a,r2a,dwa,r1b,r2b,dwb,kab,kba,w1_vector, t0, tp), ...
        tspan([1 length(tspan)]), z_i); % t_cont and z_cont are just (continuous) function call names.  No other use.
    z_returned([1 length(tspan)],:) = z_temp([1 size(z_temp,1)],:);   % Break up to 2 different calls to ode45, since ode45 gives error if tspan(1)=tspan(2) or tspan(end)=tspan(end-1), which may sometimes happen.
                % Also, by avoiding this error I still use same structure I
                % used for z_an 
    [dummy_t, z_returned(2:(length(tspan)-1),:) ] = ode45(@(t_cont,z_cont) bloch_coupled_shaped_w1_zaiss(t_cont, z_cont, r1a,r2a,dwa,r1b,r2b,dwb,kab,kba,w1_vector, t0, tp), ...
        tspan( 2:(length(tspan)-1) ), z_i);    
    
    z_num(pulse_time_seg_ind,:,j,k,l) = z_returned(2:(size(z_returned,1) - 1),:); % ignore first (t0) and last (tp1) results
    z_i = z_returned(size(z_returned,1),:); % last value is initial value of next segment
    w1_num(pulse_time_seg_ind,j,k,l) = w1_from_shape(w1_vector, t(pulse_time_seg_ind), t0, tp);
    
    % spoiler is inherent, since only the z-components of the z_i will be
    % used during the pause below
    
    % pause
    
    t0 = (tp+td)*(m-1) + tp;
    t_end = t0 + td;
    pause_time_seg_ind = logical( (t0 <= t) .* (t <= t_end + .000000001) ); % using .* as logical &&.  Added .000000001 since there seems to be a numerical issue at the last time point. Why?
    tspan = [t0; t(pause_time_seg_ind); t_end];

    z_returned = z_evol_pause3(tspan, [z_i(3) z_i(6)], r1a, r1b , kba, kab );
    z_num(pause_time_seg_ind,3,j,k,l) = z_returned(2:(size(z_returned,1) - 1),1); % ignore first (t0) and last (tp1) results
    z_num(pause_time_seg_ind,6,j,k,l) = z_returned(2:(size(z_returned,1) - 1),2);
    z_i = [0 0 z_returned(size(z_returned,1),1) 0 0 z_returned(size(z_returned,1),2)]; % last value is initial value of next segment
    w1_num(pause_time_seg_ind,j,k,l) = 0;
    
end

end
end
end

% plots
%% 


% **************** paper version plot parameters ********************

ps = 4;     % point skip in plot of analytic solns, freq dim
ps_t = 10;   % point skip in plot of analytic solns, time dim
fs =12;     % xlabel, ylabel font size
markers = 7;
mlw = 1;
LineWidth = 1.0;
axis_fontsize = 12;
legend_fontsize = 13;
sp_edge = 0.08;      % edge size of subplots
sp_sep = 0.02;     % separation of subplots


% ************* paper version fig 10 ***************************

figure(10);
clf;
for l=1:n_exch
    kba = kbas(l);      % for title 
    ab_offset = ab_offsets(l);      % change in parallel with kba to match amides, amines, hydroxyls
for k=1:n_powers
    w1_avg = w1_avgs(k);
sp_handle = subplot(n_exch,n_powers+1,(n_powers+1)*(l-1)+k);  % replace n_powers with n_powers+1 to allow for legend

tms = 1000*ts;

p_num=plot( ...
    tms(:,k), z_num(:,3,offset_res_disp(l),k,l), tms(:,k), z_num(:,6,offset_res_disp(l),k,l), ...
    'LineWidth',LineWidth ...
    );
hold on;
ax = gca;
ax.ColorOrderIndex=1;
ax.FontSize = axis_fontsize;

p_an1=plot(...
    tms(:,k), z_an(:,1,offset_res_disp(l),k,l), '.:', ...
    tms(:,k), z_an(:,2,offset_res_disp(l),k,l), '.:', ...
    'markers', markers, 'linewidth', mlw, 'MarkerIndices', 1:ps_t:length(ts) ...
    );


p_handle = [p_num;p_an1];

%plot range
xlim( [tms(1,k) tms(end,k)]);
ylim([-0.75 1.01]);

% legend, xlabels, ylabels, titles
legend_text = ...
    { ...
    '${Z_{\rm{a}}}(t)$'; ...
    '${Z_{\rm{b}}}(t)$'; ...
    ' '; ... % no labels for corresponding analytic labels
    ' ' ... 
    };
legend_order = [1 3 2 4]; % remove analytic labels and add line skips between c values

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
    xlabel(['t \rm{(ms)}'],'FontSize',fs,'FontWeight','bold')
elseif l==1 && k==2
    title(['B_{avgpower} = ' num2str(Bavgp(k)) ' \rm{\muT}'],'FontSize',fs)
elseif l==1 && k==3
    title(['B_{avgpower} = ' num2str(Bavgp(k)) ' \rm{\muT}'],'FontSize',fs)
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
set(legend_h,'position',[sp_edge+(4-1)*(sp_w+sp_sep), sp_edge+(3-2)*(sp_h+sp_sep)+0.25*sp_h, sp_w, 0.5*sp_h]);

% set figure size such that print -depsc will generate .eps that looks like
% figure
  pu = get(gcf,'PaperUnits');
  pp = get(gcf,'PaperPosition');
  set(gcf,'Units',pu,'Position',pp)

  