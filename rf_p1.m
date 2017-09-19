function p1 = rf_p1(rf_shape_file)

%  alpha      : flip angle of excitation pulse (degrees)
%  pw         : length of pulse

rf_shape = read_phased_rf_shape(rf_shape_file);
%rf_integral = abs(trapz(rf_shape));
% 
max_rf = max(abs(rf_shape));
%p1 = rf_integral/(max_rf*length(rf_shape));

p1 = mean2(rf_shape)/max_rf;  % could p1 be complex????

% was abs(mean2(rf_shape))/max_rf;  



    %simple sum gives almost same result as integral for pulses with many
    %segments.  But for hard.RF, only the mean2 gives correct result, so
    %I'll use it.




%w1_vector = w1_avg*length(rf_shape)*rf_shape/sum(rf_shape);  version in
%shaped_pulse_soln.m -> no good for phased pulses




