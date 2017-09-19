function p2 = rf_p2(rf_shape_file)



rf_shape = read_phased_rf_shape(rf_shape_file);


max_rf = max(abs(rf_shape))^2;


p2 = abs(mean2(rf_shape.*conj(rf_shape)))/max_rf;  
    %simple sum gives almost same result as integral for pulses with many
    %segments.  But for hard.RF, only the mean2 gives correct result, so
    %I'll use it.




%w1_vector = w1_avg*length(rf_shape)*rf_shape/sum(rf_shape);  version in
%shaped_pulse_soln.m -> no good for phased pulses




