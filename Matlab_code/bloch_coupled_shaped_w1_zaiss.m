function dzdt = bloch_coupled_shaped_w1_zaiss(t, z, r1a,r2a,dwa,r1b,r2b,dwb,kab,kba,w1_vector,pulse_start, tp)

w1 = w1_from_shape(w1_vector, t, pulse_start, tp);  

dzdt = bloch_coupled_zaiss(t,z,r1a,r2a,dwa,r1b,r2b,dwb,kab,kba,w1);


