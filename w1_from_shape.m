function [ w1 ] = w1_from_shape( w1_vector, t, t0, pw )
%w1_from_shape gets w1 value at a particular time (or times)

index = ceil(length(w1_vector)*(t-t0)/pw);
index( index==0 ) = 1;
index(index > length(w1_vector)) = length(w1_vector); %corrects numerical precission error

w1 = w1_vector(index);


