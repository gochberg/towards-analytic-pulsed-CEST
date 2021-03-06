function dzdt = bloch_coupled_zaiss(t,z,r1a,r2a,dwa,r1b,r2b,dwb,kab,kba,w1)

%CHANGE 9/13/2016.  Flip signs in front of dwa and dwb so consistent LHR

A = [ ...   
        -r2a-kab,   dwa,       -imag(w1),  kab,        0,          0; ...
        -dwa,        -r2a-kab,   real(w1),   0,          kab,        0; ...
        imag(w1),   -real(w1),  -r1a-kab,   0,          0,          kab; ...
        kba,        0,          0,          -r2b-kba,   dwb,       -imag(w1); ...
        0,          kba,        0,          -dwb,        -r2b-kba,   real(w1); ...
        0,          0,          kba,        imag(w1),   -real(w1),  -r1b-kba ...

    ];

% in bloch_coupled.m, I used the below, but changed -2*pi*dwa to dwa, -2*pi*dwb to dwb, w1 to -w1 to match notation of Moritz's 2013 review eqn 1,etc.
% I also switched kba and kab such that get the same coupling across a row
% A = [ ...
%         -r2a-kab,    2*pi*dwa,   0,          kba,        0,          0; ...
%         -2*pi*dwa,  -r2a-kab,    -w1,        0,          kba,        0; ...
%         0,          w1,         -r1a-kab,    0,          0,          kba; ...
%         kab,        0,          0,          -r2b-kba,    2*pi*dwb,    0; ...
%         0,          kab,        0,          -2*pi*dwb,  -r2b-kba,    -w1; ...
%         0,          0,          kab,        0,          w1,         -r1b-kba ...
% 
%     ];


B = [0;0;r1a;0;0;r1b];

dzdt = A*z + B;
