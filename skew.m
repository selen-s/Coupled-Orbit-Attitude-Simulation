%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for SLING Astrodynamics team
% Input: (w(w1,w2,w3))
% Output: 3x3 skew symmetric matrix
% This function takes a 3 entry vector and converts it into a skew
% symmetric matrix. Multiplying this matrix by a vector, v, gives the cross
% product (w)x(v).
%
% Written by: Aidan Hertel
% Last updated: 2/3/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wSkewSymMatrix = skew(w)
    wSkewSymMatrix = [0      -w(3)   w(2);
                     w(3)    0       -w(1);
                     -w(2)   w(1)    0   ];
end