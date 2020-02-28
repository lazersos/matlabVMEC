function output1 = AL31(X,ZEFF)
%AL31(X,Z) Calculates the L31 Transport Coefficient
%   This function calculates the L31 transport coefficient according to
%       Hirshman, S. Phys. Fluids 31, 3150 (1988)
%   The code takes two inputs and returns the coefficient
%       X:      Ration of trapped to passing particles
%       ZEFF:   Effective Ion Charge Number
%
% Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
% Version:       1.00

Z2 = ZEFF.*ZEFF;
d  = 1.414.*ZEFF + z2 + X.*(0.754+2.657.*ZEFF+2.*Z2) + X.*X.*(0.348+1.243.*ZEFF+Z2);
a  = 0.754 + 2.21.*ZEFF + Z2 + X.*(0.348+1.243.*ZEFF+Z2);
output1 = x.*a./d;
end

