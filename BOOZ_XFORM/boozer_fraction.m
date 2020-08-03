function [fpassing,ftrapped] = boozer_fraction(booz_data)
%BOOZER_FRACTION Calculates passing and trapped particle fractions
%   The BOOZER_FRACTION function calculates the passing and trapped
%   particle function from a data structure as returned by READ_BOOZER.
%   The method is based up on the same one use the BOOTSJ code:
%       K.C. Shaing, B.A. Carreras, N. Dominguez, V.E. Lynch, J.S. Tolliver
%           "Bootstrap current control in stellarators", 
%           Phys. Fluids B1, 1663 (1989). 
%           https://aip.scitation.org/doi/10.1063/1.858945
%   Example:
%       booz_data=read_boozer('boozmn_test.nc');
%       [fp,ft] = boozer_fraction(booz_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    1.00

%Defaults
fpassing = zeros(1,booz_data.ns);
ftrapped = zeros(1,booz_data.ns);


% Setup grid
nth = max(2.*booz_data.mpol,64);
nzt = max(2.*booz_data.ntor,32);
nla = 64;
hth = 2.*pi./(nth-1);
hzt = 2.*pi./(nzt-1);
hla = 1./(nla-1);
theta = 0:hth:2*pi;
zeta  = 0:hzt:2*pi;
lambda = 0:hla:1;


% Transform fields
b = abs(cfunct(theta,zeta,booz_data.bmnc,booz_data.xm,booz_data.xn./booz_data.nfp));
g = 1./(b.*b);
sumg = sum(g,[2 3])';
b2avg = sum(b.*b.*g,[2 3])'./sumg;
bmax=max(b,[],[2 3])';
avgbobm2 = b2avg./(bmax.*bmax);
b = b./repmat(bmax',[1 length(theta) length(zeta)]);
b(b>1) = 1;

% Do problem
for k=2:booz_data.ns
    if bmax(k) <= 0, continue; end
    avgbpov = zeros(1,nla);
    for l=1:nla
        avgbpov(l) = sum(sqrt(abs(1 - lambda(l).*b(k,:,:))).*g(k,:,:),[2 3]);
    end
    fpassing(k) = 0.75*avgbobm2(k).*sum(lambda./avgbpov).*hla.*sumg(k);
    ftrapped(k) = 1-fpassing(k);
end

return;

end

