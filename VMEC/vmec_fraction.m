function [fpassing,ftrapped] = vmec_fraction(vmec_data)
%VMEC_FRACTION Calculates passing and trapped particle fractions
%   The VMEC_FRACTION function calculates the passing and trapped
%   particle function from a data structure as returned by READ_VMEC.
%   The method is based up on the same one use the BOOTSJ code:
%       K.C. Shaing, B.A. Carreras, N. Dominguez, V.E. Lynch, J.S. Tolliver
%           "Bootstrap current control in stellarators", 
%           Phys. Fluids B1, 1663 (1989). 
%           https://aip.scitation.org/doi/10.1063/1.858945
%   Example:
%       vmec_data=read_vmec('wout_test.nc');
%       [fp,ft] = vmec_fraction(vmec_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    1.00

%Defaults
fpassing = zeros(1,vmec_data.ns);

% Setup grid
nth = max(2.*vmec_data.mpol,64);
nzt = max(2.*vmec_data.ntor,32);
nla = 64;
theta = linspace(0,2*pi,nth);
zeta = linspace(0,2*pi,nzt);
lambda = linspace(0,1,nla);


% Transform fields
b = abs(cfunct(theta,zeta,vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq./vmec_data.nfp));
g = -cfunct(theta,zeta,vmec_data.gmnc,vmec_data.xm_nyq,vmec_data.xn_nyq./vmec_data.nfp);
sumg = trapz(theta,trapz(zeta,g,3),2);
b2avg = trapz(theta,trapz(zeta,b.*b.*g,3),2)./sumg;
bmax=max(b,[],[2 3])';
avgbobm2 = b2avg./(bmax.*bmax);
b = b./repmat(bmax',[1 length(theta) length(zeta)]);
b(b>1) = 1;

% Do problem
for k=2:vmec_data.ns
    if bmax(k) <= 0, continue; end
    avgbpov = zeros(1,nla);
    for l=1:nla
        avgbpov(l) = trapz(theta,trapz(zeta,sqrt(abs(1 - lambda(l).*b(k,:,:))).*g(k,:,:),3),2);
    end
    fpassing(k) = 0.75*avgbobm2(k).*trapz(lambda,lambda./avgbpov).*sumg(k);
end
fpassing(1) = 2.*fpassing(2)-fpassing(3);
ftrapped = 1-fpassing;

return;

end

