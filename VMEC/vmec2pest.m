function u = vmec2pest(s,theta,zeta,lmns,xm,xn)
%VMEC2PEST Convert from VMEC to PEST poloidal angle
%   This subroutine converts from the VMEC polidal angle to the PEST
%   poloidal angle.  It takes a radial gridpoint (s), vmec polidal angle
%   (theta), VMEC toroidal angle (zeta), the LMNS array, the XM array and
%   XN array as input.
%
% Example usage
%     vmec_data=read_vmec('wout_test.nc');
%     u = vmec2pest(25,pi/16,0.0,vmec_data.lmns,vmec_data.xm,vmec_data.xn);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       0.1

dth = 1.0;
n1  = 0;
u  = theta;
th1 = theta;
while (abs(dth) >= 1.0E-05 && n1< 500)
    lam = sfunct(u,zeta,lmns(:,s),xm,xn);
    dlam = cfunct(u,zeta,lmns(:,s).*xm',xm,xn);
    dth = -(u + lam - th1) / (1+dlam);
    u  = u + 0.5*dth;
    n1 = n1 + 1;
end
return;
end

