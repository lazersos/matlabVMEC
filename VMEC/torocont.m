function torocont(r,z,val,s)
% TOROCONT(r,z,val,s) Contour plots of a toroidal slice.
% This function contour plots a toroidal slice for a given surface
% slice.  Here s is the index of a given toroidal surface zeta(s).
%
% TOROCONT(r,z,val,s)
% Inputs
% r:    Radial matrix
% z:    Vertical matrix
% val:  Matrix of values to plot.
% s:    Index of toroidal cut to plot.
%
% Exmaple Usage
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      r=cfunct(theta,zeta,data.rbc,data.nfp);
%      z=sfunct(theta,zeta,data.zbs,data.nfp);
%      b=cfunct(theta,zeta,data.be,data.nfp);
%      torocont(r,z,b,1);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

% Fix axis (s=1) values by setting to average of s=2 values.
val2=val; % Copy val to val2 to protect val
val2(1,:,s) = val(2,:,s);
%val2(1,:,s)=2*val(2,:,s)-1*val(3,:,s);
h=pcolor(r(:,:,s),z(:,:,s),val2(:,:,s)); % Plot Val2
set(h,'EdgeColor','none');
axis equal
return
end