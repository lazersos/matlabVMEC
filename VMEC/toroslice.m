function toroslice(r,zeta,z,s,varargin)
% TOROSLICE(r,zeta,z,s) Plots a toroidal slice for surface s.
% This function plots a toroidal slice for a given surface s:
% TOROSLICE(r,zeta,z,s,varargin)
% Inputs
% r:        Radial position r(s,theta,zeta)
% zeta:     Toroidal Index of plot
% z:        Vertical position z(s,theta,zeta)
% s:        Vector of surfaces to plot
% varargin: Anything else passed to the routine is then passed to the plot
%           plot routines.
%
% Note:  The magnetic axis is always plotted as a '+' Marker.
% Exmaple Usage (data assumed to have at least 10 flux surfaces)
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      r=cfunct(theta,zeta,data.rbc,data.nfp);
%      lambda=sfunct(theta,zeta,data.lbs,data.nfp);
%      z=sfunct(theta,zeta,data.zbs,data.nfp);
%      toroslice(r,zeta,z,[2 5 10]);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.5
ns=size(squeeze(s),2);
if (nargin==4)
    varargin{1}='Color';
    varargin{2}='blue';
end
% Handle plotting a + to denote the magnetic Axis
if (s(1) == 1)
    plot(r(s(1),1,zeta),z(s(1),1,zeta),varargin{:},'Marker','+');
else
    plot(r(s(1),:,zeta),z(s(1),:,zeta),varargin{:});
end
% Handle varargin
if ns > 1
    hold on
    for i=2:ns
        plot(r(s(i),:,zeta),z(s(i),:,zeta),varargin{:});
    end
    hold off
end
axis equal;
end