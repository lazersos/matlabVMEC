function [S11 S12 S21 S22] = calc_susceptance(data)
%CALC_SUSCEPTANCE(vmec_data) Calculate susceptance matrix values.
%   The CALC_SUSCEPTANCE routine calculates the susceptance matrix values
%   given a vmec_data structure as returned by READ_VMEC.
%
%   Usage:
%       vmec_data=read_vmec('wout_test.nc');
%       [S11 S12 S21 S22]=calc_susceptance(vmec_data);
%
%   See also read_vmec.
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Version:    1.00
%   Date:       06/28/17

% Assume VMEC for now
nu = 64;
nv = 32;
theta = 0:2*pi/(nu-1):2*pi;
zeta = 0:2*pi/(nv-1):2*pi;

% setup derivatives
for i=1:data.mnmax
    rumns(i,:) =-data.xm(i).*data.rmnc(i,:);
    zumnc(i,:) = data.xm(i).*data.zmns(i,:);
    rvmns(i,:) =-data.xn(i).*data.rmnc(i,:);
    zvmnc(i,:) = data.xn(i).*data.zmns(i,:);
    lumnc(i,:) = data.xm(i).*data.lmns(i,:);
    lvmnc(i,:) = data.xn(i).*data.lmns(i,:);
end
if data.iasym
    for i=1:data.mnmax
        rumnc(i,:) = data.xm(i).*data.rmns(i,:);
        zumns(i,:) =-data.xm(i).*data.zmnc(i,:);
        rvmnc(i,:) = data.xn(i).*data.rmns(i,:);
        zvmns(i,:) =-data.xn(i).*data.zmnc(i,:);
        lumns(i,:) =-data.xm(i).*data.lmnc(i,:);
        lvmns(i,:) =-data.xn(i).*data.lmnc(i,:);
    end
end

% Now Fourier Transform to real space
r  = cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
ru = sfunct(theta,zeta,rumns,data.xm,data.xn);
zu = cfunct(theta,zeta,zumnc,data.xm,data.xn);
rv = sfunct(theta,zeta,rvmns,data.xm,data.xn);
zv = cfunct(theta,zeta,zvmnc,data.xm,data.xn);
lu = cfunct(theta,zeta,lumnc,data.xm,data.xn);
lv = cfunct(theta,zeta,lvmnc,data.xm,data.xn);
g  = cfunct(theta,zeta,data.gmnc,data.xm,data.xn);
if data.iasym
    r  =  r + sfunct(theta,zeta,data.rmns,data.xm,data.xn);
    ru = ru + cfunct(theta,zeta,rumnc,data.xm,data.xn);
    zu = zu + sfunct(theta,zeta,zumns,data.xm,data.xn);
    rv = rv + cfunct(theta,zeta,rvmnc,data.xm,data.xn);
    zv = zv + sfunct(theta,zeta,zvmns,data.xm,data.xn);
    lu = lu + sfunct(theta,zeta,lumns,data.xm,data.xn);
    lv = lv + sfunct(theta,zeta,lvmns,data.xm,data.xn);
     g =  g + sfunct(theta,zeta,data.gmns,data.xm,data.xn);
end

% CALC S11
S11 = (ru.*ru + zu.*zu)./g;
S12 = ((ru.*rv + zu.*zv).*(1+lu) - (ru.*ru+zu.*zu).*lv)./g;
S21 = (ru.*rv + zu.*zv)./g;
S22 = ((rv.*rv + zv.*zv + r.*r).*(1+lu) - (ru.*rv+zu.*zv).*lv)./g;
S11 = trapz(theta,trapz(zeta,S11,3),2)./(4*pi*pi);
S12 = trapz(theta,trapz(zeta,S12,3),2)./(4*pi*pi);
S21 = trapz(theta,trapz(zeta,S21,3),2)./(4*pi*pi);
S22 = trapz(theta,trapz(zeta,S22,3),2)./(4*pi*pi);


end

