function [xm_vmec,xn_vmec,rmnc,zmns] = delta2vmec(xm_pg,xn_pg,deltamnc)
%DELTA2VMEC Convert P. Garabedian coordinates to VMEC
%   This subroutine converts the Paul Garabedian representation of toroidal
%   shapes
%   R+iZ = exp(iu)*\sum{\Delta_{mn}*exp(-imu+inv)}
%   to the VMEC RMNC/ZMNS representation.
%   The can handle storage of the R_major term in deltamn(0,0)
%
%   Example
%       [xn, xm] = ndgrid(-1:3,-1:4);
%       deltamn  = [ 0.15  0.00  0.10  0.01  0.00  0.01;...
%                    0.09  1.00  3.50 -0.30  0.09 -0.02;...
%                    0.00  0.01  0.00 -0.30 -0.03  0.02;...
%                    0.00 -0.02 -0.01  0.04  0.06  0.00;...
%                    0.00  0.00  0.00  0.02  0.00 -0.02];
%       xm = reshape(xm,[1 numel(xm)]);
%       xn = reshape(xn,[1 numel(xn)]);
%       deltamn = reshape(deltamn,[1 numel(deltamn)]);
%       [xm,xn,rmnc,zmns]=delta2vmec(xm,xn,deltamn');

ns=size(deltamnc,2);
% Handle axisymmetry
if all(xn_pg==0)
    nfp_pg = 1;
    nmax_pg = 0;
    nmin_pg = 0;
else
    dex  = abs(xn_pg)>0;
    nfp_pg = min(abs(xn_pg(dex)));
    nmax_pg = max(xn_pg./nfp_pg);
    nmin_pg = min(xn_pg./nfp_pg);
end
mmax_pg = max(xm_pg);
mmin_pg = min(xm_pg);

% Define VMEC harmonics
nmin_vmec = nmin_pg;
nmax_vmec = nmax_pg;
mmax_vmec = max([mmax_pg-1 1]);
mnmax_vmec = (nmax_vmec-nmin_vmec+1)*(mmax_vmec+1);
xn_vmec=zeros(1,mnmax_vmec); xm_vmec=zeros(1,mnmax_vmec);
rmnc = zeros(mnmax_vmec,ns);
zmns = zeros(mnmax_vmec,ns);
ik=1;
for m = 0:mmax_vmec
    for n = nmin_vmec:nmax_vmec
        %if (n<0 && m==0); continue; end
        xn_vmec(ik) = n;
        xm_vmec(ik) = m;
        ik = ik + 1;
    end
end

% Handle storing of R00 in delta(0,0)
mn00 = and(xn_pg==0,xm_pg==0);
mn10 = and(xn_pg==0,xm_pg==1);
rnorm = deltamnc(mn00,:)/deltamnc(mn10,:);
r00 = deltamnc(mn00,:);
deltamnc(mn00,:) = 1; % Assume PG delta(0,0) is 

% Calculate RMNC/ZMNS from deltamnc
for mpg = mmin_pg:0
    for npg = nmin_pg:nmax_pg
        mvmec = mpg+1;
        nvmec = -npg;
        dvmec = and(xn_vmec == nvmec,xm_vmec==mvmec);
        dpg = and(xn_pg == npg,xm_pg==mpg);
        rmnc(dvmec,:) = rmnc(dvmec,:) + deltamnc(dpg,:);
        zmns(dvmec,:) = zmns(dvmec,:) + deltamnc(dpg,:);
    end
end
for mpg = 1:mmax_pg
    for npg = nmin_pg:nmax_pg
        mvmec = mpg-1;
        nvmec = npg;
        dvmec = and(xn_vmec == nvmec,xm_vmec==mvmec);
        dpg = and(xn_pg == npg,xm_pg==mpg);
        rmnc(dvmec,:) = rmnc(dvmec,:) + deltamnc(dpg,:);
        zmns(dvmec,:) = zmns(dvmec,:) - deltamnc(dpg,:);
    end
end

% Compress m=0 n>1 modes to positive only (VMEC convention)
dexm = and(xm_vmec==0,xn_vmec<0);
dexp = and(xm_vmec==0,xn_vmec>0);
rmnc(dexp,:) = rmnc(dexp,:) + rmnc(dexm,:);
zmns(dexp,:) = zmns(dexp,:) - rmnc(dexm,:);
rmnc(dexm,:) = [];
zmns(dexm,:) = [];
xm_vmec(dexm) = [];
xn_vmec(dexm) = [];

% Renormalize
mnmax_vmec = size(rmnc,1);
for mn = 1:mnmax_vmec
    rmnc(mn,:) = rmnc(mn,:)*rnorm;
    zmns(mn,:) = zmns(mn,:)*rnorm;
end
mn00 = and(xn_vmec==0,xm_vmec==0);
rmnc(mn00,:) = r00;
zmns(mn00,:) = 0;
xn_vmec = xn_vmec.*nfp_pg;

return;
end

