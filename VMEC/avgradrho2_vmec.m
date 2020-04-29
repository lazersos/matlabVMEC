function avgrho2 = avgradrho2_vmec(vmec_data)
%AVGRADRHO2_VMEC Calculates the <|grad(rho)|^2> used for transport.
%   The AVGGRADRHO2_VMEC function calculates <|grad(rho)|^2> on the VMEC
%   radial grid.  It uses 128 poloidal and 64*nfp toroidal points for the
%   integrals.  The value at s=0 is extrapolated.
%
%   Example:
%       vmec_data=read_vmec('wout_test.nc');
%       avgrho2 = avgradrho2_vmec(vmec_data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

avgrho2=[];

ntheta=128;
nzeta = 64;

lasym=vmec_data.iasym==1;

theta = 0:2*pi/(ntheta-1):2*pi;
zeta  = 0:2*pi/(nzeta.*vmec_data.nfp-1):2*pi;
s     = 0:1./(vmec_data.ns):1;
rho   = sqrt(s);

rus = vmec_data.rmnc;
rvs = vmec_data.rmnc;
zuc = vmec_data.zmns;
zvc = vmec_data.zmns;
for mn=1:vmec_data.mnmax
    rus(mn,:) = -rus(mn,:).*vmec_data.xm(mn);
    zuc(mn,:) =  zuc(mn,:).*vmec_data.xm(mn);
    rvs(mn,:) = -rvs(mn,:).*vmec_data.xn(mn);
    zvc(mn,:) =  zvc(mn,:).*vmec_data.xn(mn);
end
if lasym
    ruc = vmec_data.rmns;
    rvc = vmec_data.rmns;
    zus = vmec_data.zmnc;
    zvs = vmec_data.zmnc;
    for mn=1:vmec_data.mnmax
        ruc(mn,:) =  ruc(mn,:).*vmec_data.xm(mn);
        zus(mn,:) = -zus(mn,:).*vmec_data.xm(mn);
        rvc(mn,:) =  rvc(mn,:).*vmec_data.xn(mn);
        zvs(mn,:) = -zvs(mn,:).*vmec_data.xn(mn);
    end
end

r     = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
sqrtg = cfunct(theta,zeta,vmec_data.gmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
ru    = sfunct(theta,zeta,rus,vmec_data.xm,vmec_data.xn);
rv    = sfunct(theta,zeta,rvs,vmec_data.xm,vmec_data.xn);
zu    = cfunct(theta,zeta,zuc,vmec_data.xm,vmec_data.xn);
zv    = cfunct(theta,zeta,zvc,vmec_data.xm,vmec_data.xn);
if lasym
    r     =     r + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    sqrtg = sqrtg + sfunct(theta,zeta,vmec_data.gmns,vmec_data.xm_nyq,vmec_data.xn_nyq);
    ru    =    ru + cfunct(theta,zeta,ruc,vmec_data.xm,vmec_data.xn);
    rv    =    rv + cfunct(theta,zeta,rvc,vmec_data.xm,vmec_data.xn);
    zu    =    zu + sfunct(theta,zeta,zus,vmec_data.xm,vmec_data.xn);
    zv    =    zv + sfunct(theta,zeta,zvs,vmec_data.xm,vmec_data.xn);
end

% Calc metrics
gsr = -zu.*r;
gsp = zu.*rv - ru.*zv;
gsz = ru.*r;
gs  = (gsr.*gsr+gsp.*gsp+gsz.*gsz)./(sqrtg.*sqrtg); %grad(s)^2
for i = 2: vmec_data.ns
    gs(i,:,:) = gs(i,:,:)./(4*rho(i)*rho(i));
end
vp    = sum(sum(sqrtg,2),3);
avgrho2 = sum(sum(gs.*sqrtg,3),2)./vp;
avgrho2(1) = 2*avgrho2(2)-avgrho2(3);
return
end


