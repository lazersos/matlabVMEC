function [g11, g12, g22, Bhat, abs_jac, L2, L1, dBdt, kp1] = vmec2gist(vmec_data,nu,nv)
%vmec2gist Calculates various metric quantities relate to GIST/GENE
%   The VMEC2GIST code calculates various metric quantities related to the
%   GIST/GENE interface.  The codes takes a VMEC data strcture as returned
%   by the READ_VMEC routine, and a number of poloidal and toroidal grid
%   points as input.  It then calculates various quantities on the half
%   toroidal flux grid.  For details please see:
%       https://doi.org/10.1063/1.3187907
%   this code mimics the STELLOPT impelmentation stellopt_txport.f90.
%
% Example usage
%    vmec_data=read_vmec('wout_test.nc');
%    [g11,g12,g22,Bhat,abs_jac,L2,L1,dBdt,kp1]=vmec2gist(vmec_data,64,32);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0


% Helpers
nu1 = nu-1;
nv1 = nv-1;
u=0:1.0/nu1:1;
v=0:1.0/nv1:1;
v = v./vmec_data.nfp;
theta = u.*2.*pi;
zeta  = v.*2.*pi;

% Create B derivatives
xm2d = repmat(vmec_data.xm',[1 vmec_data.ns]);
xn2d = repmat(vmec_data.xn',[1 vmec_data.ns]);
xm2dn = repmat(vmec_data.xm_nyq',[1 vmec_data.ns]);
xn2dn = repmat(vmec_data.xn_nyq',[1 vmec_data.ns]);
lumnc = vmec_data.lmns.*xm2d;
lvmnc = vmec_data.lmns.*xn2d;
bumns = -vmec_data.bmnc.*xm2dn;
bvmns = -vmec_data.bmnc.*xn2dn;

if vmec_data.iasym==1
    lumns = -vmec_data.lmnc.*xm2d;
    lvmns = -vmec_data.lmnc.*xn2d;
    bumnc = vmec_data.bmns.*xm2dn;
    bvmnc = vmec_data.bmns.*xn2dn;
end

% Transform
r    = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z    = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
l    = sfunct(theta,zeta,vmec_data.lmns,vmec_data.xm,vmec_data.xn);
lu   = cfunct(theta,zeta,lumnc,vmec_data.xm,vmec_data.xn);
lv   = cfunct(theta,zeta,lvmnc,vmec_data.xm,vmec_data.xn);
ru   = sfunct(theta,zeta,vmec_data.rumns,vmec_data.xm,vmec_data.xn);
rv   = sfunct(theta,zeta,vmec_data.rvmns,vmec_data.xm,vmec_data.xn);
zu   = cfunct(theta,zeta,vmec_data.zumnc,vmec_data.xm,vmec_data.xn);
zv   = cfunct(theta,zeta,vmec_data.zvmnc,vmec_data.xm,vmec_data.xn);
b    = cfunct(theta,zeta,vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
bu   = sfunct(theta,zeta,bumns,vmec_data.xm_nyq,vmec_data.xn_nyq);
bv   = sfunct(theta,zeta,bvmns,vmec_data.xm_nyq,vmec_data.xn_nyq);
if vmec_data.iasym==1
    r    = r  + sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z    = z  + cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
    l    = l  + cfunct(theta,zeta,vmec_data.lmnc,vmec_data.xm,vmec_data.xn);
    lu   = lu + sfunct(theta,zeta,lumns,vmec_data.xm,vmec_data.xn);
    lv   = lv + sfunct(theta,zeta,lvmns,vmec_data.xm,vmec_data.xn);
    ru   = ru + cfunct(theta,zeta,vmec_data.rumnc,vmec_data.xm,vmec_data.xn);
    rv   = rv + cfunct(theta,zeta,vmec_data.rvmnc,vmec_data.xm,vmec_data.xn);
    zu   = zu + sfunct(theta,zeta,vmec_data.zumns,vmec_data.xm,vmec_data.xn);
    zv   = zv + sfunct(theta,zeta,vmec_data.zvmns,vmec_data.xm,vmec_data.xn);
    b    = b  + sfunct(theta,zeta,vmec_data.bmns,vmec_data.xm_nyq,vmec_data.xn_nyq);
    bu   = bu + cfunct(theta,zeta,bumnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
    bv   = bv + cfunct(theta,zeta,bvmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
end

% To half grid
ns1  = vmec_data.ns-1;
ns   = vmec_data.ns;
ds   = 1./double(ns1);
rs = diff(r)./ds;
zs = diff(z)./ds;
bs = diff(b)./ds;
ls = diff(l)./ds;
rh   = (r(1:ns1,:,:)+r(2:ns,:,:)).*0.5;
lh   = (l(1:ns1,:,:)+l(2:ns,:,:)).*0.5;
bh   = (b(1:ns1,:,:)+b(2:ns,:,:)).*0.5;
ru   = (ru(1:ns1,:,:)+ru(2:ns,:,:)).*0.5;
rv   = (rv(1:ns1,:,:)+rv(2:ns,:,:)).*0.5;
zu   = (zu(1:ns1,:,:)+zu(2:ns,:,:)).*0.5;
zv   = (zv(1:ns1,:,:)+zv(2:ns,:,:)).*0.5;
lu   = (lu(1:ns1,:,:)+lu(2:ns,:,:)).*0.5;
lv   = (lv(1:ns1,:,:)+lv(2:ns,:,:)).*0.5;
bu   = (bu(1:ns1,:,:)+bu(2:ns,:,:)).*0.5;
bv   = (bv(1:ns1,:,:)+bv(2:ns,:,:)).*0.5;

% Create qprime
iotap = diff(vmec_data.iotaf)./ds;
iotah = (vmec_data.iotaf(1:ns1)+vmec_data.iotaf(2:ns)).*0.5;
qh    = 1./iotah;
qprim = -iotap.*qh.*qh;
ba = abs(vmec_data.phi(end)./(pi*vmec_data.Aminor.^2));

% Vectors
sqrtg = rh.*(ru.*zs-zu.*rs);
es1 = (-zu.*rh)./sqrtg;
es2 = (zu.*rv-ru.*zv)./sqrtg;
es3 = (ru.*rh)./sqrtg;
eu1 = (zs.*rh)./sqrtg;
eu2 = (zv.*rs-rv.*zs)./sqrtg;
eu3 = (-rs.*rh)./sqrtg;
ev1 = zeros(ns1,nu,nv);
ev2 = (ru.*zs-rs.*zu)./sqrtg;
ev3 = zeros(ns1,nu,nv);

% grad(B)
gradB1 = bs.*es1 + bu.*eu1 + bv.*ev1;
gradB2 = bs.*es2 + bu.*eu2 + bv.*ev2;
gradB3 = bs.*es3 + bu.*eu3 + bv.*ev3;

% Adjust eu by lambda
eu1 = eu1 + ls.*es1 + lu.*eu1 + lv.*ev1;
eu2 = eu2 + ls.*es2 + lu.*eu2 + lv.*ev2;
eu3 = eu3 + ls.*es3 + lu.*eu3 + lv.*ev3;

% Create alpha
thetastar = repmat(theta,[ns1 1 nv])-lh;
grada1 = repmat(qprim',[1 nu nv]) .*thetastar.*es1 + repmat(qh',[1 nu nv]).*eu1 - ev1;
grada2 = repmat(qprim',[1 nu nv]) .*thetastar.*es2 + repmat(qh',[1 nu nv]).*eu2 - ev2;
grada3 = repmat(qprim',[1 nu nv]) .*thetastar.*es3 + repmat(qh',[1 nu nv]).*eu3 - ev3;

% Metric elements
gss = es1.*es1+es2.*es2+es3.*es3;
gsa = es1.*grada1+es2.*grada2+es3.*grada3;
gst = es1.*eu1+es2.*eu2+es3.*eu3;
gaa = grada1.*grada1+grada2.*grada2+grada3.*grada3;
gat = grada1.*eu1+grada2.*eu2+grada3.*eu3;

% Jacobian
wrk1 =  es2.*grada3 - es3.*grada2;
wrk2 =  es3.*grada1 - es1.*grada3;
wrk3 =  es1.*grada2 - es2.*grada1;
jac1 = 1./(wrk1.*eu1+wrk2.*eu2+wrk3.*eu3);

% From this point things get a bit screwy alot of GIST deffinitions
Bhat = bh./ba;

s = 0:ds:1;
s = (s(1:ns1)+s(2:ns)).*0.5;
s = repmat(s',[1 nu nv]);
g11 = (gss.*vmec_data.Aminor.^2)./(4*s);
g12 = (gsa.*vmec_data.Aminor.^2).*repmat(iotah',[1 nu nv]);
g22 = (Bhat.*Bhat+g12.*g12)./g11;
abs_jac = abs(jac1.*2.*repmat(qh',[1 nu nv])./(vmec_data.Aminor.^3));

Es1 = (grada2.*eu3-grada3.*eu2).*jac1;
Es2 = (grada3.*eu1-grada1.*eu3).*jac1;
Es3 = (grada1.*eu2-grada2.*eu1).*jac1;
Ea1 = (eu2.*es3-eu3.*es2).*jac1;
Ea2 = (eu3.*es1-eu1.*es3).*jac1;
Ea3 = (eu1.*es2-eu2.*es1).*jac1;
Et1 = (es2.*grada3-es3.*grada2).*jac1;
Et2 = (es3.*grada1-es1.*grada3).*jac1;
Et3 = (es1.*grada2-es2.*grada1).*jac1;

gradB1 = gradB1./ba;
gradB2 = gradB2./ba;
gradB3 = gradB3./ba;

dBds = gradB1.*Es1+gradB2.*Es2+gradB3.*Es3;
dBda = gradB1.*Ea1+gradB2.*Ea2+gradB3.*Ea3;
dBdt = gradB1.*Et1+gradB2.*Et2+gradB3.*Et3;

c = repmat(iotah'.*iotah'.*vmec_data.Aminor.^4,[1 nu nv]);
L1 = repmat(qh',[1 nu nv]).*(dBda + c.*(gss.*gat-gsa.*gst).*dBdt./(4.*Bhat.*Bhat))./sqrt(s);
L2 = 2.*sqrt(s).*(dBds + c.*(gaa.*gst-gsa.*gat).*dBdt./(4.*Bhat.*Bhat));

pprime = repmat(diff(vmec_data.presf)'./ds,[1 nu nv]);

mu0 = pi.*4E-7;

dpdx = -4.*sqrt(s).*pprime.*mu0./(ba.^2);

kp1 = L2 - dpdx./(2.*Bhat);

end

