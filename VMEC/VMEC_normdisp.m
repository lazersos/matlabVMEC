function VMEC_normdisp( vmec_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ntheta = 89;
nzeta  = 89;
theta=0:2*pi/ntheta:2*pi;
zeta=0:2*pi/nzeta:2*pi;

xm_full = vmec_data.xm;
xn_full = vmec_data.xn;
rmnc_full = vmec_data.rmnc;
zmns_full = vmec_data.zmns;
rmnc_axi  = 0.*vmec_data.rmnc;
zmns_axi  = 0.*vmec_data.zmns;
rumns_axi  = 0.*vmec_data.rmnc;
rvmns_axi  = 0.*vmec_data.zmns;
zumnc_axi  = 0.*vmec_data.rmnc;
zvmnc_axi  = 0.*vmec_data.zmns;

rmnc_axi(xn_full == 0,:) = vmec_data.rmnc(xn_full == 0,:);
zmns_axi(xn_full == 0,:) = vmec_data.zmns(xn_full == 0,:);
rumns_axi(xn_full == 0,:) = vmec_data.rumns(xn_full == 0,:);
rvmns_axi(xn_full == 0,:) = vmec_data.rvmns(xn_full == 0,:);
zumnc_axi(xn_full == 0,:) = vmec_data.zumnc(xn_full == 0,:);
zvmnc_axi(xn_full == 0,:) = vmec_data.zvmnc(xn_full == 0,:);

r_full = cfunct(theta,zeta,rmnc_full,xm_full,xn_full);
z_full = sfunct(theta,zeta,zmns_full,xm_full,xn_full);

r_axi = cfunct(theta,zeta,rmnc_axi,xm_full,xn_full);
z_axi = sfunct(theta,zeta,zmns_axi,xm_full,xn_full);

for v=1:nzeta+1
    %a_full(:,:,v)  = r_full(:,:,v)-100.;
    a_full(:,:,v)  = r_full(:,:,v)-r_full(1,1,v);
    %a_axi(:,:,v) = r_axi(:,:,v)-100.;
    a_axi(:,:,v) = r_axi(:,:,v)-r_axi(1,1,v);
end
rho_axi = sqrt(a_axi.*a_axi+z_axi.*z_axi);
rho_full= sqrt(a_full.*a_full+z_full.*z_full);
dr_real=a_full-a_axi;
%plot(a_axi(:,1,1));

plot(rho_full(:,1,45),dr_real(:,1,45)./dr_real(vmec_data.ns,1,45));
end
