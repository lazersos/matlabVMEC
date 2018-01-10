function plot_displacement( vmec_data,factor )
%PLOT_DISPLACEMENT(vmec_data,factor) Make displacement plot from VMEC
% This funciton creats a 3 pannel plot of the VMEC edge displacement of a
% nominally axisymmetric equilibria.  We assume the n=0 components of
% equilibrium define the 'axisymmetric' state and n/=0 components are the
% perturbation.  Plots are then constructed.  Factor is a scaling factor
% applied to exagerate the displacements in the cross section plot only.
%
% Example usage
%      vmec_data=read_vmec('wout.test');
%      plot_displacement(vmec_data,10);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

theta=0:2*pi./359:2*pi;
zeta =0:2*pi./359:2*pi;

rmnc_axi = vmec_data.rmnc;
zmns_axi = vmec_data.zmns;
rmnc_axi(vmec_data.xn ~= 0,:) = 0.0*rmnc_axi(vmec_data.xn ~= 0,:);
zmns_axi(vmec_data.xn ~= 0,:) = 0.0*zmns_axi(vmec_data.xn ~= 0,:);
rmnc_3d = vmec_data.rmnc;
zmns_3d = vmec_data.zmns;
rmnc_3d(vmec_data.xn ~= 0,:) = factor*rmnc_3d(vmec_data.xn ~= 0,:);
zmns_3d(vmec_data.xn ~= 0,:) = factor*zmns_3d(vmec_data.xn ~= 0,:);

r=cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z=sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
raxi=cfunct(theta,zeta,rmnc_axi,vmec_data.xm,vmec_data.xn);
zaxi=sfunct(theta,zeta,zmns_axi,vmec_data.xm,vmec_data.xn);
r3d=cfunct(theta,zeta,rmnc_3d,vmec_data.xm,vmec_data.xn);
z3d=sfunct(theta,zeta,zmns_3d,vmec_data.xm,vmec_data.xn);
if (vmec_data.iasym)
    rmns_axi = vmec_data.rmns;
    zmnc_axi = vmec_data.zmnc;
    rmns_axi(vmec_data.xn ~= 0,:) = 0.0*rmns_axi(vmec_data.xn ~= 0,:);
    zmnc_axi(vmec_data.xn ~= 0,:) = 0.0*zmnc_axi(vmec_data.xn ~= 0,:);
    rmns_3d = vmec_data.rmns;
    zmnc_3d = vmec_data.zmnc;
    rmns_3d(vmec_data.xn ~= 0,:) = factor*rmns_3d(vmec_data.xn ~= 0,:);
    zmnc_3d(vmec_data.xn ~= 0,:) = factor*zmnc_3d(vmec_data.xn ~= 0,:);
    r=r+sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z=z+cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
    raxi=raxi+sfunct(theta,zeta,rmns_axi,vmec_data.xm,vmec_data.xn);
    zaxi=zaxi+cfunct(theta,zeta,zmnc_axi,vmec_data.xm,vmec_data.xn);
    r3d=r3d+sfunct(theta,zeta,rmns_3d,vmec_data.xm,vmec_data.xn);
    z3d=z3d+cfunct(theta,zeta,zmnc_3d,vmec_data.xm,vmec_data.xn);
end
fig=figure('Position',[ 1 1 1024 786],'Color','white');
subplot(2,2,[1 2])
set(gca,'FontSize',18,'FontWeight','bold');
ns=vmec_data.ns;
%ns = 75;
a1 = raxi-raxi(1,1,1);
a2 = zaxi-zaxi(1,1,1);
theta0 = atan2(a2(ns,:,1),a1(ns,:,1));
rho0 = sqrt(a1.*a1+a2.*a2);
a1 = r - raxi(1,1,1);
a2 = z - zaxi(1,1,1);
for v = 1: length(zeta)
    theta2 = atan2(a2(ns,:,v),a1(ns,:,v));
    x = a1(ns,:,v);
    y = a2(ns,:,v);
    r_temp = pchip(theta2,x,theta0);
    z_temp = pchip(theta2,y,theta0);
    rho2(:,v) = sqrt(r_temp.*r_temp+z_temp.*z_temp);
end
rho_act = rho2-squeeze(rho0(ns,:,:));
pixplot(zeta,theta,squeeze(rho_act)'.*100);
xlabel('Phi');
ylabel('Theta');
%title(['Max displacement: ' num2str(max(max(rho_act)).*100,'%5.3f') ' [cm]; q = ' num2str(1./vmec_data.iotaf(ns),'%5.3f') ]);
title(['Max displacement: ' num2str(max(max(rho_act(1:180,:))).*100,'%5.3f') ' [cm]; q = ' num2str(1./vmec_data.iotaf(ns),'%5.3f') ]);
ha = colorbar;
ylabel(ha,'[cm]','FontSize',18);
set(ha,'FontSize',18,'FontWeight','bold');
subplot(2,2,4);
set(gca,'FontSize',18,'FontWeight','bold');
temp(vmec_data.ns,:,:) = squeeze(rho_act).*100;
isotoro(r3d,z3d,zeta,ns,temp);
ha = colorbar;
set(ha,'FontSize',18,'FontWeight','bold');
title('Exaggerated Surface');
subplot(2,2,3);
set(gca,'FontSize',18,'FontWeight','bold');
nsurf=vmec_data.ns;
plot(r(nsurf,:,1),z(nsurf,:,1),'k');
hold on;
plot(r3d(nsurf,:,1),z3d(nsurf,:,1),'r');
title(['Perturbed Surface (phi=0, ' num2str(factor,'%3d') 'X)']);
xlabel('R [m]');
ylabel('Z [m]');
hold off;
axis equal

end

