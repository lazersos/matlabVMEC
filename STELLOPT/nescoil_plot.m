function [ output_args ] = nescoil_plot( data )
%NESCOIL_PLOT Summary of this function goes here
%   Detailed explanation goes here


% Initialize
ntheta = data.nu;
nzeta  = data.nv*data.nfp;
ntheta = 360;
nzeta = 360;
theta = 0:2*pi/(ntheta-1):2*pi;
zeta = 0:2*pi/(nzeta-1):2*pi;
nred = data.nv;

% Transform
rp=cfunct(theta,zeta,data.rmnc_plas,data.xm_plas,data.xn_plas);
zp=sfunct(theta,zeta,data.zmns_plas,data.xm_plas,data.xn_plas);
rs=cfunct(theta,zeta,data.rmnc_surf,data.xm_surf,data.xn_surf);
zs=sfunct(theta,zeta,data.zmns_surf,data.xm_surf,data.xn_surf);
ps=sfunct(theta,zeta,data.phimns,data.xm_pot,data.xn_pot);
ju=cfunct(theta,zeta,data.phimns.*data.xm_pot',data.xm_pot,data.xn_pot);
jv=cfunct(theta,zeta,data.phimns.*data.xn_pot',data.xm_pot,data.xn_pot);

% Smooth
for i=1:size(ps,3)
    ps(1,:,i) = smooth(ps(1,:,i),'rlowess');
end
for i=1:size(ps,2)
    ps(1,i,:) = smooth(ps(1,i,:),'rlowess');
end
min_val = min(ps(1,:,1));
max_val = max(ps(1,:,1));
%ps (ps < min_val ) = min_val;
%ps (ps > max_val ) = max_val;

% Plot 3D
subplot(2,1,1);
isotoro(rp(:,:,1:nred),zp(:,:,1:nred),zeta(1:nred),1);
isotoro(rs(:,:,1:nred),zs(:,:,1:nred),zeta(1:nred),1,ps(:,:,1:nred));
camlight left;
title('Plasma and Coil Surface');

% Plot 2D
subplot(2,1,2);
%edge = ps(1,:,1);
%maxval = 2*std(edge);
%ps(ps > maxval) = maxval;
%ps(ps < -maxval) = -maxval;
[C,h]=contour(zeta,theta,squeeze(ps));
pixplot(zeta(1:nred),theta,squeeze(ps(:,:,1:nred))');
xtemp = xlim;
ytemp = ylim;
verts=C;
nels = verts(2,1);
level(1) = verts(1,1);
j = 1;
C = [];
d1 = j+1;
d2 = j+nels;
C{1} = verts(:,d1:d2);
j = d2+1;
n = 2;
hold on;
while j < size(verts,2)
    nels = verts(2,j);
    level(n) = verts(1,j);
    d1 = j+1;
    d2 = j+nels;
    C{n} = verts(:,d1:d2);
    plot(verts(1,d1:d2),verts(2,d1:d2),'k');
    n=n+1;
    j = d2+1;
end
xlim(xtemp);
ylim(ytemp);
title('Current Potential');
xlabel('Toroidal Angle');
ylabel('Poloidal Angle');
set(gca,'XTick',0:2*pi/3/data.nfp:2*pi./data.nfp);
set(gca,'XTickLabelMode','auto');
set(gca,'YTick',min(theta):max(theta)/3:max(theta));
set(gca,'YTickLabelMode','auto');


end

