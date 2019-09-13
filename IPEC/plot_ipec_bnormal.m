function [ output_args ] = plot_ipec_bnormal( ipec_data,n)
%plot_ipec_bnormal(ipec_data,n) Make a B-normal plot
%   Detailed explanation goes here

if ~strcmp(ipec_data.type,'IPEC_XBNORMAL_FUN')
    disp('Wrong data type!');
    return;
end

nzeta = 90;
zeta = 0:2*pi/(nzeta-1):2*pi;
coszeta  = cos(n*zeta);
sinzeta  = sin(n*zeta);
r = ipec_data.r(:,end);
z = ipec_data.z(:,end);
r0 = ipec_data.r(1,1);
z0 = ipec_data.z(1,1);
bnr = ipec_data.bnr(:,end);
bni = ipec_data.bni(:,end);
% We need to resort the data
theta = atan2(z-z0,r-r0);
[theta, dex]= sort(theta);
r = r(dex);
z = z(dex);
bnr = bnr(dex)*1.0E4;
bni = bni(dex)*1.0E4;
% Transform
bn  = bnr*coszeta+bni*sinzeta;
% Adjust amplitude
b_max=max(max(max(bn)));
b_min=min(min(min(bn)));
if (b_max > 100 || b_min < -100)
    bn(bn>100) = 100;
    bn(bn<-100) = -100;
end
% For scaling
%bn(bn>20) = 20;
%bn(bn<-20) = -20;
% Now plot
set(gcf,'Position',[1 1 1024 768],'Color','white');
subplot(1,2,1);
pixplot(zeta,theta,bn');
set(gca,'FontSize',18);
set(gca,'XTick',[0 pi 1.95*pi]);
set(gca,'XTickLabel',{'-180','0','180'});
set(gca,'YTick',[-0.95*pi 0 0.95*pi]);
set(gca,'YTickLabel',{'-180','0','180'});
xlabel('Toroidal Angle');
ylabel('Poloidal Angle');
subplot(1,2,2);
r = repmat(r,1,nzeta);
z = repmat(z,1,nzeta);
n1 = size(r,1); n2 = size(r,2);
r2 = zeros([1,n1,n2]);
z2 = zeros([1,n1,n2]);
b2 = zeros([1,n1,n2]);
r2(1,:,:) = r;
z2(1,:,:) = z;
b2(1,:,:) = bn;
isotoro(r2,z2,zeta,1,b2);
zoom(1.5);
colorbar off;
set(gca,'FontSize',18);
axis off;
title('');
ha=colorbar;
set(ha,'FontSize',18);
ylabel(ha,'B.n [G]');









end

