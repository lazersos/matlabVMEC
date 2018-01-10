function varargout= VMEC_jcon( data , theta, zeta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu0=4*pi*1.0E-07;


% Create Poloidal Derivatives
% Create Toroidal Derivatives
% Create Radial Derivatives
rho = data.phi./data.phi(data.ns);
for mn=1:data.mnmax
    dbudumns(mn,:)=-data.xm(mn).*data.bsubumnc(mn,:);
    dbvdumns(mn,:)=-data.xm(mn).*data.bsubvmnc(mn,:);
    dbudvmns(mn,:)=-data.xn(mn).*data.bsubumnc(mn,:);
    dbvdvmns(mn,:)=-data.xn(mn).*data.bsubvmnc(mn,:);
    bu_spl=pchip(rho,data.bsubumnc(mn,:));
    bv_spl=pchip(rho,data.bsubvmnc(mn,:));
    bu_prime_spl=spline_deriv(bu_spl,1);
    bv_prime_spl=spline_deriv(bv_spl,1);
    dbudsmnc(mn,:) = ppval(bu_prime_spl,rho);
    dbvdsmnc(mn,:) = ppval(bv_prime_spl,rho);
    %hold on
    %plot(rho,dbudsmnc(mn,:),'r')
    %plot(rho,data.bsubumnc(mn,:),'ok')
    %plot(0:0.0001:1,ppval(bu_spl,0:0.0001:1),'r')
    %hold off
end

% Now Fourier Tranfrom to real space
%dbudu=sfunct(theta,zeta,dbudumns,data.xm,data.xn);
dbvdu=sfunct(theta,zeta,dbvdumns,data.xm,data.xn);
dbudv=sfunct(theta,zeta,dbudvmns,data.xm,data.xn);
%dbvdv=sfunct(theta,zeta,dbvdvmns,data.xm,data.xn);
%dbudu=sfunct(theta,zeta,dbudumns,data.xm,data.xn);
g=cfunct(theta,zeta,data.gmnc,data.xm,data.xn);
curru=cfunct(theta,zeta,data.currumnc,data.xm,data.xn);
currv=cfunct(theta,zeta,data.currvmnc,data.xm,data.xn);
%dbuds=cfunct(theta,zeta,dbudsmnc,data.xm,data.xn);
%dbvds=cfunct(theta,zeta,dbvdsmnc,data.xm,data.xn);

% Now create Current densities
jsups=(dbvdu-dbudv)./(g.*mu0);
jsupu=curru./(g.*mu0);
jsupv=currv./(g.*mu0);

varargout{1}=jsups;
varargout{2}=jsupu;
varargout{3}=jsupv;

%jsupu=(dbsdv-dbvds)./g;
%jsupv=(dbuds-dbsdu)./g;






end

