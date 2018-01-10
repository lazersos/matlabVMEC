function [ output_args ] = torocut( vmec_data, z0)
%TOROCUT(vmec_data,z) Plots cuts of the flux surfaces at a given elevation.
%   Detailed explanation goes here


phimin = 0.0;
phimax = 2.*pi;
nphi = 360.;
phi=phimin:(phimax-phimin)/(nphi-1):phimax;
theta=[0 pi+pi/10.];
rho=round(2:(vmec_data.ns-2)/9:vmec_data.ns);
r=cfunct(theta,phi,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
if vmec_data.iasym
    r=r+sfunct(theta,phi,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
end
cosphi=cos(phi);
sinphi=sin(phi);
r(:,:,:)=0.0;
% Now we need to find zeros of function
if vmec_data.iasym
    for u=1:nphi
        for s=1:length(rho)
            za = @(x1) sfunct(x1,phi(u),vmec_data.zmns(:,rho(s)),vmec_data.xm,vmec_data.xn);
            zb = @(x2) cfunct(x2,phi(u),vmec_data.zmnc(:,rho(s)),vmec_data.xm,vmec_data.xn);
            %ra = @(x1) cfunct(x1,phi(u),vmec_data.rmnc(:,rho(s)),vmec_data.xm,vmec_data.xn);
            %rb = @(x1) sfunct(x1,phi(u),vmec_data.rmns(:,rho(s)),vmec_data.xm,vmec_data.xn);
            %fr  = @(x3) ra(x3)+rb(x3);
            fz  = @(x3) za(x3)+zb(x3)-z0;
            f2 = @(x4) -fz(x4);
            min=fminbnd(fz,-2*pi,0);
            max=fminbnd(f2,0,2*pi);
            if fz(min) < 0
                theta1=fzero(fz,[min max]);
                theta2=fzero(fz,[max min+2*pi]);
                r_temp=cfunct([theta1 theta2],phi(u),vmec_data.rmnc(:,rho(s)),vmec_data.xm,vmec_data.xn)+...
                    sfunct([theta1 theta2],phi(u),vmec_data.rmns(:,rho(s)),vmec_data.xm,vmec_data.xn);
                r(rho(s),1,u)=r_temp(1);
                r(rho(s),2,u)=r_temp(2);
            else
                r(rho(s),1,u)=0.0;
                r(rho(s),2,u)=0.0;
            end
        end
    end
else
    for u=1:nphi
        for s=1:length(rho)
            f = @(x) sfunct(x,phi(u),vmec_data.zmns(:,rho(s)),vmec_data.xm,vmec_data.xn)-z0;
            theta1=fzero(f,[-pi/2 pi/2]);
            theta2=fzero(f,[pi/2 3*pi/2]);
            r_temp=cfunct([theta1 theta2],phi(u),vmec_data.rmnc(:,rho(s)),vmec_data.xm,vmec_data.xn);
            r(rho(s),1,u)=r_temp(1);
            r(rho(s),2,u)=r_temp(2);
        end
    end
end
    % Now make plots
    plot(squeeze(r(1,1,:)).*cosphi(:),squeeze(r(1,1,:)).*sinphi(:),'r--');
    hold on
    for i=1:length(rho)-1
        plot(squeeze(r(rho(i),1,:)).*cosphi(:),squeeze(r(rho(i),1,:)).*sinphi(:),'r');
        plot(squeeze(r(rho(i),2,:)).*cosphi(:),squeeze(r(rho(i),2,:)).*sinphi(:),'r');
    end
    i=length(rho);
    plot(squeeze(r(rho(i),1,:)).*cosphi(:),squeeze(r(rho(i),1,:)).*sinphi(:),'r','LineWidth',2);
    plot(squeeze(r(rho(i),2,:)).*cosphi(:),squeeze(r(rho(i),2,:)).*sinphi(:),'r','LineWidth',2);

end

