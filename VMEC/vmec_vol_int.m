function [bx by bz] = vmec_vol_int( data , xs,ys,zs)
%VMEC_VOL_INT Preformed volume integral over 
%   Detailed explanation goes here

% Setup for quick calculation
nu=4*data.mpol+1;
nv=4*data.ntor*data.nfp+1;
nu4=4*nu;
nv4=4*nv;
theta=0:2*pi/nu4:2*pi;
zeta=0:2*pi/nv4:2*pi;
cosz=cos(zeta);
sinz=sin(zeta);
r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
ru=sfunct(theta,zeta,data.rumns,data.xm,data.xn);
rv=sfunct(theta,zeta,data.rvmns,data.xm,data.xn);
z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
zu=cfunct(theta,zeta,data.zumnc,data.xm,data.xn);
zv=cfunct(theta,zeta,data.zvmnc,data.xm,data.xn);
g=cfunct(theta,zeta,data.gmnc,data.xm,data.xn);
ju=cfunct(theta,zeta,data.currumnc,data.xm,data.xn);
jv=cfunct(theta,zeta,data.currvmnc,data.xm,data.xn);
if data.iasym
    r=data.r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
    ru=data.ru+cfunct(theta,zeta,data.rumnc,data.xm,data.xn);
    rv=data.rv+cfunct(theta,zeta,data.rvmnc,data.xm,data.xn);
    z=data.z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
    zu=data.zu+sfunct(theta,zeta,data.zumns,data.xm,data.xn);
    zv=data.zv+sfunct(theta,zeta,data.zvmns,data.xm,data.xn);
    g=data.g+sfunct(theta,zeta,data.gmns,data.xm,data.xn);
    ju=data.ju+sfunct(theta,zeta,data.currumns,data.xm,data.xn);
    jv=data.jv+sfunct(theta,zeta,data.currvmns,data.xm,data.xn);
end
ju = ju./g;
jv = jv./g;
jr = ju.*ru+jv.*rv;
jphi = jv.*r;
jz = ju.*zu+jv.*zv;
jx = 0.0.*jz;
jy = 0.0.*jz;
for v = 1:nv4
    jx(:,:,v) = jr(:,:,v).*cosz(v)-jphi(:,:,v).*sinz(v);
    jy(:,:,v) = jr(:,:,v).*sinz(v)+jphi(:,:,v).*cosz(v);
end
rminor = r - r(1,1,1);
th     = atan2(z,rminor);
th(th < 0.0) = th(th < 0.0) + 2*pi;
dtheta=zeros(data.ns,nu4+1,nv4+1);
dtheta(:,2:nu4,:) = 0.5.*(th(:,3:nu4+1,:)-th(:,1:nu4-1,:));
dtheta(:,1,:) = 0.5.*(atan2(z(:,2,:),rminor(:,2,:)) - atan2(z(:,nu4,:),rminor(:,nu4,:)));
dzeta = 2*pi/(nu4-1);

bx=0; by=0; bz=0;
for s = 1:data.ns
    tx=0; ty=0; tz=0;
    for u = 1:nu4-1
        for v = 1:nv4-1
            dx = xs - r(s,u,v).*cosz(v);
            dy = ys - r(s,u,v).*sinz(v);
            dz = zs - z(s,u,v);
            rho  = sqrt(dx*dx+dy*dy+dz*dz);
            atopx = jy(s,u,v)*dz - jz(s,u,v)*dy;
            atopy = jz(s,u,v)*dx - jx(s,u,v)*dz;
            atopz = jx(s,u,v)*dy - jy(s,u,v)*dx;
            tx = tx + atopx*dtheta(s,u,v)*dzeta/(rho*rho*rho);
            ty = ty + atopy*dtheta(s,u,v)*dzeta/(rho*rho*rho);
            tz = tz + atopz*dtheta(s,u,v)*dzeta/(rho*rho*rho);
        end
    end
    bx = bx+tx*data.vp(s)/(4*pi*pi)/(data.ns-1);
    by = by+ty*data.vp(s)/(4*pi*pi)/(data.ns-1);
    bz = bz+tz*data.vp(s)/(4*pi*pi)/(data.ns-1);
end
bx = 1.0E-7 * bx ;
by = 1.0E-7 * by ;
bz = 1.0E-7 * bz ;





end

