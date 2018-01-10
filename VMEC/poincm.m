function poincm_VMEC(cdata,varargin)
%POINCM_VMEC(cdata) Creates a Poincare plot from VMEC and Coil data
%   POINCM_VMEC(vdata,cdata) Creates a Poincare plot from VMEC and coil
%   data.

rmin=3.8;
rmax=4.5;
nr=20;
zstart=0;
phistart=0;
maxiter=2000;
dt=0.01;
datatype='coils';
% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'VMEC'}
                outputtype=varargin{i};
                i=i+1;
                vdata=varargin{i};
            case {'EXTCUR'}
                i=i+1;
                extcur=varargin{i};
        end
    end
end
%Begin Program
%Setup Coil Parameters
nels=max(size(cdata));
coils.x=cdata.verts(1,:);
coils.y=cdata.verts(2,:);
coils.z=cdata.verts(3,:);
coils.dx=(cdata.verts(1,2:nels)-cdata.verts(1,1:nels-1)).*cdata.verts(4,1:nels-1);
coils.dy=(cdata.verts(2,2:nels)-cdata.verts(2,1:nels-1)).*cdata.verts(4,1:nels-1);
coils.dz=(cdata.verts(3,2:nels)-cdata.verts(3,1:nels-1)).*cdata.verts(4,1:nels-1);
coils.vx=cdata.verts(2,1:nels-1).*coils.dz-cdata.verts(3,1:nels-1).*coils.dy;
coils.vy=cdata.verts(3,1:nels-1).*coils.dx-cdata.verts(1,1:nels-1).*coils.dz;
coils.vz=cdata.verts(1,1:nels-1).*coils.dy-cdata.verts(2,1:nels-1).*coils.dx;
%Initialize positions
cx=rmin:(rmax-rmin)/(nr-1):rmax;
cy=zeros(1,nr);
cz=zeros(1,nr);
r0=sqrt(cx.*cx+cy.*cy);
plot(r0,cz,'.r');
%Begin iteration
step=0;
while (step < maxiter)
    for i=1:nr
        [bx,by,bz]=biotcoils(cx(i),cy(i),cz(i),coils);
        b=sqrt(bx*bx+by*by+bz*bz);
        sx=bx/b;
        sy=by/b;
        sz=bz/b;
        cx(i)=cx(i)+sx*dt;
        cy(i)=cy(i)+sy*dt;
        cz(i)=cz(i)+sz*dt;
        r=sqrt(cx(i)*cx(i)+cy(i)*cy(i));
        zeta=atand(cy(i)/cx(i));
        if zeta > zetamax
            zeta=zeta-zetamax;
            cx(i)=r*cosd(zeta);
            cy(i)=r*sind(zeta);
            r=sqrt(cx(i)*cx(i)+cy(i)*cy(i));
            hold on
            plot(r,cz(i),'.');
            hold off
        end 
    end
end
end

function [bx,by,bz]=biotcoils(x,y,z,coils)
% Calculates Magnetic field at a point due to coils
fac=1d-7;
nels=max(size(coils.x));
x1=x-coils.x;
y1=y-coils.y;
z1=z-coils.z;
rw=sqrt(x1.*x1+y1.*y1+z1.*z1);
fa=(rw(2:nels)+rw(1:nels-1))./...
    (rw(2:nels)*rw(1:nels-1)*(...
    rw(2:nels)*rw(1:nels-1)+x1(2:nels)*x1(1:nels-1)+...
    y1(2:nels)*y1(1:nels-1)+z1(2:nels)*z1(1:nels-1)));
ax=sum(fa(1:nels-1)*coils.dx(1:nels-1));
ay=sum(fa(1:nels-1)*coils.dy(1:nels-1));
az=sum(fa(1:nels-1)*coils.dz(1:nels-1));
bx=fac*(sum(fa(1:nels-1)*coils.vx(1:nels-1))-y.*az+z.*ay);
by=fac*(sum(fa(1:nels-1)*coils.vy(1:nels-1))-z.*ax+x.*az);
bz=fac*(sum(fa(1:nels-1)*coils.vz(1:nels-1))-x.*ay+y.*ax);
end
