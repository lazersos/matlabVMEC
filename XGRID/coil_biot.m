function [bx by bz] = coil_biot(coildata,x,y,z,extcur)
%COIL_BIOT Calculates the field at a point in space from a coil structure
%   COIL_BIOT(coildata,x,y,z,extcur) Calculates the field due to a coil
%   system at a point in space.  The algorithm is based upon:
%   J.D. Hanson and S.P. Hirshman,"Compact expressions for the Biot-Savart
%   fields of a filamentary segment." Phys. Plasmas, Vol 9. No. 10, 2002
%   The x, y, and z inputs can now be vectors. They must not be 0.
%
%   Usage:
%       coil_data=read_coils('coils.test');
%       extcur=[1.2e4 1.2e4 1.2e4 -3.5e3 1.1e6 -2.5e4];
%       coil_data=coil_biot_prep(coil_data);
%       [bx by bz]=coil_biot(coildata,5,0,0,extcur);
%
%   See also read_coils, plot_coils, coil_biot_prep.
%
%   Written by:     S. Lazerson (lazerson@pppl.gov)
%   Version:        2.0
%   Date:           12/2/20

bx=0.0;
by=0.0;
bz=0.0;
fac=1.00e-7; % J in A then B in T
if ~isfield(coildata,'vx')
    coildata=coil_biot_prep(coildata);
end
if (max(size(extcur)) ~= max(coildata.vert(5,:)))
    disp(['ERROR: Number of current groups must equal '...
        num2str(max(coildata.vert(5,:))) '!']);
    bx=-1;
    by=-1;
    bz=-1;
    return
end
% x-x'
x1=x-coildata.xw;
y1=y-coildata.yw;
z1=z-coildata.zw;
nels=coildata.nels;
npts=length(x);
%Geometric Quantities
rw=sqrt(x1.*x1+y1.*y1+z1.*z1);
fa=(rw(:,2:nels)+rw(:,1:nels-1))./...
    (rw(:,2:nels).*rw(:,1:nels-1).*(...
    rw(:,2:nels).*rw(:,1:nels-1)+x1(:,2:nels).*x1(:,1:nels-1)...
    +y1(:,2:nels).*y1(:,1:nels-1)+z1(:,2:nels).*z1(:,1:nels-1)));
%Add currents
jtemp=coildata.vert(4,:);
jtemp(jtemp>0) = 1;
for i=1:max(size(extcur))
    jtemp(coildata.vert(5,:)==i)=extcur(i)*jtemp(coildata.vert(5,:)==i);
end
dx = coildata.dx(1:nels-1).*jtemp(1:nels-1);
dy = coildata.dy(1:nels-1).*jtemp(1:nels-1);
dz = coildata.dz(1:nels-1).*jtemp(1:nels-1);
ax = fa*dx';
ay = fa*dy';
az = fa*dz';
%ax=sum(fa(1:nels-1).*coildata.dx(1:nels-1).*jtemp(1:nels-1));
%ay=sum(fa(1:nels-1).*coildata.dy(1:nels-1).*jtemp(1:nels-1));
%az=sum(fa(1:nels-1).*coildata.dz(1:nels-1).*jtemp(1:nels-1));
bx = fac.*(fa*(coildata.vx.*jtemp(1:nels-1))'-y.*az+z.*ay);
by = fac.*(fa*(coildata.vy.*jtemp(1:nels-1))'-z.*ax+x.*az);
bz = fac.*(fa*(coildata.vz.*jtemp(1:nels-1))'-x.*ay+y.*ax);
%bx=fac.*(sum(fa(1:nels-1).*coildata.vx(1:nels-1).*jtemp(1:nels-1))-y.*az+z.*ay);
%by=fac.*(sum(fa(1:nels-1).*coildata.vy(1:nels-1).*jtemp(1:nels-1))-z.*ax+x.*az);
%bz=fac.*(sum(fa(1:nels-1).*coildata.vz(1:nels-1).*jtemp(1:nels-1))-x.*ay+y.*ax);
return
end

