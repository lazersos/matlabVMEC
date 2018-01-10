function [ax ay az] = coil_biot_vec(coildata,x,y,z,extcur)
%COIL_BIOT_VEC Calculates the vector potential at a point in space from a coil structure
%   COIL_BIOT_VEC(coildata,x,y,z,extcur) Calculates the vector potential due to a coil
%   system at a point in space.  The algorithm is based upon:
%   J.D. Hanson and S.P. Hirshman,"Compact expressions for the Biot-Savart
%   fields of a filamentary segment." Phys. Plasmas, Vol 9. No. 10, 2002
%
%   Usage:
%       coil_data=read_coils('coils.test');
%       extcur=[1.2e4 1.2e4 1.2e4 -3.5e3 1.1e6 -2.5e4];
%       coil_data=coil_biot_prep(coil_data);
%       [ax ay az]=coil_biot_vec(coildata,5,0,0,extcur);
%
%   See also read_coils, plot_coils, coil_biot_prep.
%
%   Written by:     S. Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           12/21/10

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
%Geometric Quantities
rw=sqrt(x1.*x1+y1.*y1+z1.*z1);
fa=(rw(2:nels)+rw(1:nels-1))./...
    (rw(2:nels).*rw(1:nels-1).*(...
    rw(2:nels).*rw(1:nels-1)+x1(2:nels).*x1(1:nels-1)...
    +y1(2:nels).*y1(1:nels-1)+z1(2:nels).*z1(1:nels-1)));
%Add currents
jtemp=coildata.vert(4,:);
for i=1:max(size(extcur))
    jtemp(coildata.vert(5,:)==i)=extcur(i)*jtemp(coildata.vert(5,:)==i);
end
ax=sum(fa(1:nels-1).*coildata.dx(1:nels-1).*jtemp(1:nels-1));
ay=sum(fa(1:nels-1).*coildata.dy(1:nels-1).*jtemp(1:nels-1));
az=sum(fa(1:nels-1).*coildata.dz(1:nels-1).*jtemp(1:nels-1));
return
end