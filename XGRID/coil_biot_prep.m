function coildata=coil_biot_prep(coildata)
%COIL_BIOT_PREP(coildata) Precalculator for COIL_BIOT.
%   COIL_BIOT_PREP(coildata) Calculates helper variables for the COIL_BIOT
%   routine.  It also normalizes all current values >0 to 1.
%
%   Usage:
%       coil_data=read_coils('coil.test');
%       coil_data=coil_biot_prep(coil_data);
%
%   See also read_coils, plot_coils, coil_biot.
%
%   Written by:     S. Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           12/21/10


% First sort the data by current group
coildata.vert=sortrows(coildata.vert',5)';
nels=max(size(coildata.vert));
% Now remove any existing currents
coildata.vert(4,(abs(coildata.vert(4,:))>0))=1;
coildata.xw=coildata.vert(1,:);
coildata.yw=coildata.vert(2,:);
coildata.zw=coildata.vert(3,:);
coildata.jnorm=coildata.vert(4,:)>0;
coildata.dx=(coildata.xw(2:nels)-coildata.xw(1:nels-1)).*coildata.jnorm(1:nels-1);
coildata.dy=(coildata.yw(2:nels)-coildata.yw(1:nels-1)).*coildata.jnorm(1:nels-1);
coildata.dz=(coildata.zw(2:nels)-coildata.zw(1:nels-1)).*coildata.jnorm(1:nels-1);
coildata.vx=coildata.yw(1:nels-1).*coildata.dz-coildata.zw(1:nels-1).*coildata.dy;
coildata.vy=coildata.zw(1:nels-1).*coildata.dx-coildata.xw(1:nels-1).*coildata.dz;
coildata.vz=coildata.xw(1:nels-1).*coildata.dy-coildata.yw(1:nels-1).*coildata.dx;
coildata.nels=nels;
end