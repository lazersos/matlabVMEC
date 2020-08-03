function [raxis, zaxis] = beams3d_magaxis(beam_data)
%BEASM3D_MAGAXIS Extracts the magnetic axis
%   The BEASM3D_MAGAXIS function returns R,Z position of the magnetic axis
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [raxis,zaxis] = beams3d_magaxis(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


raxis=[];
zaxis=[];
% Assume we want zeta=0 plane
phidex=1;

raxis=zeros(1,beam_data.nphi);
zaxis=zeros(1,beam_data.nphi);
% Make 2D
for i=1:beam_data.nphi-1
    S2D = squeeze(beam_data.S_ARR(:,i,:));
    smin = min(min(S2D));
    [row,col] = find(S2D==smin);
    raxis(i) = beam_data.raxis(row);
    zaxis(i) = beam_data.zaxis(col); 
end
raxis(end) = raxis(1);
zaxis(end) = zaxis(1);

return;

end

