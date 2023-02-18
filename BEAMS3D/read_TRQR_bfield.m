function [xaxis, yaxis, zaxis, BX,BY,BZ] = read_TRQR_bfield(x0,fbx,fby,fbz)
%READ_TRQR_BFIELD Reads the B-field for the TRQR code
% This funciton reads the magnetic field inputs files (asde.au1-9) for the
% TRQR code.  It returns the components of the magnetic field and grid
% axes.  The code takes as input the X0 offset of the grid and three
% factors (fbx,fby,fbz) which scale the field (usually 1).  Thes parameter
% can be found at line 73 of the file trq_asd.dat.  
%
% Example usage
%      [xaxis, yaxis, zaxis, BX,BY,BZ] = read_TRQR_bfield(2591,1,1,1);
%
% Note about coordiantes (W7-X Profi)
%   xaxis is along beamline axis toward machine x=0 is magnet entrance
%   yaxis is in the perpendicular direction
%   zaxis is traditional z direction
%   To get profi coodinates use x0 = 0.  Note models are often over 1
%   quadrant.
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

file_names=dir('asde.au*');
% Do first file separately
i=1;
file_name = ['asde.au' num2str(i,'%i')];
fid=fopen(file_name,'r');
temp=fscanf(fid,'%f %f %f',3);
x = temp(1); y = temp(2); z = temp(3);
temp=fscanf(fid,'%i %i %i',3);
nx = temp(1); ny = temp(2); nz = temp(3);
temp=fscanf(fid,'%f %f %f',3);
dx = temp(1); dy = temp(2); dz = temp(3);
temp = fscanf(fid,'%f %f %f %f %f %f %f %f %f',[9 inf]);
fclose(fid);
BX=temp(1,:);
BY=temp(2,:);
BZ=temp(3,:);
% Now read the rest
for i = 2:length(file_names)
    file_name = ['asde.au' num2str(i,'%i')];
    fid=fopen(file_name,'r');
    temp = fscanf(fid,'%f %f %f %f %f %f %f %f %f',[9 inf]);
    fclose(fid);
    BX = [BX temp(1,:)];
    BY = [BY temp(2,:)];
    BZ = [BZ temp(3,:)];
end
% The files are linearly sorted but in ny,nx,nz order
BX = permute(reshape(BX,[nz,nx,ny]).*fbx,[2 3 1]);
BY = permute(reshape(BY,[nz,nx,ny]).*fby,[2 3 1]);
BZ = permute(reshape(BZ,[nz,nx,ny]).*fbz,[2 3 1]);
xbstar = x0+x;
ybstar = y;
zbstar = z;
xbend = xbstar+(nx-1)*dx;
ybend = ybstar+(ny-1)*dy;
zbend = zbstar+(nz-1)*dz;
xaxis = xbstar:dx:xbend;
yaxis = ybstar:dy:ybend;
zaxis = zbstar:dz:zbend;

