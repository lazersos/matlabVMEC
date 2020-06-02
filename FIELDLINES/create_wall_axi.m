function wall_data=create_wall_axi(rin,zin,ntor,name)
%CREATE_WALL_AXI Create a 3D wall from an axisymmetric trace
%   The CREATE_WALL_AXI function creates a wall structure from an
%   axisymmetric trace defined by rin and zin.  Here rin and zin define the
%   R/Z datapoints of the axisymmetric shape and ntor defines the number of
%   toroidal slices.  The rin and zin arrays must terminate where they
%   start.  The name parameter is writen to the machine name field.
%
%   Example:
%       rvec = 0.2:0.01:0.4;
%       zvec = -0.25:0.01:0.25;
%       r=[rvec zvec.*0.0+max(rvec) fliplr(rvec)  zvec.*0+min(rvec)];
%       z=[rvec.*0+min(zvec) zvec rvec.*0+max(zvec) fliplr(zvec)];
%       wall_data = create_wall_axi(r,z,255,'TEST VESSEL');
%       plot_divertor(wall_data,'wire');
%       write_wall(wall_data,'test_vessel.dat');
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           4/28/20

wall_data=[];
dphi = 2*pi/(ntor-1); 
phi = 0:dphi:2*pi;
l = length(rin);
n = 1;
vertex=[];
faces=[];
% this is the cheap way
for i=1:length(phi)
    for j=1:l
        x = rin(j).*cos(phi(i));
        y = rin(j).*sin(phi(i));
        vertex=[vertex; x y zin(j)];
    end
end
n=1;
for i=1:length(phi)-1 % last point is first
    for j=1:l-1 % last point is first
        d1 = n;
        d2 = n+1;
        d3 = n+l;
        d4 = d2+l;
        faces=[faces; d1 d3 d2];
        faces=[faces; d3 d4 d2];
        n = n + 1;
    end
    n=n+1;
    % Add last duplicate point
end
% check with plot
[faces, vertex]=reducepatch(faces,vertex,1);
wall_data=[];
wall_data.machine=name;
wall_data.date=datestr(now,'mm-dd-yy');
wall_data.datatype='limiter_trimesh';
wall_data.nvertex = size(vertex,1);
wall_data.nfaces = size(faces,1);
wall_data.coords = vertex';
wall_data.faces = faces';
end

