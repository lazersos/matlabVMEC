function wall3 = wall_add(wall1,wall2)
%WALL_ADD(wall1,wall2) Combine two wall models and return new model
%   WALL_ADD takes two wall structures as read by READ_WALL and returns a
%   third wall structure which contains the elements of the first two.
%
%   Example:
%       ves_data=read_wall('wall_vessel.dat');
%       div_data=read_wall('wall_divertor.dat');
%       full_wall=wall_add(ves_data,div_data);
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           10/22/2019

if isempty(wall1) && isempty(wall2)
    wall3=[];
    return;
elseif isempty(wall1)
    wall3=wall2;
    return;
elseif isempty(wall2)
    wall3=wall1;
    return;
else
    wall3=wall1;
end
wall3.coords = [wall3.coords wall2.coords];
wall3.faces  = [wall3.faces wall2.faces+wall3.nvertex];
wall3.nfaces = size(wall3.faces,2);
wall3.nvertex= size(wall3.coords,2);
return;

end

