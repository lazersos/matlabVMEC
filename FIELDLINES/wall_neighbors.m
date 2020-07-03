function data = wall_neighbors(faces,subdex)
%WALL_NEIGHBORS Calculates the faces which have shared verticies
%   The WELL_NEIGHBORS function takes a (3xN) array of indicies (faces) and
%   an array of sub indexes (subdex).  It returns a cell array of length
%   equal to the length of subdex. In each cell is a list of the face
%   indices which share at least one common index. If an empty array is
%   passed to subdex then the entire length of faces is considered.
%
% Example usage
%      data=read_wall('test.dat');
%      neighbors = wall_neighbors(data.faces,[]);
%
% Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
% Version:       1.00

data=cell(1,size(faces,2)); % Needs to be a cell array
if isempty(subdex)
    subdex = 1:size(faces,2);
end
for i=subdex
    %disp(num2str(i));
    d1 = faces(1,i);
    d2 = faces(2,i);
    d3 = faces(3,i);
    dex1 = ceil(find(faces==d1)/3)';
    dex2 = ceil(find(faces==d2)/3)';
    dex3 = ceil(find(faces==d3)/3)';
    dex = [dex1 dex2 dex3];
    data{i} = unique(dex);
end
return;
end

