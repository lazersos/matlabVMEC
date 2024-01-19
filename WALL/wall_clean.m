function wall_out = wall_clean(wall_in)
%WALL_CLEAN(wall) Cleans a wall model removing bad cells.
%   WALL_CLEAN looks for bad triangles (zero size) and removes them.
%
%   Example:
%       ves_bad=read_wall('wall_vessel.dat');
%       ves_clean=wall_clean(ves_bad);
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           04/15/2021

wall_out = [];
if isempty(wall_in)
    return;
else
    wall_out = wall_in;
end
% Area Check
area = sum(wall_out.FN.^2);
toosmall = find(area<=0);
% Inverse denom check
toosmall = [toosmall find(wall_out.invDenom == Inf)];
toosmall = unique(toosmall);
if ~isempty(toosmall)
    disp(strcat('Found ',num2str(length(toosmall),' %i'),' bad triangles'));
    wall_out.faces(:,toosmall)=[];
    wall_out.A(:,toosmall)=[];
    wall_out.V0(:,toosmall)=[];
    wall_out.V1(:,toosmall)=[];
    wall_out.FN(:,toosmall)=[];
    wall_out.DOT00(:,toosmall)=[];
    wall_out.DOT01(:,toosmall)=[];
    wall_out.DOT11(:,toosmall)=[];
    wall_out.invDenom(:,toosmall)=[];
    wall_out.d(:,toosmall)=[];
    wall_out.nfaces = size(wall_out.faces,2);
end

return;

end

