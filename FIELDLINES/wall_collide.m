function [XW, YW, ZW, ik] = wall_collide(wall,x0,y0,z0,x1,y1,z1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
XW=[];
YW=[];
ZW=[];
ik=-1;
if ~isfield(wall,'FN')
    wall = wall_collide_prep(wall);
end

r0 = zeros(size(wall.FN));
dr = zeros(size(wall.FN));
r0(1,:) = x0;
r0(2,:) = y0;
r0(3,:) = z0;
dr(1,:) = x1-x0;
dr(2,:) = y1-y0;
dr(3,:) = z1-z0;
alpha = dot(wall.FN,dr);
beta  = dot(wall.FN,r0);
t     = (wall.d - beta)./alpha;
t(alpha ==0)=500;
t(t<=0) = 500;
t(t>1) = 500;
if all(and(t>1,t<=0)), return; end
t2 = repmat(t,[3 1]);
V2 = r0 + t2.*dr-wall.A;
DOT02 = dot(wall.V0,V2);
DOT12 = dot(wall.V1,V2);
alpha = ((wall.DOT11.*DOT02)-(wall.DOT01.*DOT12)).*wall.invDenom;
beta = ((wall.DOT00.*DOT12)-(wall.DOT01.*DOT02)).*wall.invDenom;
dex1 = alpha<0;
dex2 = beta<0;
dex3 = (alpha+beta)>1;
%alpha(alpha<0) = 2;
%beta(beta<0) = 2;
t(dex1) = 1E5;
t(dex2) = 1E5;
t(dex3) = 1E5;
dex = min(t);
if dex < 2
    ik = find(t==dex);
    XW = V2(1,ik)+wall.A(1,ik);
    YW = V2(2,ik)+wall.A(2,ik);
    ZW = V2(3,ik)+wall.A(3,ik);
end
return;


end

function wall_out = wall_collide_prep(wall_in)
wall_out = wall_in;
vertex = wall_in.coords;
face   = wall_in.faces;
dex1 = face(1,:);
dex2 = face(2,:);
dex3 = face(3,:);
wall_out.A  = vertex(:,dex1);
wall_out.V0 = vertex(:,dex3)-vertex(:,dex1);
wall_out.V1 = vertex(:,dex2)-vertex(:,dex1);
wall_out.FN = cross(wall_out.V1,wall_out.V0);
wall_out.DOT00=dot(wall_out.V0,wall_out.V0);
wall_out.DOT01=dot(wall_out.V0,wall_out.V1);
wall_out.DOT11=dot(wall_out.V1,wall_out.V1);
wall_out.invDenom = 1.0./(wall_out.DOT00.*wall_out.DOT11-wall_out.DOT01.*wall_out.DOT01);
wall_out.d = dot(wall_out.FN,wall_out.A);
return;
end


