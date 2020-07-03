function wall_out = truncate_wall(wall_data,rmax,rmin,zmax,zmin)
%TRUNCATE_WALL Truncates a wall strcuture
%   The TRUNCATE_WALL function truncates a triangulated wall by rmax, rmin,
%   zmax and zmin.  The truncation method removes triangles completely
%   outside the domain.  Triangles with one vertex outside the domain are
%   converted into three triangles.  Triangles with two verticies outside
%   the domain are shrunk.  A new wall structure is returned.
%
%   Example
%       wall_data=read_wall('test.dat');
%       wall_out = truncate_wall(wall_data,6.5,3.5,1.5,-1.5);
%
%  See also read_wall, reduce_wall, plot_wall.
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           6/17/20


wall_out=[];
faces=wall_data.faces;
verts=wall_data.coords;

d1=faces(1,:);
d2=faces(2,:);
d3=faces(3,:);

xyz1=verts(:,d1);
xyz2=verts(:,d2);
xyz3=verts(:,d3);
r1=sqrt(xyz1(1,:).^2+xyz1(2,:).^2);
r2=sqrt(xyz2(1,:).^2+xyz2(2,:).^2);
r3=sqrt(xyz3(1,:).^2+xyz3(2,:).^2);
z1=xyz1(3,:);
z2=xyz2(3,:);
z3=xyz3(3,:);

% Remove bad points
rmaxbad = [r1>rmax;r2>rmax;r3>rmax];
rminbad = [r1<rmin;r2<rmin;r3<rmin];
zmaxbad = [z1>zmax;z2>zmax;z3>zmax];
zminbad = [z1<zmin;z2<zmin;z3<zmin];
lrmax = sum(rmaxbad)==3;
lrmin = sum(rminbad)==3;
lzmax = sum(zmaxbad)==3;
lzmin = sum(zminbad)==3;
lbad = (lrmax+lrmin+lzmax+lzmin)>0; % If any one is bad
xyz1=xyz1(:,~lbad);
xyz2=xyz2(:,~lbad);
xyz3=xyz3(:,~lbad);

% Now recalce the helpers
r1=sqrt(xyz1(1,:).^2+xyz1(2,:).^2);
r2=sqrt(xyz2(1,:).^2+xyz2(2,:).^2);
r3=sqrt(xyz3(1,:).^2+xyz3(2,:).^2);
z1=xyz1(3,:);
z2=xyz2(3,:);
z3=xyz3(3,:);
rmaxbad = [r1>rmax;r2>rmax;r3>rmax];
rminbad = [r1<rmin;r2<rmin;r3<rmin];
zmaxbad = [z1>zmax;z2>zmax;z3>zmax];
zminbad = [z1<zmin;z2<zmin;z3<zmin];

xyznew1=[];
xyznew2=[];
xyznew3=[];
% Handle rmax
maskd = sum(rmaxbad)==2;
maska = sum(rmaxbad)>0;
mask1 = rmaxbad(1,:)==1;
mask2 = rmaxbad(2,:)==1;
mask3 = rmaxbad(3,:)==1;
mask1(maskd) = 0;
mask2(maskd) = 0;
mask3(maskd) = 0;
maska(maskd) = 0;
dex1=find(mask1==1);
dex2=find(mask2==1);
dex3=find(mask3==1);
dexa=find(maska==1);
% xyz1
for i=1:length(dex1)
    %pt1
    dx = xyz1(:,dex1(i))-xyz2(:,dex1(i));
    dr = r1(dex1(i))-r2(dex1(i));
    l = (rmax-r2(dex1(i)))./dr;
    xnew1=xyz2(:,dex1(i))+l*dx;
    %pt2
    dx = xyz1(:,dex1(i))-xyz3(:,dex1(i));
    dr = r1(dex1(i))-r3(dex1(i));
    l = (rmax-r3(dex1(i)))./dr;
    xnew2=xyz3(:,dex1(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz2(:,dex1(i))+xyz3(:,dex1(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz2(:,dex1(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz3(:,dex1(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% xyz2
for i=1:length(dex2)
    %pt1
    dx = xyz2(:,dex2(i))-xyz1(:,dex2(i));
    dr = r2(dex2(i))-r1(dex2(i));
    l = (rmax-r1(dex2(i)))./dr;
    xnew1=xyz1(:,dex2(i))+l*dx;
    %pt2
    dx = xyz2(:,dex2(i))-xyz3(:,dex2(i));
    dr = r2(dex2(i))-r3(dex2(i));
    l = (rmax-r3(dex2(i)))./dr;
    xnew2=xyz3(:,dex2(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex2(i))+xyz3(:,dex2(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex2(i))];
    xyznew2=[xyznew2 xnew1];
    xyznew3=[xyznew3 xnew3];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xyz3(:,dex2(i))];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xnew3];
end
% xyz3
for i=1:length(dex3)
    %pt1
    dx = xyz3(:,dex3(i))-xyz1(:,dex3(i));
    dr = r3(dex3(i))-r1(dex3(i));
    l = (rmax-r1(dex3(i)))./dr;
    xnew1=xyz1(:,dex3(i))+l*dx;
    %pt2
    dx = xyz3(:,dex3(i))-xyz2(:,dex3(i));
    dr = r3(dex3(i))-r2(dex3(i));
    l = (rmax-r2(dex3(i)))./dr;
    xnew2=xyz2(:,dex3(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex3(i))+xyz2(:,dex3(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex3(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz2(:,dex3(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% Now we do the 2point removal
% Redefine dexn so we can use
maskd = sum(rmaxbad)==2;
mask1 = rmaxbad(1,:)==1;
mask2 = rmaxbad(2,:)==1;
mask3 = rmaxbad(3,:)==1;
mask1(~maskd) = 0;
mask2(~maskd) = 0;
mask3(~maskd) = 0;
dex1=find(and(mask2==1,mask3==1));
dex2=find(and(mask1==1,mask3==1));
dex3=find(and(mask1==1,mask2==1));
% Keep point 1
for i=1:length(dex1)
    %pt1
    dx = xyz2(:,dex1(i))-xyz1(:,dex1(i));
    dr = r2(dex1(i))-r1(dex1(i));
    l = (rmax-r1(dex1(i)))./dr;
    xnew1=xyz1(:,dex1(i))+l*dx;
    %pt2
    dx = xyz3(:,dex1(i))-xyz1(:,dex1(i));
    dr = r3(dex1(i))-r1(dex1(i));
    l = (rmax-r1(dex1(i)))./dr;
    xnew2=xyz1(:,dex1(i))+l*dx;
    xyz2(:,dex1(i)) = xnew1;
    xyz3(:,dex1(i)) = xnew2;
end
% Keep point 2
for i=1:length(dex2)
    %pt1
    dx = xyz1(:,dex2(i))-xyz2(:,dex2(i));
    dr = r1(dex2(i))-r2(dex2(i));
    l = (rmax-r2(dex2(i)))./dr;
    xnew1=xyz2(:,dex2(i))+l*dx;
    %pt2
    dx = xyz3(:,dex2(i))-xyz2(:,dex2(i));
    dr = r3(dex2(i))-r2(dex2(i));
    l = (rmax-r2(dex2(i)))./dr;
    xnew2=xyz2(:,dex2(i))+l*dx;
    xyz1(:,dex2(i)) = xnew1;
    xyz3(:,dex2(i)) = xnew2;
end
% Keep point 3
for i=1:length(dex3)
    %pt1
    dx = xyz1(:,dex3(i))-xyz3(:,dex3(i));
    dr = r1(dex3(i))-r3(dex3(i));
    l = (rmax-r3(dex3(i)))./dr;
    xnew1=xyz3(:,dex3(i))+l*dx;
    %pt2
    dx = xyz2(:,dex3(i))-xyz3(:,dex3(i));
    dr = r2(dex3(i))-r3(dex3(i));
    l = (rmax-r3(dex3(i)))./dr;
    xnew2=xyz3(:,dex3(i))+l*dx;
    xyz1(:,dex3(i)) = xnew1;
    xyz2(:,dex3(i)) = xnew2;
end
% Now reconstruct new array
xyz1(:,dexa)=[];xyz2(:,dexa)=[];xyz3(:,dexa)=[];
xyz1=[xyz1, xyznew1];
xyz2=[xyz2, xyznew2];
xyz3=[xyz3, xyznew3];

% Now recalce the helpers
r1=sqrt(xyz1(1,:).^2+xyz1(2,:).^2);
r2=sqrt(xyz2(1,:).^2+xyz2(2,:).^2);
r3=sqrt(xyz3(1,:).^2+xyz3(2,:).^2);
z1=xyz1(3,:);
z2=xyz2(3,:);
z3=xyz3(3,:);
% Remove bad points
rminbad = [r1<rmin;r2<rmin;r3<rmin];
zmaxbad = [z1>zmax;z2>zmax;z3>zmax];
zminbad = [z1<zmin;z2<zmin;z3<zmin];
lrmin = sum(rminbad)==3;
lzmax = sum(zmaxbad)==3;
lzmin = sum(zminbad)==3;
lbad = (lrmin+lzmax+lzmin)>0; % If any one is bad
xyz1=xyz1(:,~lbad);
xyz2=xyz2(:,~lbad);
xyz3=xyz3(:,~lbad);
r1=sqrt(xyz1(1,:).^2+xyz1(2,:).^2);
r2=sqrt(xyz2(1,:).^2+xyz2(2,:).^2);
r3=sqrt(xyz3(1,:).^2+xyz3(2,:).^2);
rminbad = [r1<rmin;r2<rmin;r3<rmin];

%%%%%%%%%%%%%%%%  RMIN
xyznew1=[];
xyznew2=[];
xyznew3=[];
maskd = sum(rminbad)==2;
maska = sum(rminbad)>0;
mask1 = rminbad(1,:)==1;
mask2 = rminbad(2,:)==1;
mask3 = rminbad(3,:)==1;
mask1(maskd) = 0;
mask2(maskd) = 0;
mask3(maskd) = 0;
maska(maskd) = 0;
dex1=find(mask1==1);
dex2=find(mask2==1);
dex3=find(mask3==1);
dexa=find(maska==1);
% xyz1
for i=1:length(dex1)
    %pt1
    dx = xyz1(:,dex1(i))-xyz2(:,dex1(i));
    dr = r1(dex1(i))-r2(dex1(i));
    l = (rmin-r2(dex1(i)))./dr;
    xnew1=xyz2(:,dex1(i))+l*dx;
    %pt2
    dx = xyz1(:,dex1(i))-xyz3(:,dex1(i));
    dr = r1(dex1(i))-r3(dex1(i));
    l = (rmin-r3(dex1(i)))./dr;
    xnew2=xyz3(:,dex1(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz2(:,dex1(i))+xyz3(:,dex1(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz2(:,dex1(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz3(:,dex1(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% xyz2
for i=1:length(dex2)
    %pt1
    dx = xyz2(:,dex2(i))-xyz1(:,dex2(i));
    dr = r2(dex2(i))-r1(dex2(i));
    l = (rmin-r1(dex2(i)))./dr;
    xnew1=xyz1(:,dex2(i))+l*dx;
    %pt2
    dx = xyz2(:,dex2(i))-xyz3(:,dex2(i));
    dr = r2(dex2(i))-r3(dex2(i));
    l = (rmin-r3(dex2(i)))./dr;
    xnew2=xyz3(:,dex2(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex2(i))+xyz3(:,dex2(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex2(i))];
    xyznew2=[xyznew2 xnew1];
    xyznew3=[xyznew3 xnew3];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xyz3(:,dex2(i))];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xnew3];
end
% xyz3
for i=1:length(dex3)
    %pt1
    dx = xyz3(:,dex3(i))-xyz1(:,dex3(i));
    dr = r3(dex3(i))-r1(dex3(i));
    l = (rmin-r1(dex3(i)))./dr;
    xnew1=xyz1(:,dex3(i))+l*dx;
    %pt2
    dx = xyz3(:,dex3(i))-xyz2(:,dex3(i));
    dr = r3(dex3(i))-r2(dex3(i));
    l = (rmin-r2(dex3(i)))./dr;
    xnew2=xyz2(:,dex3(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex3(i))+xyz2(:,dex3(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex3(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz2(:,dex3(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% Now we do the 2point removal
% Redefine dexn so we can use
maskd = sum(rminbad)==2;
mask1 = rminbad(1,:)==1;
mask2 = rminbad(2,:)==1;
mask3 = rminbad(3,:)==1;
mask1(~maskd) = 0;
mask2(~maskd) = 0;
mask3(~maskd) = 0;
dex1=find(and(mask2==1,mask3==1));
dex2=find(and(mask1==1,mask3==1));
dex3=find(and(mask1==1,mask2==1));
% Keep point 1
for i=1:length(dex1)
    %pt1
    dx = xyz2(:,dex1(i))-xyz1(:,dex1(i));
    dr = r2(dex1(i))-r1(dex1(i));
    l = (rmin-r1(dex1(i)))./dr;
    xnew1=xyz1(:,dex1(i))+l*dx;
    %pt2
    dx = xyz3(:,dex1(i))-xyz1(:,dex1(i));
    dr = r3(dex1(i))-r1(dex1(i));
    l = (rmin-r1(dex1(i)))./dr;
    xnew2=xyz1(:,dex1(i))+l*dx;
    xyz2(:,dex1(i)) = xnew1;
    xyz3(:,dex1(i)) = xnew2;
end
% Keep point 2
for i=1:length(dex2)
    %pt1
    dx = xyz1(:,dex2(i))-xyz2(:,dex2(i));
    dr = r1(dex2(i))-r2(dex2(i));
    l = (rmin-r2(dex2(i)))./dr;
    xnew1=xyz2(:,dex2(i))+l*dx;
    %pt2
    dx = xyz3(:,dex2(i))-xyz2(:,dex2(i));
    dr = r3(dex2(i))-r2(dex2(i));
    l = (rmin-r2(dex2(i)))./dr;
    xnew2=xyz2(:,dex2(i))+l*dx;
    xyz1(:,dex2(i)) = xnew1;
    xyz3(:,dex2(i)) = xnew2;
end
% Keep point 3
for i=1:length(dex3)
    %pt1
    dx = xyz1(:,dex3(i))-xyz3(:,dex3(i));
    dr = r1(dex3(i))-r3(dex3(i));
    l = (rmin-r3(dex3(i)))./dr;
    xnew1=xyz3(:,dex3(i))+l*dx;
    %pt2
    dx = xyz2(:,dex3(i))-xyz3(:,dex3(i));
    dr = r2(dex3(i))-r3(dex3(i));
    l = (rmin-r3(dex3(i)))./dr;
    xnew2=xyz3(:,dex3(i))+l*dx;
    xyz1(:,dex3(i)) = xnew1;
    xyz2(:,dex3(i)) = xnew2;
end
% Now reconstruct new array
xyz1(:,dexa)=[];xyz2(:,dexa)=[];xyz3(:,dexa)=[];
xyz1=[xyz1, xyznew1];
xyz2=[xyz2, xyznew2];
xyz3=[xyz3, xyznew3];

% Now recalc the helpers
z1=xyz1(3,:);
z2=xyz2(3,:);
z3=xyz3(3,:);
% Remove bad points
zmaxbad = [z1>zmax;z2>zmax;z3>zmax];
zminbad = [z1<zmin;z2<zmin;z3<zmin];
lzmax = sum(zmaxbad)==3;
lzmin = sum(zminbad)==3;
lbad = (lzmax+lzmin)>0; % If any one is bad
xyz1=xyz1(:,~lbad);
xyz2=xyz2(:,~lbad);
xyz3=xyz3(:,~lbad);
z1=xyz1(3,:);
z2=xyz2(3,:);
z3=xyz3(3,:);
zmaxbad = [z1>zmax;z2>zmax;z3>zmax];
zminbad = [z1<zmin;z2<zmin;z3<zmin];

%%%%%%%%%%%%%%%%  ZMAX
xyznew1=[];
xyznew2=[];
xyznew3=[];
maskd = sum(zmaxbad)==2;
maska = sum(zmaxbad)>0;
mask1 = zmaxbad(1,:)==1;
mask2 = zmaxbad(2,:)==1;
mask3 = zmaxbad(3,:)==1;
mask1(maskd) = 0;
mask2(maskd) = 0;
mask3(maskd) = 0;
maska(maskd) = 0;
dex1=find(mask1==1);
dex2=find(mask2==1);
dex3=find(mask3==1);
dexa=find(maska==1);
% xyz1
for i=1:length(dex1)
    %pt1
    dx = xyz1(:,dex1(i))-xyz2(:,dex1(i));
    dz = z1(dex1(i))-z2(dex1(i));
    l = (zmax-z2(dex1(i)))./dz;
    xnew1=xyz2(:,dex1(i))+l*dx;
    %pt2
    dx = xyz1(:,dex1(i))-xyz3(:,dex1(i));
    dz = z1(dex1(i))-z3(dex1(i));
    l = (zmax-z3(dex1(i)))./dz;
    xnew2=xyz3(:,dex1(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz2(:,dex1(i))+xyz3(:,dex1(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz2(:,dex1(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz3(:,dex1(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% xyz2
for i=1:length(dex2)
    %pt1
    dx = xyz2(:,dex2(i))-xyz1(:,dex2(i));
    dz = z2(dex2(i))-z1(dex2(i));
    l = (zmax-z1(dex2(i)))./dz;
    xnew1=xyz1(:,dex2(i))+l*dx;
    %pt2
    dx = xyz2(:,dex2(i))-xyz3(:,dex2(i));
    dz = z2(dex2(i))-z3(dex2(i));
    l = (zmax-z3(dex2(i)))./dz;
    xnew2=xyz3(:,dex2(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex2(i))+xyz3(:,dex2(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex2(i))];
    xyznew2=[xyznew2 xnew1];
    xyznew3=[xyznew3 xnew3];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xyz3(:,dex2(i))];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xnew3];
end
% % xyz3
for i=1:length(dex3)
    %pt1
    dx = xyz3(:,dex3(i))-xyz1(:,dex3(i));
    dz = z3(dex3(i))-z1(dex3(i));
    l = (zmax-z1(dex3(i)))./dz;
    xnew1=xyz1(:,dex3(i))+l*dx;
    %pt2
    dx = xyz3(:,dex3(i))-xyz2(:,dex3(i));
    dz = z3(dex3(i))-z2(dex3(i));
    l = (zmax-z2(dex3(i)))./dz;
    xnew2=xyz2(:,dex3(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex3(i))+xyz2(:,dex3(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex3(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz2(:,dex3(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% Now we do the 2point removal
% Redefine dexn so we can use
maskd = sum(zmaxbad)==2;
mask1 = zmaxbad(1,:)==1;
mask2 = zmaxbad(2,:)==1;
mask3 = zmaxbad(3,:)==1;
mask1(~maskd) = 0;
mask2(~maskd) = 0;
mask3(~maskd) = 0;
dex1=find(and(mask2==1,mask3==1));
dex2=find(and(mask1==1,mask3==1));
dex3=find(and(mask1==1,mask2==1));
% Keep point 1
for i=1:length(dex1)
    %pt1
    dx = xyz2(:,dex1(i))-xyz1(:,dex1(i));
    dz = z2(dex1(i))-z1(dex1(i));
    l = (zmax-z1(dex1(i)))./dz;
    xnew1=xyz1(:,dex1(i))+l*dx;
    %pt2
    dx = xyz3(:,dex1(i))-xyz1(:,dex1(i));
    dz = z3(dex1(i))-z1(dex1(i));
    l = (zmax-z1(dex1(i)))./dz;
    xnew2=xyz1(:,dex1(i))+l*dx;
    xyz2(:,dex1(i)) = xnew1;
    xyz3(:,dex1(i)) = xnew2;
end
% Keep point 2
for i=1:length(dex2)
    %pt1
    dx = xyz1(:,dex2(i))-xyz2(:,dex2(i));
    dz = z1(dex2(i))-z2(dex2(i));
    l = (zmax-z2(dex2(i)))./dz;
    xnew1=xyz2(:,dex2(i))+l*dx;
    %pt2
    dx = xyz3(:,dex2(i))-xyz2(:,dex2(i));
    dz = z3(dex2(i))-z2(dex2(i));
    l = (zmax-z2(dex2(i)))./dz;
    xnew2=xyz2(:,dex2(i))+l*dx;
    xyz1(:,dex2(i)) = xnew1;
    xyz3(:,dex2(i)) = xnew2;
end
% Keep point 3
for i=1:length(dex3)
    %pt1
    dx = xyz1(:,dex3(i))-xyz3(:,dex3(i));
    dz = z1(dex3(i))-z3(dex3(i));
    l = (zmax-z3(dex3(i)))./dz;
    xnew1=xyz3(:,dex3(i))+l*dx;
    %pt2
    dx = xyz2(:,dex3(i))-xyz3(:,dex3(i));
    dz = z2(dex3(i))-z3(dex3(i));
    l = (zmax-z3(dex3(i)))./dz;
    xnew2=xyz3(:,dex3(i))+l*dx;
    xyz1(:,dex3(i)) = xnew1;
    xyz2(:,dex3(i)) = xnew2;
end
% Now reconstruct new array
xyz1(:,dexa)=[];xyz2(:,dexa)=[];xyz3(:,dexa)=[];
xyz1=[xyz1, xyznew1];
xyz2=[xyz2, xyznew2];
xyz3=[xyz3, xyznew3];




% Now recalc the helpers
z1=xyz1(3,:);
z2=xyz2(3,:);
z3=xyz3(3,:);
% Remove bad points
zmaxbad = [z1>zmax;z2>zmax;z3>zmax];
zminbad = [z1<zmin;z2<zmin;z3<zmin];
lzmax = sum(zmaxbad)==3;
lzmin = sum(zminbad)==3;
lbad = (lzmax+lzmin)>0; % If any one is bad
xyz1=xyz1(:,~lbad);
xyz2=xyz2(:,~lbad);
xyz3=xyz3(:,~lbad);
z1=xyz1(3,:);
z2=xyz2(3,:);
z3=xyz3(3,:);
zminbad = [z1<zmin;z2<zmin;z3<zmin];

%%%%%%%%%%%%%%%%  ZMIN
xyznew1=[];
xyznew2=[];
xyznew3=[];
maskd = sum(zminbad)==2;
maska = sum(zminbad)>0;
mask1 = zminbad(1,:)==1;
mask2 = zminbad(2,:)==1;
mask3 = zminbad(3,:)==1;
mask1(maskd) = 0;
mask2(maskd) = 0;
mask3(maskd) = 0;
maska(maskd) = 0;
dex1=find(mask1==1);
dex2=find(mask2==1);
dex3=find(mask3==1);
dexa=find(maska==1);
% xyz1
for i=1:length(dex1)
    %pt1
    dx = xyz1(:,dex1(i))-xyz2(:,dex1(i));
    dz = z1(dex1(i))-z2(dex1(i));
    l = (zmin-z2(dex1(i)))./dz;
    xnew1=xyz2(:,dex1(i))+l*dx;
    %pt2
    dx = xyz1(:,dex1(i))-xyz3(:,dex1(i));
    dz = z1(dex1(i))-z3(dex1(i));
    l = (zmin-z3(dex1(i)))./dz;
    xnew2=xyz3(:,dex1(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz2(:,dex1(i))+xyz3(:,dex1(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz2(:,dex1(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz3(:,dex1(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% xyz2
for i=1:length(dex2)
    %pt1
    dx = xyz2(:,dex2(i))-xyz1(:,dex2(i));
    dz = z2(dex2(i))-z1(dex2(i));
    l = (zmin-z1(dex2(i)))./dz;
    xnew1=xyz1(:,dex2(i))+l*dx;
    %pt2
    dx = xyz2(:,dex2(i))-xyz3(:,dex2(i));
    dz = z2(dex2(i))-z3(dex2(i));
    l = (zmin-z3(dex2(i)))./dz;
    xnew2=xyz3(:,dex2(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex2(i))+xyz3(:,dex2(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex2(i))];
    xyznew2=[xyznew2 xnew1];
    xyznew3=[xyznew3 xnew3];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xyz3(:,dex2(i))];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew2];
    xyznew3=[xyznew3 xnew3];
end
% % xyz3
for i=1:length(dex3)
    %pt1
    dx = xyz3(:,dex3(i))-xyz1(:,dex3(i));
    dz = z3(dex3(i))-z1(dex3(i));
    l = (zmin-z1(dex3(i)))./dz;
    xnew1=xyz1(:,dex3(i))+l*dx;
    %pt2
    dx = xyz3(:,dex3(i))-xyz2(:,dex3(i));
    dz = z3(dex3(i))-z2(dex3(i));
    l = (zmin-z2(dex3(i)))./dz;
    xnew2=xyz2(:,dex3(i))+l*dx;
    %pt3
    xnew3=0.5.*(xyz1(:,dex3(i))+xyz2(:,dex3(i)));
    % Now create new triangles
    xyznew1=[xyznew1 xyz1(:,dex3(i))];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew1];
    xyznew1=[xyznew1 xnew3];
    xyznew2=[xyznew2 xyz2(:,dex3(i))];
    xyznew3=[xyznew3 xnew2];
    xyznew1=[xyznew1 xnew1];
    xyznew2=[xyznew2 xnew3];
    xyznew3=[xyznew3 xnew2];
end
% Now we do the 2point removal
% Redefine dexn so we can use
maskd = sum(zminbad)==2;
mask1 = zminbad(1,:)==1;
mask2 = zminbad(2,:)==1;
mask3 = zminbad(3,:)==1;
mask1(~maskd) = 0;
mask2(~maskd) = 0;
mask3(~maskd) = 0;
dex1=find(and(mask2==1,mask3==1));
dex2=find(and(mask1==1,mask3==1));
dex3=find(and(mask1==1,mask2==1));
% Keep point 1
for i=1:length(dex1)
    %pt1
    dx = xyz2(:,dex1(i))-xyz1(:,dex1(i));
    dz = z2(dex1(i))-z1(dex1(i));
    l = (zmin-z1(dex1(i)))./dz;
    xnew1=xyz1(:,dex1(i))+l*dx;
    %pt2
    dx = xyz3(:,dex1(i))-xyz1(:,dex1(i));
    dz = z3(dex1(i))-z1(dex1(i));
    l = (zmin-z1(dex1(i)))./dz;
    xnew2=xyz1(:,dex1(i))+l*dx;
    xyz2(:,dex1(i)) = xnew1;
    xyz3(:,dex1(i)) = xnew2;
end
% Keep point 2
for i=1:length(dex2)
    %pt1
    dx = xyz1(:,dex2(i))-xyz2(:,dex2(i));
    dz = z1(dex2(i))-z2(dex2(i));
    l = (zmin-z2(dex2(i)))./dz;
    xnew1=xyz2(:,dex2(i))+l*dx;
    %pt2
    dx = xyz3(:,dex2(i))-xyz2(:,dex2(i));
    dz = z3(dex2(i))-z2(dex2(i));
    l = (zmin-z2(dex2(i)))./dz;
    xnew2=xyz2(:,dex2(i))+l*dx;
    xyz1(:,dex2(i)) = xnew1;
    xyz3(:,dex2(i)) = xnew2;
end
% Keep point 3
for i=1:length(dex3)
    %pt1
    dx = xyz1(:,dex3(i))-xyz3(:,dex3(i));
    dz = z1(dex3(i))-z3(dex3(i));
    l = (zmin-z3(dex3(i)))./dz;
    xnew1=xyz3(:,dex3(i))+l*dx;
    %pt2
    dx = xyz2(:,dex3(i))-xyz3(:,dex3(i));
    dz = z2(dex3(i))-z3(dex3(i));
    l = (zmin-z3(dex3(i)))./dz;
    xnew2=xyz3(:,dex3(i))+l*dx;
    xyz1(:,dex3(i)) = xnew1;
    xyz2(:,dex3(i)) = xnew2;
end
% Now reconstruct new array
xyz1(:,dexa)=[];xyz2(:,dexa)=[];xyz3(:,dexa)=[];
xyz1=[xyz1, xyznew1];
xyz2=[xyz2, xyznew2];
xyz3=[xyz3, xyznew3];

%Convert back to wall format
x1x2x3=[xyz1(1,:); xyz2(1,:); xyz3(1,:)];
y1y2y3=[xyz1(2,:); xyz2(2,:); xyz3(2,:)];
z1z2z3=[xyz1(3,:); xyz2(3,:); xyz3(3,:)];
fig=figure;
p_obj=patch('Xdata',x1x2x3,'YData',y1y2y3,'ZData',z1z2z3,'FaceColor','none');
wall_out.date=datestr(now,'mm-dd-yyyy');
wall_out.machine=['Created by truncate_wall' ...
    ' RMIN: ' num2str(rmin,'%5.2f') ...
    ' RMAX: ' num2str(rmax,'%5.2f') ...
    ' ZMIN: ' num2str(zmin,'%5.2f') ...
    ' ZMAX: ' num2str(zmax,'%5.2f') ];
wall_out.coords=p_obj.Vertices';
wall_out.faces=p_obj.Faces';
wall_out.datatype='limiter_trimesh';
wall_out.nvertex=size(wall_out.coords,2);
wall_out.nfaces=size(wall_out.faces,2);
delete(p_obj);
close(fig);

end

