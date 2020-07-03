function wall_out = reduce_wall(wall_data,dr)
%REDUCE_WALL Reduces triangles to a given length dr
%   The REDUCE_WALL routine reduces a wall model using successive bisection
%   of the longest leg of each triangle.  Triangles with legs longer than
%   dr are successvily divided in two along their longest leg until all
%   legs have lengths below dr. A new wall structure is returned.
%
%   Example
%       wall_data=read_wall('test.dat');
%       wall_out = reduce_wall(wall_data,0.05);
%
%  See also read_wall, truncate_wall, plot_wall.
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

dx21=xyz2-xyz1;
dx31=xyz3-xyz1;
dx23=xyz2-xyz3;

%N=cross(dx21,dx31);
%A=sqrt(sum(N.*N));
d21=sqrt(sum(dx21.*dx21));
d31=sqrt(sum(dx31.*dx31));
d23=sqrt(sum(dx23.*dx23));
[dmax,imax]=max([d21;d31;d23]);
lbig = dmax>dr;

% Loop over in a stupid way
while any(lbig)
    imax(~lbig) = 0; % if you're not big go away
    mask1  = imax ==1;
    mask2  = imax ==2;
    mask3  = imax ==3;
    
    % Do dx21
    xyznew1=[];
    xyznew2=[];
    xyznew3=[];
    d1  = find(mask1==1);
    for i=1:length(d1)
        j=d1(i);
        x0    = xyz1(:,j);
        dx    = dx21(:,j);
        xf    = xyz3(:,j);
        nstep=2;
        dl    = dx./nstep;
        for k=1:nstep
            xnew = x0+dl;
            xyznew1=[xyznew1 x0];
            xyznew2=[xyznew2 xnew];
            xyznew3=[xyznew3 xf];
            x0 = xnew;
        end
    end
     
    % Do dx31
    d2  = find(mask2==1);
    for i=1:length(d2)
        j=d2(i);
        x0    = xyz1(:,j);
        dx    = dx31(:,j);
        xf    = xyz2(:,j);
        nstep=2;
        dl    = dx./nstep;
        for k=1:nstep
            xnew = x0+dl;
            xyznew1=[xyznew1 x0];
            xyznew2=[xyznew2 xf];
            xyznew3=[xyznew3 xnew];
            x0 = xnew;
        end
    end
    
    % Do dx23
    d3  = find(mask3==1);
    for i=1:length(d3)
        j=d3(i);
        x0    = xyz3(:,j);
        dx    = dx23(:,j);
        xf    = xyz1(:,j);
        nstep=2;
        dl    = dx./nstep;
        for k=1:nstep
            xnew = x0+dl;
            xyznew1=[xyznew1 xf];
            xyznew2=[xyznew2 xnew];
            xyznew3=[xyznew3 x0];
            x0 = xnew;
        end
    end
       
    % Now update the arrays
    dt=[d1 d2 d3];
    xyz1(:,dt)=[];
    xyz2(:,dt)=[];
    xyz3(:,dt)=[];
    xyz1=[xyz1, xyznew1];
    xyz2=[xyz2, xyznew2];
    xyz3=[xyz3, xyznew3];
    
    % Recalc helpers
    dx21=xyz2-xyz1;
    dx31=xyz3-xyz1;
    dx23=xyz2-xyz3;
    d21=sqrt(sum(dx21.*dx21));
    d31=sqrt(sum(dx31.*dx31));
    d23=sqrt(sum(dx23.*dx23));
    [dmax,imax]=max([d21;d31;d23]);
    lbig = dmax>dr;
      
end 

%Convert back to wall format
x1x2x3=[xyz1(1,:); xyz2(1,:); xyz3(1,:)];
y1y2y3=[xyz1(2,:); xyz2(2,:); xyz3(2,:)];
z1z2z3=[xyz1(3,:); xyz2(3,:); xyz3(3,:)];
fig=figure;
p_obj=patch('Xdata',x1x2x3,'YData',y1y2y3,'ZData',z1z2z3,'FaceColor','none');
wall_out.date=datestr(now,'mm-dd-yyyy');
wall_out.machine=['Created by reduce_wall' ...
    ' dr: ' num2str(dr,'%5.2f') ];
wall_out.coords=p_obj.Vertices';
wall_out.faces=p_obj.Faces';
wall_out.datatype='limiter_trimesh';
wall_out.nvertex=size(wall_out.coords,2);
wall_out.nfaces=size(wall_out.faces,2);
delete(p_obj);
close(fig);

end

