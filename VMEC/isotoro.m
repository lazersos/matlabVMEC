function [varargout]=isotoro(r,z,zeta,s,varargin)
% ISOTORO(r,z,zeta,s,[color]) Plots a 3d isosurface of the
% flux surface indexed by s.  The surface will be plotted in an existing
% axis object (or a new one if none exist).
%
% hpatch=ISOTORO(r,z,zeta,[s,color]) Plots multiple 3d
% isosurfaces of the flux surface indexed by s.  The surfaces will be
% plotted in an existing axis object (or a new one if none exist).  A
% surface will be plotted for each value in [s].  The alpha channel of
% each surface will be chosen so the inner surfaces are more visible
% than the outter.  Will return handles to the patch surfaces it plots.
%
% This function plots a 3d isosurface of the flux surface s:
% ISOTORO(r,z,zeta,s)
% Inputs
% r:        Radial position r(s,theta,zeta)
% z:        Vertical position z(s,theta,zeta)
% theta:    Magnetic polodial angle (theta)
% s:        Vector of surfaces to plot
% color:    Array of colors to plot on surface.
%
% Exmaple Usage (data assumed to have at least 10 flux surfaces)
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
%      z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
%      hpatch=isotoro(r,z,zeta,[2 10]);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       2.05
if nargin > 5
    disp('-- Error too many arguments');
    return
elseif nargin >4
    mincolor=min(min(min(varargin{1})));
    maxcolor=max(max(max(varargin{1})));
    cmap=colormap;
    csize=size(cmap,1);
    colorbar;
    caxis([mincolor maxcolor]);
    if size(r,3) == 1
        for i=1:size(r,3)
            new_color(:,:,i)=varargin{1}(:,:,1);
        end
    else
        new_color=varargin{1};
    end 
else
    %cfacedata=[1 0 0];
end

loutputtoobj = 0; % Will output STL files if set to 1

maxzeta=max(zeta);
ns=size(squeeze(s),2);
ntheta=size(r,2);
nzeta=size(r,3);
% Handle plotting 2D equilibria
if nzeta==1
    nzeta=60;
    zeta=0:2*pi/59:2*pi;
    for i=1:nzeta
        new_r(:,:,i)=r(:,:,1);
        new_z(:,:,i)=z(:,:,1);
    end
else
    new_r=r;
    new_z=z;
end
nvertex=ntheta*nzeta;
vertex=zeros(nvertex,3,ns);
faces=zeros(ntheta*(nzeta-1),3,ns);
cfacedata=zeros(nvertex,3,ns);
% Now we calculate the vertex and face data for the patch surface
for k=1:ns
    ivertex = 1;
    ifaces = 1;
    for j=1:nzeta
        for i=1:ntheta-1
            % X Position
            vertex(ivertex,1,k)=new_r(s(k),i,j)*cos(zeta(j));
            % Y Position
            vertex(ivertex,2,k)=new_r(s(k),i,j)*sin(zeta(j));
            % Z Position
            vertex(ivertex,3,k)=new_z(s(k),i,j);
            if (j==nzeta)
                
            elseif (i==ntheta)
                faces(ifaces,1,k)=ivertex;
                faces(ifaces,2,k)=ivertex-ntheta+1;
                faces(ifaces,3,k)=ivertex+1;
                
            else
                faces(ifaces,1,k)=ivertex;
                faces(ifaces,2,k)=ivertex+1; %theta+1
                faces(ifaces,3,k)=ivertex+ntheta+1; %zeta+1,theta+1
            end
            ifaces=ifaces+1;
            if (j==nzeta)
                
            elseif (i==ntheta)
                faces(ifaces,1,k)=ivertex;
                faces(ifaces,3,k)=ivertex+1;
                faces(ifaces,2,k)=ivertex+ntheta;
            else
                faces(ifaces,1,k)=ivertex;
                faces(ifaces,2,k)=ivertex+ntheta+1; %zeta+1,theta+1
                faces(ifaces,3,k)=ivertex+ntheta;
            end
            ifaces=ifaces+1;
            ivertex=ivertex+1;
            
        end
        % X Position
        vertex(ivertex,1,k)=new_r(s(k),1,j)*cos(zeta(j));
        % Y Position
        vertex(ivertex,2,k)=new_r(s(k),1,j)*sin(zeta(j));
        % Z Position
        vertex(ivertex,3,k)=new_z(s(k),1,j);
        ivertex=ivertex+1;
    end
end
% Get Color Information
if nargin > 4
    for k=1:ns
        ivertex=1;
        for j=1:nzeta
            for i=1:ntheta
                % Color Data
                temp=fix(csize*...
                    (new_color(s(k),i,j)-mincolor)...
                    /(maxcolor-mincolor));
                if temp==0
                    temp=1;
                end
                cfacedata(ivertex,1,k)=cmap(temp,1);
                cfacedata(ivertex,2,k)=cmap(temp,2);
                cfacedata(ivertex,3,k)=cmap(temp,3);
                % Increment ivertex
                ivertex=ivertex+1;
            end
        end
    end
end
% Get rid of the superflous faces if the torus isn't a full torus
%if (maxzeta ~= 2*pi)
%    faces=faces(1:(nzeta-1)*ntheta,:,:);
%end

% check for 3D printer output
scale = 1.0;
if (loutputtoobj == 1)
    maker_scale=150; %mm
    xmax = max(vertex(:,1,ns));
    ymax = max(vertex(:,2,ns));
    zmax = max(vertex(:,3,ns));
    xmin = min(vertex(:,1,ns));
    ymin = min(vertex(:,2,ns));
    zmin = min(vertex(:,3,ns));
    xlen = xmax-xmin;
    ylen = ymax-ymin;
    zlen = zmax-zmin;
    scale = min([maker_scale/xlen, maker_scale/ylen, maker_scale/zlen]);
    x0   = xmin+0.5*xlen;
    y0   = ymin+0.5*ylen;
    z0   = zmin;
end

% Handle plotting a single surface
if ns==1
    hpatch=patch('Vertices',vertex,'Faces',faces,'FaceVertexCData',cfacedata);
    if nargin > 4
        set(hpatch,'EdgeColor','none','FaceColor','interp','CDataMapping','direct');
    else        
        set(hpatch,'EdgeColor','none','FaceColor','red');
    end
    
    if (loutputtoobj ==1)
        vertex2(:,1) = vertex(:,1)-x0;
        vertex2(:,2) = vertex(:,2)-y0;
        vertex2(:,3) = vertex(:,3)-z0;
        %scale = 0.1350;
        stlwrite2('isotoro.stl',faces(:,[2 1 3]),vertex2*scale);
    end
    
    
else %Multiple surfaces
    hpatch=zeros(ns);
    hold on
    for i=1:ns
        if (loutputtoobj == 1)
            vertex2(:,1,i) = vertex(:,1,i)-x0;
            vertex2(:,2,i) = vertex(:,2,i)-y0;
            vertex2(:,3,i) = vertex(:,3,i)-z0;
            stlwrite(['isotoro_ns' num2str(i,'%2.2d')  '.stl'],faces(:,[2 1 3],i),vertex2(:,:,i)*scale);
        end
        hpatch(i)=patch('Vertices',vertex(:,:,i),'Faces',faces(:,:,i),'FaceVertexCData',cfacedata(:,:,i));
        set(hpatch(i),'EdgeColor','none',...
            'FaceAlpha',0.5,'FaceColor','red');
        alpha(hpatch(i),0.6*(ns-i+1)/ns);
    end
    hold off
end
% Clean up the plot with some settings. 
set(gcf,'Renderer','OpenGL'); lighting gouraud; %camlight left;
% Note:  We use zbuffer as it supports camlight (old)
%set(gcf,'Renderer','zbuffer'); lighting phong; camlight left;
view(3);
title('Flux Surfaces');
axis equal
% Output the patch surfaces
varargout{1}=hpatch;
end
