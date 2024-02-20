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
% zeta:     Magnetic toroidal angle (zeta)
% s:        Vector of surfaces to plot
% color:    Array of colors to plot on surface.
% 'STL':    Will output to an STL file.
%
% Exmaple Usage (data assumed to have at least 10 flux surfaces)
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
%      z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
%      hpatch=isotoro(r,z,zeta,[2 10]);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       2.05

new_color=[];
loutputtoobj = 0; % Will output STL files if set to 1

if nargin > 4
    j=1;
    while j<=length(varargin)
        if isstr(varargin{j})
            switch varargin{j}
                case{'STL','stl'}
                    loutputtoobj=1;
            end
        elseif isnumeric(varargin{j})
            mincolor=min(min(min(varargin{j})));
            maxcolor=max(max(max(varargin{j})));
            cmap=colormap;
            csize=size(cmap,1);
            caxis([mincolor maxcolor]);
            if size(r,3) == 1
                for i=1:size(r,3)
                    new_color(:,:,i)=varargin{j}(:,:,1);
                end
            else
                new_color=varargin{j};
            end
        end
        j=j+1;
    end
    %cfacedata=[1 0 0];
end

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
cfacedata=zeros(nvertex,ns);
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
if ~isempty(new_color)
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
                cfacedata(ivertex,k)=new_color(s(k),i,j);
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
scale = 1.0; x0=0; y0=0; z0=0;

% Handle plotting a single surface
if ns==1
    try
        app = evalin('base','app');
        hpatch=patch(app.UIAxes, 'Vertices',vertex,'Faces',faces, 'FaceVertexCData',cfacedata, 'FaceColor','none', 'EdgeColor','flat');
        %set(hpatch, 'FaceAlpha',0.5, 'EdgeColor','black');
%         hpatch.AmbientStrength = 0.55;
%         hpatch.DiffuseStrength = 0.7;
%         hpatch.SpecularStrength = 0;
        assignin('base','xData',squeeze(vertex(:,1,1))); %Expects there to only be one plotted flux surface
        assignin('base', 'yData',squeeze(vertex(:,2,1)));
        assignin('base','zData',squeeze(vertex(:,3,1)));
        assignin('base','Faces',faces);
    catch
        hpatch=patch('Vertices',vertex,'Faces',faces,'FaceVertexCData',cfacedata);
    end
    if ~isempty(new_color)
        set(hpatch,'EdgeColor','none','FaceColor','interp','CDataMapping','scaled');
    else        
        set(hpatch,'EdgeColor','none','FaceColor','red');
    end
    
    if (loutputtoobj ==1)
        vertex2(:,1) = vertex(:,1)-x0;
        vertex2(:,2) = vertex(:,2)-y0;
        vertex2(:,3) = vertex(:,3)-z0;
        %scale = 0.1350;
        TR = triangulation(faces(:,[2 1 3]),vertex2*scale);
        stlwrite(TR, 'isotoro.stl');
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
        hpatch(i)=patch('Vertices',vertex(:,:,i),'Faces',faces(:,:,i),'FaceVertexCData',cfacedata(:,i));
        set(hpatch(i),'EdgeColor','none',...
            'FaceAlpha',0.5,'FaceColor','red');
        alpha(hpatch(i),0.6*(ns-i+1)/ns);
    end
    hold off
end
% Clean up the plot with some settings. 
%set(gca,'LineStyle','none');
set(gcf,'Renderer','OpenGL'); lighting gouraud; %camlight left;
% Note:  We use zbuffer as it supports camlight (old)
%set(gcf,'Renderer','zbuffer'); lighting phong; camlight left;
view(3);
title('Flux Surfaces');
axis equal
% Output the patch surfaces
varargout{1}=hpatch;
end
