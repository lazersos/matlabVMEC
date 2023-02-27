function [ output_args ] = plot_fieldlines(data,varargin)
%PLOT_FIELDLINES(data,[plottype])  Plots the data from read_mgrid
%   The PLOT_FIELDLINES routine plots data read by READ_FIELDLINES.  There
%   are various plotting options.
%   Options:
%       'basic':        Poincare plot on the first cutplane.
%       '3D':           Same as 'basic' but in 3D space.
%       'color':        Specify color ('color','r')
%       'iota':         Color according to local iota value
%       'cutplane':     Poincare plot on specific cutplane ('cutplane',5)
%       'camera':       Make a camera image by binning poincare points.
%       'camera_AEV30': W7-X AEV30 view (from AEQ21)
%       'strike_2D':    Strike pattern
%       'skip':         Skip this many fieldlines in plot
%       'wall_strike':  Strucutre strike heat map
%       'camview':      Use current view to construct a camer view
%       'phi':          Plot on interpolated value of phi 2D ('phi',0.51)
%       'phi3D':        Plot on interpolated value of phi 3D ('phi3D',0.51)
%       'movie':        Make cutplane movie (.mp4)
%
%   Usage:
%       line_data=read_fieldlines('fieldlines_test.h5');
%       plot_fieldlines(line_data);
%
%   See also read_fieldlines.
%
%   Created by: S. Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    1.22
%   Date:       11/19/14


% Initialize some variables
plottype=0;
phi0=0;
nphi=1;
npoinc = data.npoinc;
nsteps = data.nsteps;
nlines = data.nlines;
line_color='k';
camera = [];
skip=1;
liota=0;

% Handle varargin
if nargin > 1
    i = 1;
    while i < nargin
        switch varargin{i}
            case 'basic'
                plottype=0;
            case '3D'
                plottype=101;
            case 'camera'
                plottype=1;
            case 'camera_frame'
                plottype=2;
            case 'camera_AEV30_old'
                plottype=3;
            case 'camera_AEV30'
                plottype=4;
            case 'camera_AEV30_3D'
                plottype=5;
            case 'cutplane'
                i=i+1;
                nphi=varargin{i};
            case 'movie'
                plottype=11;
            case 'phi'
                i=i+1;
                phi0=varargin{i};
                plottype=102;
            case 'phi3D'
                i=i+1;
                phi0=varargin{i};
                plottype=103;
            case 'strike_2D'
                plottype=6;
            case 'wall_strike'
                plottype=7;
            case 'camera_AEA30'
                plottype=8;
            case 'camera_AEA30_3D'
                plottype=9;
            case 'camview'
                plottype=10;
            case 'color'
                i=i+1;
                line_color=varargin{i};
            case 'resolution'
                i=i+1;
                camera=varargin{i};
            case 'skip'
                i=i+1;
                skip=varargin{i};
            case 'iota'
                liota = 1;
                x = data.R_lines - data.R_lines(1,:);
                y = data.Z_lines - data.Z_lines(1,:);
                theta = atan2(y,x);
                dtheta = diff(theta,[],2); %axis=0?
                dtheta(dtheta<-pi) = dtheta(dtheta<-pi)+2*pi;
                dtheta(dtheta>pi) = dtheta(dtheta>pi)-2*pi;
                theta = abs(cumsum(dtheta,2));
                iota = zeros(size(x,1),1);
                for i = 1:data.nlines
                    p = polyfit(data.PHI_lines(i,1:data.nsteps-1),theta(i,:),1);
                    iota(i) = p(1);
                end
                iota(1) = 2*iota(2)-iota(3);
        end
        i=i+1;
    end
end

switch plottype
    case{0}
        line_dex = nphi:npoinc:nsteps;
        x=data.R_lines(1:skip:nlines,line_dex);
        y=data.Z_lines(1:skip:nlines,line_dex);
        if isfield(data,'rho')
            s=ones(size(x,1),size(x,2));
            c=repmat(data.rho(1:skip:nlines),[1 size(x,2)]);
            scatter(x(:),y(:),s(:).*0.1,c(:),'.');
            caxis([0 data.rho(end-1)]);
        else
            if liota
                t = repmat(iota,size(x,2),1);
                scatter(reshape(x,[],1),reshape(y,[],1),1.0,t);
            else
                %plot(x,y,'.','Color',line_color,'MarkerSize',0.1);
                %colorder = jet(numel(line_dex));
                %colororder(colorder)
                plot(x,y,'k.','MarkerSize',1.0);
            end
        end
        if isfield(data,'Rhc_lines')
            line_dex = nphi:npoinc:size(data.Rhc_lines,2);
            for i=1:size(data.Rhc_lines,1)
                hold on
                plot(data.Rhc_lines(i,line_dex),data.Zhc_lines(i,line_dex),'.r');
                hold off
            end
            hold on
            plot(data.Rhc_lines(1,nphi),data.Zhc_lines(1,nphi),'+r');
            hold off
        end
        axis equal
        xlabel('R [m]')
        ylabel('Z [m]')
    case{101}
        line_dex = nphi:npoinc:nsteps;
        R = data.R_lines(1:skip:nlines,line_dex);
        Z = data.Z_lines(1:skip:nlines,line_dex);
        phi = data.PHI_lines(1,78);
        X = R.*cos(phi);
        Y = R.*sin(phi);
        hold on;
        plot3(X,Y,Z,'.','Color',line_color,'MarkerSize',0.1);
        %for i=1:skip:nlines
        %    hold on;
        %    %plot3(data.X_lines(i,line_dex),data.Y_lines(i,line_dex),data.Z_lines(i,line_dex),'.','Color',line_color,'MarkerSize',0.1);
        %    hold off;
        %end
        axis equal
    case{102}
        phi2 = phi0:max(data.phiaxis):max(max(data.PHI_lines));
        R=zeros(data.nlines,length(phi2));
        Z=zeros(data.nlines,length(phi2));
        for i=1:data.nlines
            n = find(data.R_lines(i,:)==0,1,'first')-1;
            phi = data.PHI_lines(i,1:n);
            n2 = find(phi2 > max(phi),1,'first')-1;
            if isempty(n2); n2=length(phi2)-1; end
            R(i,1:n2) = pchip(phi,data.R_lines(i,1:n),phi2(1:n2));
            Z(i,1:n2) = pchip(phi,data.Z_lines(i,1:n),phi2(1:n2));
        end
        x=R(1:skip:nlines,:);
        y=Z(1:skip:nlines,:);
        if isfield(data,'rho')
            n1 = size(x,1);
            n2 = size(x,2);
            s=ones(n1,n2);
            c=repmat(data.rho(1:skip:nlines),[1 n2]);
            scatter(x(:),y(:),s(:).*0.1,c(:),'.');
        else
            plot(x,y,'.','Color',line_color,'MarkerSize',0.1);
        end
        %plot(R(1:skip:data.nlines,:),Z(1:skip:nlines,:),'.','Color',line_color,'MarkerSize',0.1);
        axis equal
        x1=min(min(x(x>0)));
        x2=max(xlim);
        y1=min(min(y));
        y2=max(max(y));
        d1=x2-x1;
        d2=y2-y1;
        d=max(d1,d2).*0.5;
        c=(x2+x1).*0.5;
        d=d.*1.1;
        xlim([-d d]+c);
        ylim([-d d]);
    case{103}
        phi2 = phi0:2*pi./data.nfp:max(max(data.PHI_lines));
        R=zeros(data.nlines,length(phi2));
        Z=zeros(data.nlines,length(phi2));
        for i=1:skip:data.nlines
            n = find(data.R_lines(i,:)==0,1,'first')-1;
            phi = data.PHI_lines(i,1:n);
            n2 = find(phi2 > max(phi),1,'first')-1;
            if isempty(n2); n2=length(phi2)-1; end
            R(i,1:n2) = pchip(phi,data.R_lines(i,1:n),phi2(1:n2));
            Z(i,1:n2) = pchip(phi,data.Z_lines(i,1:n),phi2(1:n2));
        end
        X = R.*cos(phi0);
        Y = R.*sin(phi0);
        hold on;
        plot3(X',Y',Z','.','Color',line_color,'MarkerSize',0.1);
        axis equal
    case{1}
        if isempty(camera), camera=[1024 1024];end
        line_dex = nphi:npoinc:nsteps;
        r=data.R_lines(:,line_dex);
        z=data.Z_lines(:,line_dex);
        r=reshape(r,[numel(r) 1]);
        z=reshape(z,[numel(z) 1]);
        % Use vessel as constraint
        r=[4.298; r; 6.341];
        z=[-0.76; z; 0.76];
        syn=hist3([r z],camera);
        xb=linspace(min(r),max(r),size(syn,1));
        yb=linspace(min(z),max(z),size(syn,1));
        pixplot(xb,yb,syn)
        caxis([0 max(mean(syn))]);
        colormap bone;
        axis equal;
        set(gca,'Color','black');
        %axis tight;
    case{2}
        if isempty(camera), camera=[1300 1030]; end
        cx=[0.53 0.91];
        cy=[-0.19 0.19];
        line_dex = nphi:npoinc:nsteps;
        r=data.R_lines(:,line_dex);
        z=data.Z_lines(:,line_dex);
        r=reshape(r,[numel(r) 1]);
        z=reshape(z,[numel(z) 1]);
        z = z(r>0.541);
        r = r(r>0.541);
        syn=hist3([r z],camera);
        xb=linspace(min(r),max(r),size(syn,1));
        yb=linspace(min(z),max(z),size(syn,1));
        pixplot(xb,yb,syn)
        caxis([0 max(mean(syn))]);
        colormap bone;
        axis equal;
        xlim(cx);
        ylim(cy);
    case{3} % AEV30 W7-X
        if isempty(camera), camera=[1392 1024]; end
        n_cam = [-0.91206066 -0.37672560 -0.16193573];
        x_cam = [1.69477886 6.12453262 0.64880745];
        a_cam = 32.21906432;
        u_cam = [-0.24104390  0.17308427  0.95495532]; %camroll(14.59570350);
        dex = and(data.R_lines>0,data.PHI_lines <= 2*pi);
        dex = and(dex,data.R_lines>data.raxis(1));
        dex = and(dex,data.R_lines<data.raxis(end));
        dex = and(dex,data.Z_lines>data.zaxis(1));
        dex = and(dex,data.Z_lines<data.raxis(end));
        X   = data.X_lines(dex);
        Y   = data.Y_lines(dex);
        Z   = data.Z_lines(dex);
        % A failed attempt to do a weighting
        %phi = data.PHI_lines;
        %n   = data.nsteps-2;
        %for i=1:data.nlines
        %    phi(i,2:data.nsteps)=1:-1/double(n):0;
        %end
        %phi = phi(dex);
        %phi = phi(2:end);
        % New Stuff
        [x_im,  y_im] = points_to_camera(X(2:end),Y(2:end),Z(2:end),...
            'camera',camera,...
            'fov',a_cam,'camera_pos',x_cam,'camera_normal',n_cam,...
            'camera_up',u_cam);
        x_max = max(x_im);
        y_max = max(y_im);
        x_min = min(x_im);
        y_min = min(y_im);
        syn=hist3([x_im y_im],'nbins',[round(x_max-x_min) round(y_max-y_min)]);
        xb=linspace(x_min,x_max,size(syn,1));
        yb=linspace(y_min,y_max,size(syn,2));
        pixplot(xb,yb,syn)
        caxis([0 max(mean(syn))]);
        set(gcf,'Units','pixels','Position',[1 1 camera]);
        set(gca,'Units','pixels','Color','black','Position',[1 1 camera]);
        xlim([1 camera(1)]);
        ylim([1 camera(2)]);
        colormap hot;
        hold on; plot([333 1336],[376 622],'r','LineWidth',4.0); % plot z=0
    case{4} % AEV30 W7-X
        if isempty(camera), camera=[1392 1024]; end
        n_cam = [-0.91206066 -0.37672560 -0.16193573];
        x_cam = [1.69477886 6.12453262 0.64880745];
        a_cam = 32.21906432;
        u_cam = [-0.24104390  0.17308427  0.95495532]; %camroll(14.59570350);
        syn = zeros(camera);
        for i = 2:data.nsteps
            dex = and(data.R_lines(:,i)>0,data.PHI_lines(:,i) <= 2*pi);
            dex = and(dex,data.R_lines(:,i)>data.raxis(1));
            dex = and(dex,data.R_lines(:,i)<data.raxis(end));
            dex = and(dex,data.Z_lines(:,i)>data.zaxis(1));
            dex = and(dex,data.Z_lines(:,i)<data.raxis(end));
            X   = data.X_lines(dex,i);
            Y   = data.Y_lines(dex,i);
            Z   = data.Z_lines(dex,i);
            if isempty(X), continue; end
            [x_im,  y_im] = points_to_camera(X(2:end),Y(2:end),Z(2:end),...
                'camera',camera,...
                'fov',a_cam,'camera_pos',x_cam,'camera_normal',n_cam,...
                'camera_up',u_cam);
            dex   = x_im <= camera(1);
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = y_im <= camera(2);
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = x_im >= 1;
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = y_im >= 1;
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            if isempty(x_im), continue; end
            x_max = max(x_im);
            y_max = max(y_im);
            x_min = min(x_im);
            y_min = min(y_im);
            syn_temp=hist3([x_im y_im],'nbins',[round(x_max-x_min) round(y_max-y_min)]);
            xb=linspace(x_min,x_max,size(syn_temp,1));
            yb=linspace(y_min,y_max,size(syn_temp,2));
            syn(round(xb),round(yb))=syn_temp./double(i)+syn(round(xb),round(yb));
        end
        % New Stuff
        pixplot(syn)
        caxis([0 max(mean(syn))]);
        set(gcf,'Units','pixels','Position',[1 1 camera]);
        set(gca,'Units','pixels','Color','black','Position',[1 1 camera]);
        xlim([1 camera(1)]);
        ylim([1 camera(2)]);
        colormap hot;
        hold on; plot([333 1336],[376 622],'r','LineWidth',4.0); % plot z=0
    case{5} % AEV30 W7-X (3D plot)
        % OLD STUFF
        n_cam = [-0.91206066 -0.37672560 -0.16193573];
        x_cam = [1.69477886 6.12453262 0.64880745];
        a_cam = 32.21906432;
        u_cam = [-0.24104390  0.17308427  0.95495532]; %camroll(14.59570350);
        dex = and(data.R_lines>0,data.PHI_lines <= 2*pi);
        X   = data.X_lines(dex);
        Y   = data.Y_lines(dex);
        Z   = data.Z_lines(dex);
        campos(x_cam);
        camtarget(x_cam+n_cam);
        camva(a_cam);
        camup(u_cam);
        %axis equal;
        axis off;
        camproj('perspective');
        hold on;
        plot3(X,Y,Z,'.');
    case{6} % Strike_2D
        r=data.R_lines(:,2);
        z=data.Z_lines(:,2);
        r0 = 9.8;
        rr=r-r0;
        u = atan2(z,rr);
        u = u(isfinite(u));
        v = mod(data.PHI_lines(:,2),data.phiaxis(end));
        v = v(isfinite(v));
        u(u<0) = u(u<0) + 2*pi;
        v(v<0) = v(v<0) + data.phiaxis(end);
        %plot(v,u,'.')
        syn=hist3([v u],[100 100]);
        xb=linspace(min(v),max(v),size(syn,1));
        yb=linspace(min(u),max(u),size(syn,1));
        pixplot(xb,yb,syn)
        colormap hot;
        axis square;
        xlim([min(v) max(v)]);
        ylim([min(u) max(u)]);
    case{7} % Wall strike heat map
        dex1 = data.wall_faces(:,1);
        dex2 = data.wall_faces(:,2);
        dex3 = data.wall_faces(:,3);
        V0   = data.wall_vertex(dex3,:)-data.wall_vertex(dex1,:);
        V1   = data.wall_vertex(dex2,:)-data.wall_vertex(dex1,:);
        FNx  = V1(:,2).*V0(:,3)-V1(:,3).*V0(:,2);
        FNy  = V1(:,3).*V0(:,1)-V1(:,1).*V0(:,3);
        FNz  = V1(:,1).*V0(:,2)-V1(:,2).*V0(:,1);
        heat = 2*double(data.wall_strikes)./sqrt(FNx.*FNx+FNy.*FNy+FNz.*FNz);
        output_args{1}=patch('Vertices',data.wall_vertex,'Faces',data.wall_faces,'FaceVertexCData',heat,'LineStyle','none','CDataMapping','scaled','FaceColor','flat');
    case{8} % Glen's port
        figure('Color','white','Position',[1 -100 1024 768]);
        %x_cam=[-6899 5012 -1549]./1000;
        c_cam=[-6899 5012 -1449]./1000;
        n_cam=[16.957,-12.488,12.257];
        n_cam = n_cam./sqrt(sum(n_cam.^2));
        u_cam=[ -0.3501    0.2578    0.7470];
        if isempty(camera), camera=[480 440]; end
        a_cam = 4.0; %50 mm
        X = data.X_lines(:,2);
        Y = data.Y_lines(:,2);
        Z = data.Z_lines(:,2);
        P = atan2(Y,X);
        X = X(P>7*pi/10);
        Y = Y(P>7*pi/10);
        Z = Z(P>7*pi/10);
        P = atan2(Y,X);
        X = X(P<9*pi/10);
        Y = Y(P<9*pi/10);
        Z = Z(P<9*pi/10);
        [x_im,  y_im] = points_to_camera(X,Y,Z,...
            'camera',camera,...
            'fov',a_cam,'camera_pos',c_cam,'camera_normal',n_cam,...
            'camera_up',u_cam);
        dex   = x_im <= camera(1);
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        dex   = y_im <= camera(2);
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        dex   = x_im >= 1;
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        dex   = y_im >= 1;
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        x_max = max(x_im);
        y_max = max(y_im);
        x_min = min(x_im);
        y_min = min(y_im);
        syn_temp=hist3([x_im y_im],'nbins',[round(x_max-x_min) round(y_max-y_min)]);
        xb=linspace(x_min,x_max,size(syn_temp,1));
        yb=linspace(y_min,y_max,size(syn_temp,2));
        syn = zeros(camera);
        syn(round(xb),round(yb))=syn_temp./double(i)+syn(round(xb),round(yb));
        pixplot(syn)
        set(gcf,'Units','pixels','Position',[1 1 camera]);
        set(gca,'Units','pixels','Color','black','Position',[1 1 camera]);
        xlim([1 camera(1)]);
        ylim([1 camera(2)]);
        colormap hot;
    case{9}
        figure('Color','white','Position',[1 -100 1024 768]);
        x_cam=[-6899 5012 -1549]./1000;
        n_cam=[16.957,-12.488,12.257];
        n_cam = n_cam./sqrt(sum(n_cam.^2));
        a_cam = 100.0; %50 mm
        port_data = read_limiter('/p/w7x_sci/vessel/ComponentsDB/ComponentsDB_060616/port_w7x_m3_389.dat');
        ves_data = read_limiter('/p/w7x_sci/vessel/ComponentsDB/ComponentsDB_060616/vessel_w7x_m3_342.dat');
        dex1 = data.wall_faces(:,1);
        dex2 = data.wall_faces(:,2);
        dex3 = data.wall_faces(:,3);
        V0   = data.wall_vertex(dex3,:)-data.wall_vertex(dex1,:);
        V1   = data.wall_vertex(dex2,:)-data.wall_vertex(dex1,:);
        FNx  = V1(:,2).*V0(:,3)-V1(:,3).*V0(:,2);
        FNy  = V1(:,3).*V0(:,1)-V1(:,1).*V0(:,3);
        FNz  = V1(:,1).*V0(:,2)-V1(:,2).*V0(:,1);
        heat = 2*double(data.wall_strikes)./sqrt(FNx.*FNx+FNy.*FNy+FNz.*FNz);
        bad_vertex = [1386 1388 1422 1424];
        faces = data.wall_faces;
        for i=1:length(bad_vertex)
            dex = find(faces(:,1)==bad_vertex(i));
            for j=1:length(dex)
                faces = [faces(1:dex(j)-1,:); faces(dex(j)+1:end,:)];
                heat = [heat(1:dex(j)-1,:); heat(dex(j)+1:end,:)];
                dex = dex -1;
            end
            dex = find(faces(:,2)==bad_vertex(i));
            for j=1:length(dex)
                faces = [faces(1:dex(j)-1,:); faces(dex(j)+1:end,:)];
                heat = [heat(1:dex(j)-1,:); heat(dex(j)+1:end,:)];
                dex = dex -1;
            end
            dex = find(faces(:,3)==bad_vertex(i));
            for j=1:length(dex)
                faces = [faces(1:dex(j)-1,:); faces(dex(j)+1:end,:)];
                heat = [heat(1:dex(j)-1,:); heat(dex(j)+1:end,:)];
                dex = dex -1;
            end
        end
        output_args{1}=patch('Vertices',data.wall_vertex,'Faces',faces,'FaceVertexCData',heat,'LineStyle','none','CDataMapping','scaled','FaceColor','flat');
        colormap hot;
        campos(x_cam);
        camtarget(x_cam+n_cam);
        camva(a_cam);
        hold on;
        ha=plot_divertor(port_data,'solid');
        set(ha,'FaceColor','blue');
        ha=plot_divertor(ves_data,'solid');
        set(ha,'FaceColor','blue');
        camlight left;
        axis equal;
    case{10} % Camview
        if isempty(camera), camera=[1024 768]; end
        x_cam = campos;
        a_cam = camva;
        u_cam = camup;
        t_cam = camtarget;
        n_cam = t_cam-x_cam;
        n_cam = n_cam./sqrt(sum(n_cam.*n_cam));
        syn = zeros(camera);
        X   = data.X_lines(:,2);
        Y   = data.Y_lines(:,2);
        Z   = data.Z_lines(:,2);
        [x_temp] = points_to_camera(X(2:end),Y(2:end),Z(2:end),...
            'camera',camera,...
            'fov',a_cam,'camera_pos',x_cam,'camera_normal',n_cam,...
            'camera_up',u_cam);
        x_im = x_temp(:,1);
        y_im = x_temp(:,2);
        dex   = x_im <= camera(1);
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        dex   = y_im <= camera(2);
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        dex   = x_im >= 1;
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        dex   = y_im >= 1;
        x_im  = x_im(dex);
        y_im  = y_im(dex);
        x_max = max(x_im);
        y_max = max(y_im);
        x_min = min(x_im);
        y_min = min(y_im);
        syn_temp=hist3([x_im y_im],'nbins',[round(x_max-x_min) round(y_max-y_min)]);
        xb=linspace(x_min,x_max,size(syn_temp,1));
        yb=linspace(y_min,y_max,size(syn_temp,2));
        syn(round(xb),round(yb))=syn_temp./double(i)+syn(round(xb),round(yb));

        % Smooth
        sigma = 1; % set sigma to the value you need
        sz = 2*ceil(2.6 * sigma) + 1; % See note below
        mask = fspecial('gauss', sz, sigma);
        syn2 = conv2(syn, mask, 'same');
        syn=syn2; syn2=[];
        % New Stuff
        pixplot(syn)
        caxis([0 max(mean(syn))]);
        set(gcf,'Units','pixels','Position',[1 1 camera]);
        set(gca,'Units','pixels','Color','black','Position',[1 1 camera]);
        xlim([1 camera(1)]);
        ylim([1 camera(2)]);
        colormap hot;
        %hold on; plot([333 1336],[376 622],'r','LineWidth',4.0); % plot z=0
    case{11}
        v = VideoWriter([data.input_extension(1:end-3),'.avi'], 'Motion JPEG AVI');
        open(v);
        f=figure;
        set(gcf,'position',get(0,'ScreenSize'));
        for nphi = 1:npoinc
            line_dex = nphi:npoinc:nsteps;
            x=data.R_lines(1:skip:nlines,line_dex);
            y=data.Z_lines(1:skip:nlines,line_dex);
            plot(x,y,'.','Color',line_color,'MarkerSize',0.1);
            axis equal
            title('Poicare plot')
            xlabel('R [m]')
            ylabel('Z [m]')
            %xlim([1.5 1.9])
            %ylim([-0.2 0.3])
            xlim([4.5 6.5])
            ylim([-1.5 1.5])
            frame = getframe(f); %for writing
            writeVideo(v,frame);
        end
        close(v);
end
end
