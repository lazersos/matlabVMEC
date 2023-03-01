function plot_mgrid(data,extcur,varargin)
%PLOT_MGRID(data,extcur,[plottype])  Plots the data from read_mgrid
%   The PLOT_MGRID routine plots data read by READ_MGRID.  There are
%   various plotting options.
%   Options:
%       'basic':    Plots the total field on an mgrid plane.
%                   'cutplane' option controls which plane (default=1)
%       'total':    Plots 3 components and vector plot for each slice in
%                   the mgrid file.
%       '3dgrid':   Shows the 3D grid Planes.
%       'modB':     |B| for a given cutplane
%
%   Usage:
%       mgrid_data=read_mgrid('mgrid.test');
%       extcur=[1.2e4 1.2e4 1.2e4 -3.5e3 1.1e6 -2.5e4];
%       plot_mgrid(mgrid_data,extcur);          %Plot Total Field (phi=0)
%
%   See also read_mgrid.
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Version:    1.0
%   Date:       10/19/10

% Handle extcur
if ~(size(extcur,2) == data.nextcur)
    disp(' - ERROR: Extcur size mismatch');
    return
end
% Handle varargin
plottype=0;
cutplane=1;
if nargin > 2
    for i=1:nargin-2
        switch varargin{i}
            case 'basic'
                plottype=0;
            case 'total'
                plottype=1;
            case '3dgrid'
                plottype=2;
            case 'modB'
                plottype=3;
            case 'cutplane'
                i=i+1;
                cutplane=varargin{i};
        end
    end
end
% Setup Plotting Window
fig=figure('Position',[1 1 1920 1080]);
% Multiply B-fields by extcur
bx=0.*data.bx;
by=0.*bx;
bz=0.*bx;
bphi=0.*bx;
br=0.*bx;
for i=1:data.nextcur
    bx(:,:,:,i)=data.bx(:,:,:,i).*extcur(i);
    by(:,:,:,i)=data.by(:,:,:,i).*extcur(i);
    bz(:,:,:,i)=data.bz(:,:,:,i).*extcur(i);
    bphi(:,:,:,i)=data.bphi(:,:,:,i)*extcur(i);
    br(:,:,:,i)=data.br(:,:,:,i).*extcur(i);
end
% Handle Creating total B-field
switch plottype
    case {0,1,3}
        bxt=zeros(data.nr,data.nz,data.nphi);
        byt=zeros(data.nr,data.nz,data.nphi);
        bzt=zeros(data.nr,data.nz,data.nphi);
        brt=zeros(data.nr,data.nz,data.nphi);
        bphit=zeros(data.nr,data.nz,data.nphi);
        for i=1:data.nextcur
            bxt=bxt+bx(:,:,:,i);
            byt=byt+by(:,:,:,i);
            bzt=bzt+bz(:,:,:,i);
            brt=brt+br(:,:,:,i);
            bphit=bphit+bphi(:,:,:,i);
        end
        bminr=min(min(min(brt)));
        bminz=min(min(min(bzt)));
        bminp=min(min(min(bphit)));
        bmin=min([bminr bminp bminz]);
        bmaxr=max(max(max(brt)));
        bmaxz=max(max(max(bzt)));
        bmaxp=max(max(max(bphit)));
        bmax=max([bmaxr bmaxz bmaxp]);
end
% Now Handles plots
switch plottype
    case 1 % Pan through cuts
        startx=3.2:.1:4.6;
        starty=0.*startx;
        haxes=subplot(2,2,1);
        raxis2d=repmat(data.raxis',[1 data.nz]);
        zaxis2d=repmat(data.zaxis,[data.nr 1]);
        bmag=sqrt(brt.*brt+bphit.*bphit+bzt.*bzt);
        for i=1:data.nphi
            % First plot Toroidal B-Field
            subplot(2,2,1);
            hplot1=pcolor(raxis2d,zaxis2d,bphit(:,:,i));
            set(hplot1,'EdgeColor','none');
            xlabel('R-Axis');
            ylabel('Z-Axis');
            title('Toroidal Field');
            colorbar
            caxis([bmin bmax]);
            axis image
            % Second plot Radial Field
            subplot(2,2,2);
            hplot2=pcolor(raxis2d,zaxis2d,brt(:,:,i));
            set(hplot2,'EdgeColor','none');
            xlabel('R-Axis');
            ylabel('Z-Axis');
            title('Radial Field');
            colorbar
            caxis([bmin bmax]);
            axis image
            % Third plot vertical Field
            subplot(2,2,3);
            hplot3=pcolor(raxis2d,zaxis2d,bzt(:,:,i));
            set(hplot3,'EdgeColor','none');
            xlabel('R-Axis');
            ylabel('Z-Axis');
            title('Vertical Field');    
            colorbar
            caxis([bmin bmax]);
            axis image
            % Now plot combine field
            subplot(2,2,4);
            hplot4=pcolor(raxis2d,zaxis2d,bphit(:,:,i));
            set(hplot4,'EdgeColor','none');
            hold on
            quiver(raxis2d,zaxis2d,brt(:,:,i),bzt(:,:,i),'Color','black');
            hold off
            colorbar
            xlabel('R-Axis');
            ylabel('Z-Axis');
            title('MGRID B-Field at \phi=0');
            caxis([bmin bmax]);
            axis image
            pause(1.0);
        end
    case 2 % Plot 3D grid cutplanes
        % Because the actual grids are very fine we abstract to a
        % 20x20 grid for visualization
        rtemp=data.rmin:(data.rmax-data.rmin)/19:data.rmax;
        ztemp=data.zmin:(data.zmax-data.zmin)/19:data.zmax;
        % First we need to 2d arrays for r and z
        raxis2d=repmat(rtemp',[1 20]);
        zaxis2d=repmat(ztemp,[20 1]);
        % Next we need to map to X(nr,nphi,nz) and Y(nr,nphi,nz)
        % Z is the same on each cut plane
        x2d=zeros(20,data.nphi,20);
        y2d=zeros(20,data.nphi,20);
        for i=1:20
            for j=1:data.nphi
                x2d(i,j,:)=raxis2d(i,:).*cos(data.phi(j));
                y2d(i,j,:)=raxis2d(i,:).*sin(data.phi(j));
            end
        end
        % Now for each surface we create a create the patch object
        % The verex list will run over nr
        for i=1:data.nphi
            vertex=[0 0 0];
            for j=1:20
                for k=1:20
                    vertex=[vertex; x2d(k,i,j) y2d(k,i,j) zaxis2d(k,j)];
                end
            end
            vertex=vertex(2:size(vertex,1),:);
            faces=[1 2 2+20 1+20];
            for j=1:19
                for k=1:19
                    index=j+20*(k-1);
                    faces=[faces; index index+1 index+1+20 index+20];
                end
            end
            faces=faces(2:size(faces,1),:);
            hold on
            patch('Vertices',vertex,'Faces',faces,'FaceColor','none')
            hold off
        end
        xlim([0 data.rmax]);
        ylim([0 data.rmax]);
        zlim([data.zmin data.zmax]);
        axis square
        view(3)
    case 3 % Plot |B|
        raxis2d=repmat(data.raxis',[1 data.nz]);
        zaxis2d=repmat(data.zaxis,[data.nr 1]);
        bmag=sqrt(brt.*brt+bphit.*bphit+bzt.*bzt);
        bmin=min(min(min(bmag)));
        bmax=max(max(max(bmag)));
        hplot1=pcolor(raxis2d,zaxis2d,bmag(:,:,cutplane));
        set(hplot1,'EdgeColor','none');
        xlabel('R-Axis');
        ylabel('Z-Axis');
        title('|B| Field');
        colorbar
        %caxis([bmin bmax]);
        axis image
        
    case 0 %Total Plot on phi=0 plane only
        haxes=subplot(2,2,1);
        raxis2d=repmat(data.raxis',[1 data.nz]);
        zaxis2d=repmat(data.zaxis,[data.nr 1]);
        bmag=sqrt(brt.*brt+bphit.*bphit+bzt.*bzt);
        % First plot Toroidal B-Field
        subplot(2,2,1);
        hplot1=pcolor(raxis2d,zaxis2d,bphit(:,:,cutplane));
        set(hplot1,'EdgeColor','none');
        xlabel('R-Axis');
        ylabel('Z-Axis');
        title('Toroidal Field');
        colorbar
        caxis([bmin bmax]);
        axis image
        % Second plot Radial Field
        subplot(2,2,2);
        hplot2=pcolor(raxis2d,zaxis2d,brt(:,:,cutplane));
        set(hplot2,'EdgeColor','none');
        xlabel('R-Axis');
        ylabel('Z-Axis');
        title('Radial Field');
        colorbar
        caxis([bmin bmax]);
        axis image
        % Third plot vertical Field
        subplot(2,2,3);
        hplot3=pcolor(raxis2d,zaxis2d,bzt(:,:,cutplane));
        set(hplot3,'EdgeColor','none');
        xlabel('R-Axis');
        ylabel('Z-Axis');
        title('Vertical Field');
        colorbar
        caxis([bmin bmax]);
        axis image
        % Now plot combine field
        subplot(2,2,4);
        hplot4=pcolor(raxis2d,zaxis2d,bphit(:,:,cutplane));
        set(hplot4,'EdgeColor','none');
        hold on
        quiver(raxis2d,zaxis2d,brt(:,:,cutplane),bzt(:,:,cutplane),'Color','black');
        hold off
        colorbar
        caxis([bmin bmax]);
        xlabel('R-Axis');
        ylabel('Z-Axis');
        title('MGRID B-Field at \phi=0');
        %caxis([bmin bmax]);
        axis image
end
end

