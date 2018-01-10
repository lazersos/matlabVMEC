function make_stellopt_movie_new(varargin)
%MAKE_STELLOPT_MOVIE_new(movietype[,vessel_data]) Create stellopt movie
%   MAKE_STELLOPT_MOVIE(movietype) Creates a movie from stellopt output
%   files.  The routine supports 3 resolutions (1080, 720, 480; 1080 is
%   default).  If vessel data is provided then all VMEC plots will include
%   that data.  The plot options are:
%       'opt'      (default) R-Z plot of the VMEC vaccum surfaces at phi=0
%       'te'         Electron Temperature plot.
%       'ti'         Ion Temperature plot.
%       'vphi'       Toroidal Velocity Plot.
%       'ne'         Electron Density plot.
%       'neline'     Line integrated Electron density plot.
%       'teline'     Line integrated Electron density plot.
%       'tiline'     Line integrated Electron density plot.
%       'mag'        Magnetic diagnostics.
%       'mse'        MSE diagnostics.
%       'iota'       Transform diagnostics.
%       'coils'      EXTCUR diagnostics.
%       'bootstrap'  Bootstrap current diagnostics.
%       'helicity'   Helicity diagnostics.
%       'neo'        Neoclassical transport diagnostics.
%       'bnorm'      Normal fields from COILOPT++.
%       'D3D'        Set views suitable for D3D
%       'ITER'       Set views suitable for ITER
%       'ncsx'       Set views suitable for NCSX
%       'lhd'        Set views suitable for LHD (default)
%       'last'       Only plot last equilibria
%       'movie'      Make an AVI file (slower)
%       'realspace'  Plot profiles in real space (R) (default)
%       'zspace'     Plot profiles in real space abs(Z)
%       'fluxspace'  Plot profiles in flux space (s)
%
%   Example:
%       ves_data=read_vessel('vessel.dat');
%       make_stellopt_movie('1080','recon',ves_data);
%
%   See also plot_stellopt, read_stellopt, and read_vessel.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.5
%   Date:           1/06/13

% Handle Input Arguments
plottype='opt';
axistype='realspace';
plotves=0;
ves_data=[];
restype='1080';
warning('off');
text_handle=[];
lold_data=0;
isote=0;
lmovie=0;
llast=0;
mag_specs=[];
coil_data=[];
stel_data=[];
ylim_mag = [-1.0 1.0];
neline_plot=[3 6];
faraday_plot=[3 6];
thom_plot=[3 6];
grey_val = 0.5;
ltop_down = 0;
lrz_mag = 1;
%%%%%%  LHD
vmec_xmin=2.5; vmec_xmax=5;
vmec_ymin=-1.5; vmec_ymax=1.5;
% Set up plot offsets
dx=(vmec_xmax-vmec_xmin)/20;  % 5%
dy=(vmec_ymax-vmec_ymin)/100;  % 1%
if nargin>0
    i=1;
    while i<=nargin
        if isfield(varargin{i},'datatype')
            switch varargin{i}.datatype
                case 'vessel'
                    plotves=1;
                    ves_data=varargin{i};
                case 'fluxloop'
                    mag_specs.fluxloops=varargin{i};
                case 'coil_data'
                    coil_data=varargin{i};
                case 'stellopt_new'
                    stel_data=varargin{i};
            end
        else
            switch varargin{i}
                case {'opt','mse','te','ne','ti','mag','coils','bootstrap',...
                        'helicity','neo','balloon','neline','iota',...
                        'faraday','sxr','vphi','bnorm','teline','tiline'}
                    plottype=varargin{i};
                case {'1080','720','480'}
                    restype=varargin{i};
                case {'realspace','fluxspace','zspace','topdown'}
                    axistype=varargin{i};
                case {'old'}
                    lold_data=1;
                case {'W7X','w7x'}
                    vmec_xmin=4.5; vmec_xmax=6.5;
                    vmec_ymin=-1.00; vmec_ymax=1.00;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                    thom_plot=[2 5];
                case {'area51','AREA51'}
                    vmec_xmin=0.0; vmec_xmax=2.0;
                    vmec_ymin=-1.0; vmec_ymax=1.0;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                case {'ncsx','NCSX'}
                    vmec_xmin=0.75; vmec_xmax=2.25;
                    vmec_ymin=-0.75; vmec_ymax=0.75;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                    ylim_mag = [-1 1]*5E-3;
                    thom_plot=[1 4];
                    faraday_plot=[1 4];
                    neline_plot=[1 4];
                case {'LHD'}
                    vmec_xmin=2.5; vmec_xmax=5;
                    vmec_ymin=-1.5; vmec_ymax=1.5;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                    ylim_mag = [-1 1]*5E-3;
                    thom_plot=[3 6];
                    neline_plot=[1 4];
                    faraday_plot=[1 4];
                case {'HSX'}
                    vmec_xmin=0.85; vmec_xmax=1.55;
                    vmec_ymin=-0.35; vmec_ymax=0.35;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                    ylim_mag = [-1 1]*1E-3;
                    thom_plot=[3 6];
                    neline_plot=[1 4];
                    faraday_plot=[1 4];
                    lrz_mag = 0;
                case {'D3D','DIIID'}
                    vmec_xmin=0.0; vmec_xmax=3.2;
                    vmec_ymin=-1.6; vmec_ymax=1.6;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                    ylim_mag = [-1 1];
                    neline_plot=[3 6];
                    thom_plot=[3 6];
                case {'NSTX'}
                    vmec_xmin=-1.0; vmec_xmax=2.5;
                    vmec_ymin=-1.75; vmec_ymax=1.75;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                    ylim_mag = [-1 1];
                    neline_plot=[3 6];
                    thom_plot=[3 6];
                case {'ITER'}
                    vmec_xmin=0.5; vmec_xmax=10.5;
                    vmec_ymin=-5; vmec_ymax=5;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                    ylim_mag = [-3 0].*1E2;
                    neline_plot=[3 6];
                    thom_plot=[3 6];
                case {'QPS'}
                    vmec_xmin=0.5; vmec_xmax=1.5;
                    vmec_ymin=-0.5; vmec_ymax=0.5;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                case {'isote'}
                    isote = 1;
                case {'movie'}
                    lmovie = 1;
                case {'last'}
                    llast = 1;
            end
        end
        i=i+1;
    end
else
    disp('ERROR: Specify movie type');
    return
end

% Setup Figure Window
%   1080 1920x1080
%   720 1280x720
%   480 640x480
%   text_pos = [x y w h]
if strcmp(restype,'720')
    win_res=[1 1 1280 720];
    text_pos=[660 90 200 100];
    fontsiz1=13;
    arrowsiz=6;
    markersiz=8;
    chisqsiz=10;
elseif strcmp(restype,'480')
    win_res=[1 1 640 480];
    text_pos=[295 75 150 50];
    fontsiz1=8;
    arrowsiz=2;
    markersiz=6;
    chisqsiz=10;
else
    win_res=[1 1 1920 1080];
    text_pos=[1025 100 250 200];
    fontsiz1=18;
    arrowsiz=10;
    markersiz=10;
    chisqsiz=20;
end
fig=figure('DoubleBuffer','on','Position',win_res,...
        'Color','black','BackingStore','on','MenuBar','none',...
        'Name',plottype,'InvertHardcopy','off');
set(gca,'nextplot','replacechildren','XColor','white',...
        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
% Load data
stel_filelist=dir('stellopt.*');
if isempty(stel_data)
    stel_data=read_stellopt(stel_filelist.name);
end
wout_filename='wout_*.*.nc';
filelist=dir(wout_filename);
nfiles=length(stel_data.iter);
if nfiles > 1
    if stel_data.iter(nfiles) == stel_data.iter(nfiles-1)
        nfiles=nfiles-1;
    end
end


% Creat movie filename
temp=stel_filelist.name;
dex1=strfind(temp,'.');
dex2=length(temp);
movname=strcat(plottype,temp(dex1:dex2));
if strcmp(axistype,'fluxspace')
    movname=strcat(movname,'_flux');
end

% Setup Movie File
if (lmovie)
    mov=VideoWriter([movname '_movie.avi']);
    mov.Quality= 100;
    mov.FrameRate=1;
    open(mov);
end

% Get coils points
if ~isempty(coil_data) && 0
    dex=find(coil_data.vert(4,:)==0);
    ndex=length(dex);
    dex0(2:ndex)=dex(1:ndex-1)+1;
    dex0(1)=1;
    subplot(4,3,[3 6]);
    hold on;
    for n=1:ndex
        n2=coil_data.vert(5,dex0(n));
        x=coil_data.vert(1,dex0(n):dex(n));
        y=coil_data.vert(2,dex0(n):dex(n));
        z=coil_data.vert(3,dex0(n):dex(n));
        phi=atan2(y,z);
        phi(phi<0) = phi(phi<0)+2*pi;
        r=sqrt(x.*x+y.*y);
        len=length(phi)-1;
        phi=phi(1:len);
        r=r(1:len);
        z=z(1:len);
        phi_spl=pchip(1:len,phi);
        d1=1;
        z0=[];
        z0(1)=fzero(@(x)ppval(phi_spl,x),d1,optimset('Display','off'));
        for j=2:len-1
            temp=fzero(@(x)ppval(phi_spl,x),j,optimset('Display','off'));
            if (temp >= z0(d1)+1 || temp <= z0(d1)-1)
                d1=d1+1;
                z0(d1)=temp;
            end
        end
        for j=1:length(z0)
            coil_points(1,j)=pchip(1:len,r,z0(j));
            coil_points(2,j)=pchip(1:len,z,z0(j));
            coil_points(3,j)=coil_data.vert(5,dex0(n));
        end
            
    end        
end

% Handle only plotting last equilibria
ifirst = 1;
if (llast)
    ifirst = nfiles;
end

% Plot Files
for i=ifirst:nfiles
    switch plottype
        case 'opt'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                toroslice(r,j,z,1:data.ns,'Color','white');
                hold on
                plot(r(1,1,j),z(1,1,j),'+y','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            plot(data.phi./data.phi(data.ns),data.presf,'w');
            xlim([0 1.1]);
            %if i == 1
                yvals=ylim;
                yvals(1) = 0;
                yvals(2) = max(data.presf)*1.1;
            %end
            if (yvals(2) > yvals(1)) 
                ylim(yvals);
            end
            ylabel('Pressure');
            % Plot jcurv and iotaf
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            plot(data.phi./data.phi(data.ns),data.jcurv);
            xlim([0 1.1]);
            ylabel('<J_V>');
            %plotyy(data.phi./data.phi(data.ns),data.jcurv,data.phi./data.phi(data.ns),data.iotaf);
        case 'te'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            switch axistype
                case 'topdown'
                    ltop_down=1;
                    subplot(4,3,[3 6]);
                    hold off;
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    theta=0:2*pi/90:2*pi;
                    zeta=0:2*pi/90:2*pi;
                    r=cfunct(theta,zeta,data.rmnc(:,[1 data.ns]),data.xm,data.xn);
                    z=sfunct(theta,zeta,data.zmns(:,[1 data.ns]),data.xm,data.xn);
                    if (data.iasym == 1)
                        r=r+sfunct(theta,zeta,data.rmns(:,[1 data.ns]),data.xm,data.xn);
                        z=z+cfunct(theta,zeta,data.zmnc(:,[1 data.ns]),data.xm,data.xn);
                    end
                    plot(squeeze(r(1,1,:)).*cos(zeta'),squeeze(r(1,1,:)).*sin(zeta'),'w--');
                    hold on;
                    for u=1:length(zeta)
                        f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn)+cfunct(x,zeta(u),data.zmnc(:,data.ns),data.xm,data.xn);
                        end
                        theta1=fzero(f,[-pi/2 pi/2]);
                        theta2=fzero(f,[pi/2 3*pi/2]);
                        r_temp=cfunct([theta1 theta2],zeta(u),data.rmnc(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            r_temp=r_temp+sfunct([theta1 theta2],zeta(u),data.rmns(:,data.ns),data.xm,data.xn);
                        end
                        r(2,1,u) = r_temp(1);
                        r(2,2,u) = r_temp(2);
                    end
                    plot(squeeze(r(2,1,:)).*cos(zeta'),squeeze(r(2,1,:)).*sin(zeta'),'r-');
                    plot(squeeze(r(2,2,:)).*cos(zeta'),squeeze(r(2,2,:)).*sin(zeta'),'r-');
                    axis equal;
                    subplot(4,3,9);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    errorbar(stel_data.TE_R(i,:),stel_data.TE_target(i,:),stel_data.TE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.TE_R(i,:),stel_data.TE_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TE_equil(i,:)-stel_data.TE_target(i,:))./stel_data.TE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.TE_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.TE_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.TE_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.TE_R(i,:).*cos(stel_data.TE_PHI(i,:)),...
                        stel_data.TE_R(i,:).*sin(stel_data.TE_PHI(i,:)),'og');
                    xlim([-1 1].*vmec_xmax);
                    ylim([-1 1].*vmec_xmax);
                    axis square;
                case 'realspace'
                    errorbar(stel_data.TE_R(i,:),stel_data.TE_target(i,:),stel_data.TE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.TE_R(i,:),stel_data.TE_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TE_equil(i,:)-stel_data.TE_target(i,:))./stel_data.TE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.TE_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.TE_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.TE_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,thom_plot);
                    hold on
                    plot(stel_data.TE_R(i,:),stel_data.TE_Z(i,:),'+g');
                case 'zspace'
                    errorbar(stel_data.TE_Z(i,:),stel_data.TE_target(i,:),stel_data.TE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.TE_Z(i,:),stel_data.TE_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TE_equil(i,:)-stel_data.TE_target(i,:))./stel_data.TE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.TE_Z(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.TE_Z(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.TE_Z(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,thom_plot);
                    hold on
                    plot(stel_data.TE_R(i,:),stel_data.TE_Z(i,:),'+g');
                case 'fluxspace'
                    errorbar(stel_data.TE_S(i,:),stel_data.TE_target(i,:),stel_data.TE_sigma(i,:),'ow')
                    xlim([0 1.2]);
                    hold on
                    plot(stel_data.TE_S(i,:),stel_data.TE_equil(i,:),'or');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TE_equil(i,:)-stel_data.TE_target(i,:))./stel_data.TE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    nr=size(stel_data.TE_S,2);
                    if (max(reddex) > 0)
                        bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 nr]);
                    ylim([-5 5]);
                    subplot(4,3,thom_plot);
                    hold on
                    plot(stel_data.TE_R(i,:),stel_data.TE_Z(i,:),'+g');
            end
        case 'ti'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            switch axistype
                case 'topdown'
                    ltop_down=1;
                    subplot(4,3,[3 6]);
                    hold off;
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    theta=0:2*pi/90:2*pi;
                    zeta=0:2*pi/90:2*pi;
                    r=cfunct(theta,zeta,data.rmnc(:,[1 data.ns]),data.xm,data.xn);
                    z=sfunct(theta,zeta,data.zmns(:,[1 data.ns]),data.xm,data.xn);
                    if (data.iasym == 1)
                        r=r+sfunct(theta,zeta,data.rmns(:,[1 data.ns]),data.xm,data.xn);
                        z=z+cfunct(theta,zeta,data.zmnc(:,[1 data.ns]),data.xm,data.xn);
                    end
                    plot(squeeze(r(1,1,:)).*cos(zeta'),squeeze(r(1,1,:)).*sin(zeta'),'w--');
                    hold on;
                    for u=1:length(zeta)
                        f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn)+cfunct(x,zeta(u),data.zmnc(:,data.ns),data.xm,data.xn);
                        end
                        theta1=fzero(f,[-pi/2 pi/2]);
                        theta2=fzero(f,[pi/2 3*pi/2]);
                        r_temp=cfunct([theta1 theta2],zeta(u),data.rmnc(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            r_temp=r_temp+sfunct([theta1 theta2],zeta(u),data.rmns(:,data.ns),data.xm,data.xn);
                        end
                        r(2,1,u) = r_temp(1);
                        r(2,2,u) = r_temp(2);
                    end
                    plot(squeeze(r(2,1,:)).*cos(zeta'),squeeze(r(2,1,:)).*sin(zeta'),'r-');
                    plot(squeeze(r(2,2,:)).*cos(zeta'),squeeze(r(2,2,:)).*sin(zeta'),'r-');
                    axis equal;
                    subplot(4,3,9);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    errorbar(stel_data.TI_R(i,:),stel_data.TI_target(i,:),stel_data.TI_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.TI_R(i,:),stel_data.TI_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TI_equil(i,:)-stel_data.TI_target(i,:))./stel_data.TI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.TI_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.TI_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.TI_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.TI_R(i,:).*cos(stel_data.TI_PHI(i,:)),...
                        stel_data.TI_R(i,:).*sin(stel_data.TI_PHI(i,:)),'og');
                    xlim([-1 1].*vmec_xmax);
                    ylim([-1 1].*vmec_xmax);
                    axis square;
                case 'realspace'
                    errorbar(stel_data.TI_R(i,:),stel_data.TI_target(i,:),stel_data.TI_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.TI_R(i,:),stel_data.TI_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TI_equil(i,:)-stel_data.TI_target(i,:))./stel_data.TI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.TI_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.TI_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.TI_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.TI_R(i,:),stel_data.TI_Z(i,:),'+g');
                case 'zspace'
                    errorbar(abs(stel_data.TI_Z(i,:)),stel_data.TI_target(i,:),stel_data.TI_sigma(i,:),'ow');
                    hold on
                    plot(abs(stel_data.TI_Z(i,:)),stel_data.TI_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TI_equil(i,:)-stel_data.TI_target(i,:))./stel_data.TI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(abs(stel_data.TI_Z(i,:)),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(abs(stel_data.TI_Z(i,:)),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(abs(stel_data.TI_Z(i,:)),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.TI_R(i,:),stel_data.TI_Z(i,:),'+g');
                case 'fluxspace'
                    errorbar(stel_data.TI_S(i,:),stel_data.TI_target(i,:),stel_data.TI_sigma(i,:),'ow')
                    xlim([0 1.2]);
                    hold on
                    plot(stel_data.TI_S(i,:),stel_data.TI_equil(i,:),'or');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.TI_equil(i,:)-stel_data.TI_target(i,:))./stel_data.TI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    nr=size(stel_data.TI_S,2);
                    if (max(reddex) > 0)
                        bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 nr]);
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.TI_R(i,:),stel_data.TI_Z(i,:),'+g');
            end
        case 'vphi'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            switch axistype
                case 'topdown'
                    ltop_down=1;
                    subplot(4,3,[3 6]);
                    hold off;
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    theta=0:2*pi/90:2*pi;
                    zeta=0:2*pi/90:2*pi;
                    r=cfunct(theta,zeta,data.rmnc(:,[1 data.ns]),data.xm,data.xn);
                    z=sfunct(theta,zeta,data.zmns(:,[1 data.ns]),data.xm,data.xn);
                    if (data.iasym == 1)
                        r=r+sfunct(theta,zeta,data.rmns(:,[1 data.ns]),data.xm,data.xn);
                        z=z+cfunct(theta,zeta,data.zmnc(:,[1 data.ns]),data.xm,data.xn);
                    end
                    plot(squeeze(r(1,1,:)).*cos(zeta'),squeeze(r(1,1,:)).*sin(zeta'),'w--');
                    hold on;
                    for u=1:length(zeta)
                        f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn)+cfunct(x,zeta(u),data.zmnc(:,data.ns),data.xm,data.xn);
                        end
                        theta1=fzero(f,[-pi/2 pi/2]);
                        theta2=fzero(f,[pi/2 3*pi/2]);
                        r_temp=cfunct([theta1 theta2],zeta(u),data.rmnc(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            r_temp=r_temp+sfunct([theta1 theta2],zeta(u),data.rmns(:,data.ns),data.xm,data.xn);
                        end
                        r(2,1,u) = r_temp(1);
                        r(2,2,u) = r_temp(2);
                    end
                    plot(squeeze(r(2,1,:)).*cos(zeta'),squeeze(r(2,1,:)).*sin(zeta'),'r-');
                    plot(squeeze(r(2,2,:)).*cos(zeta'),squeeze(r(2,2,:)).*sin(zeta'),'r-');
                    axis equal;
                    subplot(4,3,9);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    errorbar(stel_data.VPHI_R(i,:),stel_data.VPHI_target(i,:),stel_data.VPHI_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.VPHI_R(i,:),stel_data.VPHI_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.VPHI_equil(i,:)-stel_data.VPHI_target(i,:))./stel_data.VPHI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.VPHI_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.VPHI_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.VPHI_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.VPHI_R(i,:).*cos(stel_data.VPHI_PHI(i,:)),...
                        stel_data.VPHI_R(i,:).*sin(stel_data.VPHI_PHI(i,:)),'og');
                    xlim([-1 1].*vmec_xmax);
                    ylim([-1 1].*vmec_xmax);
                    axis square;
                case 'realspace'
                    errorbar(stel_data.VPHI_R(i,:),stel_data.VPHI_target(i,:),stel_data.VPHI_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.VPHI_R(i,:),stel_data.VPHI_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.VPHI_equil(i,:)-stel_data.VPHI_target(i,:))./stel_data.VPHI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.VPHI_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.VPHI_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.VPHI_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.VPHI_R(i,:),stel_data.VPHI_Z(i,:),'+g');
                case 'zspace'
                    errorbar(abs(stel_data.VPHI_Z(i,:)),stel_data.VPHI_target(i,:),stel_data.VPHI_sigma(i,:),'ow');
                    hold on
                    plot(abs(stel_data.VPHI_Z(i,:)),stel_data.VPHI_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.VPHI_equil(i,:)-stel_data.VPHI_target(i,:))./stel_data.VPHI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(abs(stel_data.VPHI_Z(i,:)),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(abs(stel_data.VPHI_Z(i,:)),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(abs(stel_data.VPHI_Z(i,:)),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.VPHI_R(i,:),stel_data.VPHI_Z(i,:),'+g');
                case 'fluxspace'
                    errorbar(stel_data.VPHI_S(i,:),stel_data.VPHI_target(i,:),stel_data.VPHI_sigma(i,:),'ow')
                    xlim([0 1.2]);
                    hold on
                    plot(stel_data.VPHI_S(i,:),stel_data.VPHI_equil(i,:),'or');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.VPHI_equil(i,:)-stel_data.VPHI_target(i,:))./stel_data.VPHI_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    nr=size(stel_data.VPHI_S,2);
                    if (max(reddex) > 0)
                        bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 nr]);
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.VPHI_R(i,:),stel_data.VPHI_Z(i,:),'+g');
            end
        case 'ne'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            switch axistype
                case 'topdown'
                    ltop_down=1;
                    subplot(4,3,[3 6]);
                    hold off;
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    theta=0:2*pi/90:2*pi;
                    zeta=0:2*pi/90:2*pi;
                    r=cfunct(theta,zeta,data.rmnc(:,[1 data.ns]),data.xm,data.xn);
                    z=sfunct(theta,zeta,data.zmns(:,[1 data.ns]),data.xm,data.xn);
                    if (data.iasym == 1)
                        r=r+sfunct(theta,zeta,data.rmns(:,[1 data.ns]),data.xm,data.xn);
                        z=z+cfunct(theta,zeta,data.zmnc(:,[1 data.ns]),data.xm,data.xn);
                    end
                    plot(squeeze(r(1,1,:)).*cos(zeta'),squeeze(r(1,1,:)).*sin(zeta'),'w--');
                    hold on;
                    for u=1:length(zeta)
                        f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn)+cfunct(x,zeta(u),data.zmnc(:,data.ns),data.xm,data.xn);
                        end
                        theta1=fzero(f,[-pi/2 pi/2]);
                        theta2=fzero(f,[pi/2 3*pi/2]);
                        r_temp=cfunct([theta1 theta2],zeta(u),data.rmnc(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            r_temp=r_temp+sfunct([theta1 theta2],zeta(u),data.rmns(:,data.ns),data.xm,data.xn);
                        end
                        r(2,1,u) = r_temp(1);
                        r(2,2,u) = r_temp(2);
                    end
                    plot(squeeze(r(2,1,:)).*cos(zeta'),squeeze(r(2,1,:)).*sin(zeta'),'r-');
                    plot(squeeze(r(2,2,:)).*cos(zeta'),squeeze(r(2,2,:)).*sin(zeta'),'r-');
                    axis equal;
                    subplot(4,3,9);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    errorbar(stel_data.NE_R(i,:),stel_data.NE_target(i,:),stel_data.NE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.NE_R(i,:),stel_data.NE_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.NE_equil(i,:)-stel_data.NE_target(i,:))./stel_data.NE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.NE_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.NE_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.NE_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(f_{sim}-f_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,thom_plot);
                    hold on
                    plot(stel_data.NE_R(i,:).*cos(stel_data.NE_PHI(i,:)),...
                        stel_data.NE_R(i,:).*sin(stel_data.NE_PHI(i,:)),'og');
                    xlim([-1 1].*vmec_xmax);
                    ylim([-1 1].*vmec_xmax);
                    axis square;
                case 'realspace'
                    errorbar(stel_data.NE_R(i,:),stel_data.NE_target(i,:),stel_data.NE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.NE_R(i,:),stel_data.NE_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.NE_equil(i,:)-stel_data.NE_target(i,:))./stel_data.NE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.NE_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.NE_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.NE_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(f_{sim}-f_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,thom_plot);
                    hold on
                    plot(stel_data.NE_R(i,:),stel_data.NE_Z(i,:),'+g');
                case 'zspace'
                    errorbar(stel_data.NE_Z(i,:),stel_data.NE_target(i,:),stel_data.NE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.NE_Z(i,:),stel_data.NE_equil(i,:),'or');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.NE_equil(i,:)-stel_data.NE_target(i,:))./stel_data.NE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.NE_Z(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.NE_Z(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.NE_Z(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(f_{sim}-f_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,thom_plot);
                    hold on
                    plot(stel_data.NE_R(i,:),stel_data.NE_Z(i,:),'+g');
                case 'fluxspace'
                    errorbar(stel_data.NE_S(i,:),stel_data.NE_target(i,:),stel_data.NE_sigma(i,:),'ow')
                    xlim([0 1.2]);
                    hold on
                    plot(stel_data.NE_S(i,:),stel_data.NE_equil(i,:),'or');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.NE_equil(i,:)-stel_data.NE_target(i,:))./stel_data.NE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    nr=size(stel_data.NE_S,2);
                    if (max(reddex) > 0)
                        bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(f_{sim}-f_{exp})/\sigma');
                    xlim([0 nr]);
                    ylim([-5 5]);
                    subplot(4,3,thom_plot);
                    hold on
                    plot(stel_data.NE_R(i,:),stel_data.NE_Z(i,:),'+g');
            end
        case 'neline'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            errorbar(stel_data.NELINE_target(i,:),stel_data.NELINE_sigma(i,:),'ow')
            hold on
            plot(stel_data.NELINE_equil(i,:),'or');
            hold off
            xlabel('Chord');
            if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.NELINE_equil(i,:)-stel_data.NELINE_target(i,:))./stel_data.NELINE_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.NELINE_equil,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(f_{sim}-f_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
            subplot(4,3,neline_plot);
            hold on
            for j=1:size(stel_data.NELINE_R0,2)
                if (abs(thom_wegt(j)) > 2.0)
                    temp_color='r';
                elseif(abs(thom_wegt(j))>1.0)&&(abs(thom_wegt(j))<=2.0)
                    temp_color='y';
                else
                    temp_color='g';
                end
                plot([stel_data.NELINE_R0(i,j) stel_data.NELINE_R1(i,j)],[stel_data.NELINE_Z0(i,j) stel_data.NELINE_Z1(i,j)],'Color',temp_color,'LineWidth',2.0);
            end
            hold off
        case 'teline'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            errorbar(stel_data.TELINE_target(i,:),stel_data.TELINE_sigma(i,:),'ow')
            hold on
            plot(stel_data.TELINE_equil(i,:),'or');
            hold off
            xlabel('Chord');
            if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.TELINE_equil(i,:)-stel_data.TELINE_target(i,:))./stel_data.TELINE_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.TELINE_equil,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(f_{sim}-f_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
            subplot(4,3,neline_plot);
            hold on
            for j=1:size(stel_data.TELINE_R0,2)
                if (abs(thom_wegt(j)) > 2.0)
                    temp_color='r';
                elseif(abs(thom_wegt(j))>1.0)&&(abs(thom_wegt(j))<=2.0)
                    temp_color='y';
                else
                    temp_color='g';
                end
                plot([stel_data.TELINE_R0(i,j) stel_data.TELINE_R1(i,j)],[stel_data.TELINE_Z0(i,j) stel_data.TELINE_Z1(i,j)],'Color',temp_color,'LineWidth',2.0);
            end
            hold off
        case 'tiline'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            errorbar(stel_data.TILINE_target(i,:),stel_data.TILINE_sigma(i,:),'ow')
            hold on
            plot(stel_data.TILINE_equil(i,:),'or');
            hold off
            xlabel('Chord');
            if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.TILINE_equil(i,:)-stel_data.TILINE_target(i,:))./stel_data.TILINE_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.TILINE_equil,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(f_{sim}-f_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
            subplot(4,3,neline_plot);
            hold on
            for j=1:size(stel_data.TILINE_R0,2)
                if (abs(thom_wegt(j)) > 2.0)
                    temp_color='r';
                elseif(abs(thom_wegt(j))>1.0)&&(abs(thom_wegt(j))<=2.0)
                    temp_color='y';
                else
                    temp_color='g';
                end
                plot([stel_data.TILINE_R0(i,j) stel_data.TILINE_R1(i,j)],[stel_data.TILINE_Z0(i,j) stel_data.TILINE_Z1(i,j)],'Color',temp_color,'LineWidth',2.0);
            end
            hold off
        case 'sxr'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            errorbar(stel_data.SXR_target(i,:),stel_data.SXR_sigma(i,:),'ow')
            hold on
            plot(stel_data.SXR_equil(i,:),'or');
            hold off
            xlabel('Chord');
            if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.SXR_equil(i,:)-stel_data.SXR_target(i,:))./stel_data.SXR_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.SXR_equil,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(f_{sim}-f_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
            subplot(4,3,neline_plot);
            hold on
            for j=1:size(stel_data.SXR_R0,2)
                if (abs(thom_wegt(j)) > 2.0)
                    temp_color='r';
                elseif(abs(thom_wegt(j))>1.0)&&(abs(thom_wegt(j))<=2.0)
                    temp_color='y';
                else
                    temp_color='g';
                end
                plot([stel_data.SXR_R0(i,j) stel_data.SXR_R1(i,j)],[stel_data.SXR_Z0(i,j) stel_data.SXR_Z1(i,j)],'Color',temp_color,'LineWidth',2.0);
            end
            hold off
        case 'faraday'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            errorbar(stel_data.FARADAY_target(i,:),stel_data.FARADAY_sigma(i,:),'ow')
            xlim([0 size(stel_data.FARADAY_target,2)+1]);
            hold on
            plot(stel_data.FARADAY_equil(i,:),'or');
            hold off
            xlabel('Chord');
            %if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.FARADAY_equil(i,:)-stel_data.FARADAY_target(i,:))./stel_data.FARADAY_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.FARADAY_equil,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(f_{sim}-f_{exp})/\sigma');
            xlim([0 size(stel_data.FARADAY_target,2)+1]);
            ylim([-5 5]);
            subplot(4,3,faraday_plot);
            hold on
            for j=1:size(stel_data.FARADAY_R0,2)
                if (abs(thom_wegt(j)) > 2.0)
                    temp_color='r';
                elseif(abs(thom_wegt(j))>1.0)&&(abs(thom_wegt(j))<=2.0)
                    temp_color='y';
                else
                    temp_color='g';
                end
                plot([stel_data.FARADAY_R0(i,j) stel_data.FARADAY_R1(i,j)],[stel_data.FARADAY_Z0(i,j) stel_data.FARADAY_Z1(i,j)],'Color',temp_color,'LineWidth',2.0);
            end
            hold off
        case 'mse'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D current
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('jet');
            end         
            % Plot Current profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            switch axistype
                case 'topdown'
                    ltop_down=1;
                    subplot(4,3,[3 6]);
                    hold off;
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    theta=0:2*pi/90:2*pi;
                    zeta=0:2*pi/90:2*pi;
                    r=cfunct(theta,zeta,data.rmnc(:,[1 data.ns]),data.xm,data.xn);
                    z=sfunct(theta,zeta,data.zmns(:,[1 data.ns]),data.xm,data.xn);
                    if (data.iasym == 1)
                        r=r+sfunct(theta,zeta,data.rmns(:,[1 data.ns]),data.xm,data.xn);
                        z=z+cfunct(theta,zeta,data.zmnc(:,[1 data.ns]),data.xm,data.xn);
                    end
                    plot(squeeze(r(1,1,:)).*cos(zeta'),squeeze(r(1,1,:)).*sin(zeta'),'w--');
                    hold on;
                    for u=1:length(zeta)
                        f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            f = @(x) sfunct(x,zeta(u),data.zmns(:,data.ns),data.xm,data.xn)+cfunct(x,zeta(u),data.zmnc(:,data.ns),data.xm,data.xn);
                        end
                        theta1=fzero(f,[-pi/2 pi/2]);
                        theta2=fzero(f,[pi/2 3*pi/2]);
                        r_temp=cfunct([theta1 theta2],zeta(u),data.rmnc(:,data.ns),data.xm,data.xn);
                        if (data.iasym == 1)
                            r_temp=r_temp+sfunct([theta1 theta2],zeta(u),data.rmns(:,data.ns),data.xm,data.xn);
                        end
                        r(2,1,u) = r_temp(1);
                        r(2,2,u) = r_temp(2);
                    end
                    plot(squeeze(r(2,1,:)).*cos(zeta'),squeeze(r(2,1,:)).*sin(zeta'),'r-');
                    plot(squeeze(r(2,2,:)).*cos(zeta'),squeeze(r(2,2,:)).*sin(zeta'),'r-');
                    axis equal;
                    subplot(4,3,9);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    errorbar(stel_data.MSE_R(i,:),stel_data.MSE_target(i,:),stel_data.MSE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.MSE_R(i,:),stel_data.MSE_equil(i,:),'or');
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.MSE_equil(i,:)-stel_data.MSE_target(i,:))./stel_data.MSE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.MSE_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.MSE_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.MSE_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(f_{sim}-f_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.MSE_R(i,:).*cos(stel_data.MSE_PHI(i,:)),...
                        stel_data.MSE_R(i,:).*sin(stel_data.MSE_PHI(i,:)),'og');
                    xlim([-1 1].*vmec_xmax);
                    ylim([-1 1].*vmec_xmax);
                    axis square;
                case 'realspace'
                    errorbar(stel_data.MSE_R(i,:),stel_data.MSE_target(i,:),stel_data.MSE_sigma(i,:),'ow');
                    hold on
                    plot(stel_data.MSE_R(i,:),stel_data.MSE_equil(i,:),'or');
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.MSE_equil(i,:)-stel_data.MSE_target(i,:))./stel_data.MSE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(stel_data.MSE_R(i,:),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(stel_data.MSE_R(i,:),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(stel_data.MSE_R(i,:),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(f_{sim}-f_{exp})/\sigma');
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.MSE_R(i,:),stel_data.MSE_Z(i,:),'+g');
                case 'fluxspace'
                    fac=180/pi;
                    errorbar(stel_data.MSE_S(i,:),fac*stel_data.MSE_target(i,:),fac*stel_data.MSE_sigma(i,:),'ow')
                    xlim([0 1.2]);
                    %yrange=ylim;
                    hold on
                    plot(stel_data.MSE_S(i,:),fac*stel_data.MSE_equil(i,:),'or');
                    hold off
                    %ylim(yrange);
                    xlabel('Normalized Flux');
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                    thom_wegt=(stel_data.MSE_equil(i,:)-stel_data.MSE_target(i,:))./stel_data.MSE_sigma(i,:);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    nr=size(stel_data.MSE_S,2);
                    if (max(reddex) > 0)
                        bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 nr]);
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(stel_data.MSE_R(i,:),stel_data.MSE_Z(i,:),'+g');
            end
        case 'iota'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D current
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Iota profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            errorbar(stel_data.IOTA_S(i,:),stel_data.IOTA_target(i,:),stel_data.IOTA_sigma(i,:),'ow')
            xlim([0 1.2]);
            yrange=ylim;
            hold on
            plot(stel_data.IOTA_S(i,:),stel_data.IOTA_equil(i,:),'or');
            hold off
            ylim(yrange);
            xlabel('Normalized Flux');
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.IOTA_equil(i,:)-stel_data.IOTA_target(i,:))./stel_data.IOTA_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.IOTA_S,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(\iota_{sim}-\iota_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
            subplot(4,3,[3 6]);
            hold on
        case 'mag'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Plot VMEC data
            for j=1:2
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                toroslice(r,j,z,1:data.ns-1,'Color','w');
                hold on
                plot(r(1,1,j),z(1,1,j),'+w','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
            end
            subplot(4,3,[3 6]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            cla;  % So isotoro doesn't overplot itself
            theta=0:2*pi/71:2*pi;
            zeta=0:2*pi/71:2*pi;
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            hplas=isotoro(r,z,zeta,data.ns);
            camlight left;
            set(hplas,'FaceAlpha',0.5,'FaceColor','red');
            set(gcf,'Renderer','zbuffer');
            hold off
            % Plot Magnetic data
            subplot(4,3,9); cla;
            subplot(4,3,12); cla;
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            nr=0;
            if isfield(stel_data,'B_PROBES_target')
                nr=size(stel_data.B_PROBES_target,2);
                errorbar(1:nr,stel_data.B_PROBES_target(i,:),stel_data.B_PROBES_sigma(i,:),'ow')
                hold on
                plot(1:nr,stel_data.B_PROBES_equil(i,:),'or');
                hold off
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                ylim(ylim_mag);
                xlabel('Index');
                subplot(4,3,12);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                thom_wegt=(stel_data.B_PROBES_equil(i,:)-stel_data.B_PROBES_target(i,:))./stel_data.B_PROBES_sigma(i,:);
                greydex = abs(stel_data.B_PROBES_sigma(i,:)) >= 1.0E10;
                reddex=(abs(thom_wegt)>2.0).*(abs(stel_data.B_PROBES_sigma(i,:)) < 1.0E10);
                yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0).*(abs(stel_data.B_PROBES_sigma(i,:)) < 1.0E10);
                gredex=(abs(thom_wegt)<=1.0).*(abs(stel_data.B_PROBES_sigma(i,:)) < 1.0E10);
                if (max(reddex) > 0)
                    bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                    hold on
                end
                if (max(yeldex) > 0)
                    bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                    hold on
                end
                if (max(gredex) > 0)
                    bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                    hold off
                end
                hold off
                title('Weighted Deviation','Color','white')
                xlabel('Index');
                ylabel('(f_{sim}-f_{exp})/\sigma');
                xlim([0 nr]);
                ylim([-5 5]);
                subplot(4,3,[3 6]);
                hold on
                for n=1:nr
                    color_temp=reddex(n)*[1 0 0]+yeldex(n)*[1 1 0] + gredex(n)*[0 1 0] + greydex(n)*[1 1 1].*grey_val;
                    plot3(stel_data.B_PROBES_X(i,n),stel_data.B_PROBES_Y(i,n),stel_data.B_PROBES_Z(i,n),'+','Color',color_temp,'MarkerSize',18,'LineWidth',2);
                end
                hold off
                if lrz_mag
                    subplot(4,3,[1 4]);
                    hold on
                    for n=1:nr
                        color_temp=reddex(n)*[1 0 0]+yeldex(n)*[1 1 0] + gredex(n)*[0 1 0] + greydex(n)*[1 1 1].*grey_val;
                        r_temp=sqrt(stel_data.B_PROBES_X(i,n).*stel_data.B_PROBES_X(i,n)+stel_data.B_PROBES_Y(i,n).*stel_data.B_PROBES_Y(i,n));
                        z_temp=stel_data.B_PROBES_Z(i,n);
                        plot(r_temp,z_temp,'+','Color',color_temp,'MarkerSize',14,'LineWidth',2);
                    end
                    hold off
                end
            end
            if isfield(stel_data,'FLUXLOOPS_target')
                subplot(4,3,9);
                nr2=size(stel_data.FLUXLOOPS_target,2)+nr;
                dex=size(stel_data.FLUXLOOPS_target,2);
                hold on
                errorbar(nr+1:nr2,stel_data.FLUXLOOPS_target(i,:),stel_data.FLUXLOOPS_sigma(i,:),'ow')
                plot(nr+1:nr2,stel_data.FLUXLOOPS_equil(i,:),'or');
                ylim(ylim_mag);
                xlim([0 nr2]);
                hold off
                xlabel('Index');
                subplot(4,3,12);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                thom_wegt=(stel_data.FLUXLOOPS_equil(i,:)-stel_data.FLUXLOOPS_target(i,:))./stel_data.FLUXLOOPS_sigma(i,:);
                greydex = abs(stel_data.FLUXLOOPS_sigma(i,:)) >= 1.0E10;
                reddex=(abs(thom_wegt)>2.0).*(abs(stel_data.FLUXLOOPS_sigma(i,:)) < 1.0E10);
                yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0).*(abs(stel_data.FLUXLOOPS_sigma(i,:)) < 1.0E10);
                gredex=(abs(thom_wegt)<=1.0).*(abs(stel_data.FLUXLOOPS_sigma(i,:)) < 1.0E10);
                hold on;
                if (max(reddex) > 0)
                    bar(nr+1:nr2,thom_wegt.*reddex,'r','EdgeColor','red');
                    hold on
                end
                if (max(yeldex) > 0)
                    bar(nr+1:nr2,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                    hold on
                end
                if (max(gredex) > 0)
                    bar(nr+1:nr2,thom_wegt.*gredex,'g','EdgeColor','green');
                    hold off
                end
                hold off; hold off; hold off;
                title('Weighted Deviation','Color','white')
                xlabel('Index');
                ylabel('(p_{sim}-p_{exp})/\sigma');
                xlim([0 nr2]);
                ylim([-5 5]);
                if ~isempty(mag_specs)
                    subplot(4,3,[3 6]);
                    hold on
                    if isfield(mag_specs,'fluxloops')
                        for n=1:size(stel_data.FLUXLOOPS_target,2)
                            color_temp=reddex(n)*[1 0 0]+yeldex(n)*[1 1 0] + gredex(n)*[0 1 0] + greydex(n)*[1 1 1].*grey_val;
                            plot3(mag_specs.fluxloops.loops{n}(1,:),...
                                mag_specs.fluxloops.loops{n}(2,:),...
                                mag_specs.fluxloops.loops{n}(3,:),'Color',color_temp);
                        end
                    end
                    hold off
                    if lrz_mag
                        subplot(4,3,[1 4]);
                        hold on
                        for n=1:size(stel_data.FLUXLOOPS_target,2)
                            color_temp=reddex(n)*[1 0 0]+yeldex(n)*[1 1 0] + gredex(n)*[0 1 0] + greydex(n)*[1 1 1].*grey_val;
                            r_temp=sqrt(mag_specs.fluxloops.loops{n}(1,1).*mag_specs.fluxloops.loops{n}(1,1)+mag_specs.fluxloops.loops{n}(2,1).*mag_specs.fluxloops.loops{n}(2,1));
                            z_temp=mag_specs.fluxloops.loops{n}(3,1);
                            line_temp=2.0;
                            marker_temp='o';
                            if color_temp(1) == grey_val, line_temp=0.5; marker_temp='x'; end
                            plot(r_temp,z_temp,marker_temp,'Color',color_temp,'MarkerSize',14,'LineWidth',line_temp);
                        end
                        hold off
                    end
                end
            end
        case 'coils'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                toroslice(r,j,z,1:data.ns-1,'Color','w');
                hold on
                plot(r(1,1,j),z(1,1,j),'+w','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
            end
            % Fix stel_data.EXTCUR_equil
            EXTCUR_equil = data.extcur;
            dex=size(stel_data.EXTCUR_target,2);
            EXTCUR_target = EXTCUR_equil.*0.0;
            EXTCUR_sigma = EXTCUR_equil.*0.0;
            EXTCUR_target(1:dex) = stel_data.EXTCUR_target(i,1:dex);
            EXTCUR_sigma(1:dex) = stel_data.EXTCUR_sigma(i,1:dex);
            EXTCUR_sigma(EXTCUR_sigma == 0) = 1.0;
            % Plot the external currents
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            dex=length(EXTCUR_target);
            errorbar(1:dex,EXTCUR_target,EXTCUR_sigma,'ow')
            hold on
            plot(1:dex,EXTCUR_equil,'or');
            hold off
            xlabel('Index');
            xlim([1 dex]);
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(EXTCUR_equil-EXTCUR_target)./EXTCUR_sigma;
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr = length(EXTCUR_target);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(p_{sim}-p_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
            % Handle 3D plot of diagnostics
            if ~isempty(coil_data)
                subplot(4,3,[3 6]);
                cla;  % So isotoro doesn't overplot itself
                theta=0:2*pi/71:2*pi;
                zeta=0:2*pi/71:2*pi;
                r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
                z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
                if (data.iasym == 1)
                    r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                    z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
                end
                hplas=isotoro(r,z,zeta,data.ns);
                camlight left;
                set(hplas,'FaceAlpha',0.5,'FaceColor','red');
                dex=find(coil_data.vert(4,:)==0);
                ndex=length(dex);
                dex0(2:ndex)=dex(1:ndex-1)+1;
                dex0(1)=1;
                subplot(4,3,[3 6]);
                hold on;
                for n=1:ndex
                    n2=coil_data.vert(5,dex0(n));
                    x=coil_data.vert(1,dex0(n):dex(n));
                    y=coil_data.vert(2,dex0(n):dex(n));
                    z=coil_data.vert(3,dex0(n):dex(n));
                    color_temp=reddex(n2)*[1 0 0]+yeldex(n2)*[1 1 0] + gredex(n2)*[0 1 0];
                    plot3(x,y,z,'Color',color_temp');
                end
                xlim([-vmec_xmax vmec_xmax]);
                ylim([-vmec_xmax vmec_xmax]);
                zlim([vmec_ymin vmec_ymax]);
                xlabel('X [m]');
                ylabel('Y [m]');
                zlabel('Z [m]');
                hold off;
            end
        case 'bnorm'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                toroslice(r,j,z,1:data.ns,'Color','white');
                hold on
                plot(r(1,1,j),z(1,1,j),'+y','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot BNORM profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            [xq,yq] = meshgrid(0:1/128:1,0:1/128:1);
            [BNEQ]=griddata(stel_data.COIL_BNORM_U,stel_data.COIL_BNORM_V,stel_data.COIL_BNORM_BNEQ,xq,yq);
            [BNF]=griddata(stel_data.COIL_BNORM_U,stel_data.COIL_BNORM_V,stel_data.COIL_BNORM_BNF,xq,yq);
            pixplot(0:1/128:1,0:1/128:1,BNF-BNEQ);
            set(gca,'XTick',[],'YTick',[0 0.5 0.98],'YTickLabels',{'0' '0.5' '1'});
            ylabel('V');
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            pixplot(0:1/128:1,0:1/128:1,log(abs(BNF-BNEQ)./abs(BNEQ)));
            xlabel('U');
            set(gca,'XTick',[0 0.25 0.5 0.75 0.98],'XTickLabels',{'0' '0.25' '0.5' '0.75' '1'},'YTick',[0 0.5 0.9],'YTickLabels',{'0' '0.5' '1'});
        case 'bootstrap'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D current
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('jet');
            end
            % Plot Current profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            plot(stel_data.BOOTSTRAP_s(i,:),stel_data.BOOTSTRAP_boot_jdotb(i,:),'or');
            hold on
            %plot(stel_data.BOOTSTRAP_s(i,:),stel_data.BOOTSTRAP_avg_jdotb(i,:),'b');
            plot(stel_data.BOOTSTRAP_s(i,:),stel_data.BOOTSTRAP_beam_jdotb(i,:),'+c');
            plot(data.phi./data.phi(data.ns),data.jdotb,'b');
            %plot(stel_data.BOOTSTRAP_s(i,:),stel_data.BOOTSTRAP_jBbs(i,:).*stel_data.BOOTSTRAP_facnu(i,:),'w');
            errorbar(stel_data.BOOTSTRAP_s(i,:),stel_data.BOOTSTRAP_jBbs(i,:).*stel_data.BOOTSTRAP_facnu(i,:),stel_data.BOOTSTRAP_sigma(i,:),'ow');
            hold off
            xlabel('Normalized Flux');
            xlim([0 1.1]);
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.BOOTSTRAP_equil(i,:)-stel_data.BOOTSTRAP_target(i,:))./stel_data.BOOTSTRAP_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.BOOTSTRAP_s,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(p_{sim}-p_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
        case 'helicity'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D current
            if (isfield(data,'protmnc'))
                p=cfunct(theta,zeta,data.protmnc,data.xm,data.xn);
                if (data.iasym == 1)
                    p=p+sfunct(theta,zeta,data.protmns,data.xm,data.xn);
                end
            else
                p=zeros(data.ns,361,3);
                for j=1:data.ns, p(j,:,:)=data.presf(j); end
            end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                toroslice(r,j,z,1:data.ns-1,'Color','w');
                hold on
                plot(r(1,1,j),z(1,1,j),'+w','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('jet');
            end
            % Create Helicity Map
            m1=min(stel_data.HELICITY_m(i,:));
            m2=max(stel_data.HELICITY_m(i,:));
            n1=min(stel_data.HELICITY_n(i,:));
            n2=max(stel_data.HELICITY_n(i,:));
            l1=size(stel_data.HELICITY_m,2);
            dm=m2-m1+1;
            dn=n2-n1+1;
            %dl=l1/(dm*dn);
            helicity=zeros(dm,dn);
            dexm=stel_data.HELICITY_m(1,:) == stel_data.HELICITY_m(1,1);
            dexn=stel_data.HELICITY_n(1,:) == stel_data.HELICITY_n(1,1);
            dl=find(dexn.*dexm);
            dl=dl(2)-1;
            for k=1:dl
                m=stel_data.HELICITY_m(i,k)-m1+1;
                n=stel_data.HELICITY_n(i,k)-n1+1;
                helicity(m,n)=max(stel_data.HELICITY_equil(i,k:dl:l1));
            end
            % Plot Current profile
            subplot(4,3,[9 12]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            pixplot(m1:m2,n1:n2,helicity);
            xlabel('m');
            ylabel('n');
            colorbar;
            if (i==1)
                hel_clim=caxis;
            end
            caxis(hel_clim);
        case 'neo'
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filelist(i).name);
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D current
            p=[];
            p=zeros(data.ns,361,3);
            %for j=1:size(stel_data.NEO_k,2)
            %    p(stel_data.NEO_k(i,j),:,:)=stel_data.NEO_equil(i,j);
            %end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',markersiz,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','center','VerticalAlignment','bottom');
                if data.iasym
                   plot(vmec_xmin+dy,z(1,1,j),'<','Color','yellow','MarkerSize',arrowsiz,'MarkerFaceColor','yellow');
                   text(vmec_xmin+2*dy,z(1,1,j),num2str(z(1,1,j),'%4.3f'),'Color','white','FontSize',fontsiz1,'HorizontalAlignment','left','VerticalAlignment','BaseLine');
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filelist(i).name,'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',fontsiz1+4,'VerticalAlignment','top'); end
                colormap('jet');
            end
            % Plot Current profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            %Find VMEC DOMAIN
            plot(stel_data.NEO_k(i,:),stel_data.NEO_equil(i,:),'w');
            xlim([0 data.ns]);
            xlabel('Normalized Flux');
            %Plot Pressure Residuals
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
            thom_wegt=(stel_data.NEO_equil(i,:)-stel_data.NEO_target(i,:))./stel_data.NEO_sigma(i,:);
            reddex=abs(thom_wegt)>2.0;
            yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
            gredex=abs(thom_wegt)<=1.0;
            nr=size(stel_data.NEO_k,2);
            if (max(reddex) > 0)
                bar(1:nr,thom_wegt.*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(1:nr,thom_wegt.*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(1:nr,thom_wegt.*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Index');
            ylabel('(p_{sim}-p_{exp})/\sigma');
            xlim([0 nr]);
            ylim([-5 5]);
    end
% Handle plotting the separatrix if present
    if isfield(stel_data,'SEPARATRIX_chisq')
        if (ltop_down) 
            subplot(4,3,thom_plot);
            ctemp = colormap;
            clen  = size(ctemp,1)-1;
            r_temp=stel_data.SEPARATRIX_R(i,:);
            phi_temp=stel_data.SEPARATRIX_PHI(i,:);
            z_temp=stel_data.SEPARATRIX_Z(i,:);
            chi_temp=stel_data.SEPARATRIX_chisq(i,:);
            chi_max=max(chi_temp(1,:));
            hold on
            for k=1:length(r_temp)
                itemp=round(clen.*chi_temp(1,k)./chi_max)+1;
                plot(r_temp(1,k).*cos(phi_temp(1,k)),r_temp(1,k).*sin(phi_temp(1,k)),'.','Color',ctemp(itemp,:),'MarkerSize',markersiz);
            end
            hold off
        elseif (size(stel_data.SEPARATRIX_R,2) < 5)
            subplot(4,3,thom_plot);
            ctemp = colormap;
            clen  = size(ctemp,1)-1;
            r_temp=stel_data.SEPARATRIX_R(i,:);
            z_temp=stel_data.SEPARATRIX_Z(i,:);
            chi_temp=stel_data.SEPARATRIX_chisq(i,:);
            chi_max=max(chi_temp(1,:));
            hold on
            for k=1:length(r_temp)
                itemp=round(clen.*chi_temp(1,k)./chi_max)+1;
                plot(r_temp(1,k),z_temp(1,k),'.','Color',ctemp(itemp,:),'MarkerSize',markersiz);
            end
            hold off
        else
            subplot(4,3,[1 4]);
            ctemp = colormap('Cool');
            clen  = size(ctemp,1)-1;
            dex = (stel_data.SEPARATRIX_PHI(i,:) == 0.0);
            r_temp=stel_data.SEPARATRIX_R(i,dex);
            z_temp=stel_data.SEPARATRIX_Z(i,dex);
            chi_temp=stel_data.SEPARATRIX_chisq(i,dex);
            chi_max=max(chi_temp(1,:));
            hold on
            for k=1:length(r_temp)-1
                itemp=round(clen.*chi_temp(1,k)./chi_max)+1;
                plot(r_temp(1,k),z_temp(1,k),'.','Color',ctemp(itemp,:),'MarkerSize',markersiz);
            end
            hold off
            subplot(4,3,[2 5]);
            for j=1:size(stel_data.SEPARATRIX_PHI,2);
                if stel_data.SEPARATRIX_PHI(i,j) > 2*pi/data.nfp/4
                    break;
                end
            end
            dex = (stel_data.SEPARATRIX_PHI(i,:) == stel_data.SEPARATRIX_PHI(i,j));
            if j>1
                dex2 = (stel_data.SEPARATRIX_PHI(i,:) == stel_data.SEPARATRIX_PHI(i,j-1));
            else
                dex2 = dex;
            end
            if ~isempty(dex)
                r_temp=0.5*(stel_data.SEPARATRIX_R(i,dex)+stel_data.SEPARATRIX_R(i,dex2));
                z_temp=0.5*(stel_data.SEPARATRIX_Z(i,dex)+stel_data.SEPARATRIX_Z(i,dex2));
                chi_temp=0.5*(stel_data.SEPARATRIX_chisq(i,dex)+stel_data.SEPARATRIX_chisq(i,dex2));
                chi_max=max(0.5*(stel_data.SEPARATRIX_chisq(i,dex)+stel_data.SEPARATRIX_chisq(i,dex2)));
                hold on
                for k=1:length(r_temp)-1
                    itemp=round(clen.*chi_temp(1,k)./chi_max)+1;
                    plot(r_temp(1,k),z_temp(1,k),'.','Color',ctemp(itemp,:),'MarkerSize',markersiz);
                end
                hold off
            end
            subplot(4,3,[3 6]);
            for j=1:size(stel_data.SEPARATRIX_PHI,2);
                if stel_data.SEPARATRIX_PHI(i,j) > 2*pi/data.nfp/2
                    break;
                end
            end
            dex = (stel_data.SEPARATRIX_PHI(i,:) == stel_data.SEPARATRIX_PHI(i,j));
            if j>1
                dex2 = (stel_data.SEPARATRIX_PHI(i,:) == stel_data.SEPARATRIX_PHI(i,j-1));
            else
                dex2 = dex;
            end
            if ~isempty(dex)
                r_temp=0.5*(stel_data.SEPARATRIX_R(i,dex)+stel_data.SEPARATRIX_R(i,dex2));
                z_temp=0.5*(stel_data.SEPARATRIX_Z(i,dex)+stel_data.SEPARATRIX_Z(i,dex2));
                chi_temp=0.5*(stel_data.SEPARATRIX_chisq(i,dex)+stel_data.SEPARATRIX_chisq(i,dex2));
                chi_max=max(0.5*(stel_data.SEPARATRIX_chisq(i,dex)+stel_data.SEPARATRIX_chisq(i,dex2)));
                hold on
                for k=1:length(r_temp)-1
                    itemp=round(clen.*chi_temp(1,k)./chi_max)+1;
                    plot(r_temp(1,k),z_temp(1,k),'.','Color',ctemp(itemp,:),'MarkerSize',markersiz);
                end
                hold off
            end
        end
    end
    % Plot Output
    subplot(4,3,[7 8 10 11]);
    cla;
    axis normal;
    set(gca,'YLimMode','auto');
    set(gca,'nextplot','replacechildren','XColor','white',...
        'YColor','white','ZColor','white','Color','black','FontSize',fontsiz1);
    set(gca,'Yscale','log');
    chisq_total = sum(((stel_data.TARGETS-stel_data.VALS)./stel_data.SIGMAS).^2,2);
    plot(stel_data.iter,chisq_total,'o','Color','white','MarkerSize',10,'LineWidth',2.0);
    legtext{1}='Total';
    note_string={['\chi^2_{Total}=' num2str(chisq_total(i),'%7.4e')]};
    if isfield(stel_data,'ASPECT_chisq')
        note_string=[note_string; ['\chi^2_{Aspect}=' num2str(stel_data.ASPECT_chisq(i),'%7.4e')]];
        legtext=[legtext 'Aspect Ratio'];
        hold on
        plot(stel_data.iter,stel_data.ASPECT_chisq,'.r','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'ASPECT_MAX_chisq')
        note_string=[note_string; ['\chi^2_{Aspect(max)}=' num2str(stel_data.ASPECT_MAX_chisq(i),'%7.4e')]];
        legtext=[legtext 'Aspect Ratio (max)'];
        hold on
        plot(stel_data.iter,stel_data.ASPECT_MAX_chisq,'.g','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'BETA_chisq')
        note_string=[note_string; ['\chi^2_{Beta}=' num2str(stel_data.BETA_chisq(i),'%7.4e')]];
        legtext=[legtext 'Beta'];
        hold on
        plot(stel_data.iter,stel_data.BETA_chisq,'.g','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'BETAPOL_chisq')
        note_string=[note_string; ['\chi^2_{Beta-pol}=' num2str(stel_data.BETAPOL_chisq(i),'%7.4e')]];
        legtext=[legtext 'Beta (pol)'];
        hold on
        plot(stel_data.iter,stel_data.BETAPOL_chisq,'og','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'BETATOR_chisq')
        note_string=[note_string; ['\chi^2_{Beta-tor}=' num2str(stel_data.BETATOR_chisq(i),'%7.4e')]];
        legtext=[legtext 'Beta (tor)'];
        hold on
        plot(stel_data.iter,stel_data.BETATOR_chisq,'sg','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'CURTOR_chisq')
        note_string=[note_string; ['\chi^2_{I-tor}=' num2str(stel_data.CURTOR_chisq(i),'%7.4e')]];
        legtext=[legtext 'Toroidal Current'];
        hold on
        plot(stel_data.iter,stel_data.CURTOR_chisq,'.c','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'PHIEDGE_chisq')
        note_string=[note_string; ['\chi^2_{\Phi}=' num2str(stel_data.PHIEDGE_chisq(i),'%7.4e')]];
        legtext=[legtext 'Phiedge'];
        hold on
        plot(stel_data.iter,stel_data.PHIEDGE_chisq,'.m','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'VOLUME_chisq')
        note_string=[note_string; ['\chi^2_{Volume}=' num2str(stel_data.VOLUME_chisq(i),'%7.4e')]];
        legtext=[legtext 'Volume'];
        hold on
        plot(stel_data.iter,stel_data.VOLUME_chisq,'.y','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'WP_chisq')
        note_string=[note_string; ['\chi^2_{Energy}=' num2str(stel_data.WP_chisq(i),'%7.4e')]];
        legtext=[legtext 'Stored Energy'];
        hold on
        plot(stel_data.iter,stel_data.WP_chisq,'.b','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'RBTOR_chisq')
        note_string=[note_string; ['\chi^2_{RBtor}=' num2str(stel_data.RBTOR_chisq(i),'%7.4e')]];
        legtext=[legtext 'R*Btor'];
        hold on
        plot(stel_data.iter,stel_data.RBTOR_chisq,'.w','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'R0_chisq')
        note_string=[note_string; ['\chi^2_{R0}=' num2str(stel_data.R0_chisq(i),'%7.4e')]];
        legtext=[legtext 'R0'];
        hold on
        plot(stel_data.iter,stel_data.R0_chisq,'^w','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'Z0_chisq')
        note_string=[note_string; ['\chi^2_{Z0}=' num2str(stel_data.Z0_chisq(i),'%7.4e')]];
        legtext=[legtext 'Z0'];
        hold on
        plot(stel_data.iter,stel_data.Z0_chisq,'<w','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'NE_chisq')
        note_string=[note_string; ['\chi^2_{Ne}=' num2str(sum(stel_data.NE_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Elec. Den.'];
        hold on
        plot(stel_data.iter,sum(stel_data.NE_chisq,2),'+r','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'IOTA_chisq')
        note_string=[note_string; ['\chi^2_{iota}=' num2str(sum(stel_data.IOTA_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Iota'];
        hold on
        plot(stel_data.iter,sum(stel_data.IOTA_chisq,2),'+r','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'PRESS_chisq')
        note_string=[note_string; ['\chi^2_{press}=' num2str(sum(stel_data.PRESS_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Pressure'];
        hold on
        plot(stel_data.iter,sum(stel_data.PRESS_chisq,2),'+g','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'JDOTB_chisq')
        note_string=[note_string; ['\chi^2_{J.B}=' num2str(sum(stel_data.JDOTB_chisq(i,:)),'%7.4e')]];
        legtext=[legtext '<J.B>'];
        hold on
        plot(stel_data.iter,sum(stel_data.JDOTB_chisq,2),'sr','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'JTOR_chisq')
        note_string=[note_string; ['\chi^2_{JCURV}=' num2str(sum(stel_data.JTOR_chisq(i,:)),'%7.4e')]];
        legtext=[legtext '<JCURV>'];
        hold on
        plot(stel_data.iter,sum(stel_data.JTOR_chisq,2),'sr','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'NELINE_chisq')
        note_string=[note_string; ['\chi^2_{Ne-line}=' num2str(sum(stel_data.NELINE_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Elec. Den. (line)'];
        hold on
        plot(stel_data.iter,sum(stel_data.NELINE_chisq,2),'sr','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'TELINE_chisq')
        note_string=[note_string; ['\chi^2_{Te-line}=' num2str(sum(stel_data.TELINE_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Elec. Temp. (line)'];
        hold on
        plot(stel_data.iter,sum(stel_data.TELINE_chisq,2),'sg','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'TILINE_chisq')
        note_string=[note_string; ['\chi^2_{Ti-line}=' num2str(sum(stel_data.TILINE_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Ion Temp. (line)'];
        hold on
        plot(stel_data.iter,sum(stel_data.TILINE_chisq,2),'sb','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'FARADAY_chisq')
        note_string=[note_string; ['\chi^2_{Faraday}=' num2str(sum(stel_data.FARADAY_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Faraday Rot.'];
        hold on
        plot(stel_data.iter,sum(stel_data.FARADAY_chisq,2),'sc','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'SXR_chisq')
        note_string=[note_string; ['\chi^2_{Soft X-Ray}=' num2str(sum(stel_data.SXR_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Soft X-Ray'];
        hold on
        plot(stel_data.iter,sum(stel_data.SXR_chisq,2),'sm','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'TE_chisq')
        note_string=[note_string; ['\chi^2_{Te}=' num2str(sum(stel_data.TE_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Elec. Temp.'];
        hold on
        plot(stel_data.iter,sum(stel_data.TE_chisq,2),'+g','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'TI_chisq')
        note_string=[note_string; ['\chi^2_{Ti}=' num2str(sum(stel_data.TI_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Ion Temp.'];
        hold on
        plot(stel_data.iter,sum(stel_data.TI_chisq,2),'+b','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'VPHI_chisq')
        note_string=[note_string; ['\chi^2_{VPHI}=' num2str(sum(stel_data.VPHI_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'Ion Vel.'];
        hold on
        plot(stel_data.iter,sum(stel_data.VPHI_chisq,2),'+c','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'MSE_chisq')
        note_string=[note_string; ['\chi^2_{MSE}=' num2str(sum(stel_data.MSE_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'MSE'];
        hold on
        plot(stel_data.iter,sum(stel_data.MSE_chisq,2),'or','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'B_PROBES_chisq')
        note_string=[note_string; ['\chi^2_{BPROBES}=' num2str(sum(stel_data.B_PROBES_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'BPROBES'];
        hold on
        plot(stel_data.iter,sum(stel_data.B_PROBES_chisq,2),'xr','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'FLUXLOOPS_chisq')
        note_string=[note_string; ['\chi^2_{FLUXLOOPS}=' num2str(sum(stel_data.FLUXLOOPS_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'FLUXLOOPS'];
        hold on
        plot(stel_data.iter,sum(stel_data.FLUXLOOPS_chisq,2),'xg','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'SEGROG_chisq')
        note_string=[note_string; ['\chi^2_{SEGROG}=' num2str(sum(stel_data.SEGROG_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'SEGROG'];
        hold on
        plot(stel_data.iter,sum(stel_data.SEGROG_chisq,2),'xb','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'EXTCUR_chisq')
        note_string=[note_string; ['\chi^2_{EXTCUR}=' num2str(sum(stel_data.EXTCUR_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'EXTCUR'];
        hold on
        plot(stel_data.iter,sum(stel_data.EXTCUR_chisq,2),'xc','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'SEPARATRIX_chisq')
        note_string=[note_string; ['\chi^2_{SEPARATRIX}=' num2str(sum(stel_data.SEPARATRIX_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'SEPARATRIX'];
        hold on
        plot(stel_data.iter,sum(stel_data.SEPARATRIX_chisq,2),'xm','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'LIMITER_chisq')
        note_string=[note_string; ['\chi^2_{LIMITER}=' num2str(sum(stel_data.LIMITER_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'LIMITER'];
        hold on
        plot(stel_data.iter,sum(stel_data.LIMITER_chisq,2),'xy','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'BOOTSTRAP_chisq')
        note_string=[note_string; ['\chi^2_{BOOTSTRAP}=' num2str(sum(stel_data.BOOTSTRAP_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'BOOTSTRAP'];
        hold on
        plot(stel_data.iter,sum(stel_data.BOOTSTRAP_chisq,2),'dr','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'BALLOON_chisq')
        note_string=[note_string; ['\chi^2_{BALLOON}=' num2str(sum(stel_data.BALLOON_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'BALLOON'];
        hold on
        plot(stel_data.iter,sum(stel_data.BALLOON_chisq,2),'dg','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'HELICITY_chisq')
        note_string=[note_string; ['\chi^2_{HELICITY}=' num2str(sum(stel_data.HELICITY_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'HELICITY'];
        hold on
        plot(stel_data.iter,sum(stel_data.HELICITY_chisq,2),'db','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'J_STAR_chisq')
        note_string=[note_string; ['\chi^2_{J_STAR}=' num2str(sum(stel_data.J_STAR_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'J_STAR'];
        hold on
        plot(stel_data.iter,sum(stel_data.J_STAR_chisq,2),'dc','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'NEO_chisq')
        note_string=[note_string; ['\chi^2_{NEO}=' num2str(sum(stel_data.NEO_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'NEO'];
        hold on
        plot(stel_data.iter,sum(stel_data.NEO_chisq,2),'dm','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'DKES_chisq')
        note_string=[note_string; ['\chi^2_{DKES}=' num2str(sum(stel_data.DKES_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'DKES'];
        hold on
        plot(stel_data.iter,sum(stel_data.DKES_chisq,2),'dy','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'TXPORT_chisq')
        note_string=[note_string; ['\chi^2_{TXPORT}=' num2str(sum(stel_data.TXPORT_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'TXPORT'];
        hold on
        plot(stel_data.iter,sum(stel_data.TXPORT_chisq,2),'dw','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    if isfield(stel_data,'COIL_BNORM_chisq')
        note_string=[note_string; ['\chi^2_{COILOPT++}=' num2str(sum(stel_data.COIL_BNORM_chisq(i,:)),'%7.4e')]];
        legtext=[legtext 'COILOPT++'];
        hold on
        plot(stel_data.iter,sum(stel_data.COIL_BNORM_chisq,2),'hr','MarkerSize',chisqsiz,'LineWidth',2.0);
        hold off
    end
    % Add yellow line
    set(gca,'Yscale','log');
    hold on
    plot([1 1].*stel_data.iter(i),ylim,'yellow','LineWidth',0.5);
    legtext=[legtext ' '];
    hold off
    % Adust chisq plot scale
    temp=ylim;
    if temp(1) < 1E-2, temp(1) = 1E-2; end
    if temp(2) > 1E6, temp(2) = 1E6; end
    if (temp(1) < temp(2)), ylim(temp);end
    % Tack on some Beta and toroidal current info
    note_string=[note_string; strtrim(sprintf('<\\beta>=%3.2f%%,\\beta_0=%3.2f%%',data.betatot*100.,data.betaxis*100))];
    cur_fac=1.0; cur_lab = 'A';
    if (abs(data.Itor)/1.0E6 > 1.0)
        cur_fac = 1.0E-6; cur_lab = '[MA]';
    elseif (abs(data.Itor)/1.0E3 > 1.0)
        cur_fac = 1.0E-3; cur_lab = '[kA]';
    end
    note_string=[note_string; ['I=' num2str(data.Itor*cur_fac,'%5.1f') cur_lab '  \Phi=' num2str(data.phi(data.ns),'%5.2f') '[Wb]']];
    %note_string=[note_string; ['  \Phi=' num2str(data.phi(data.ns),'%5.2f') '[Wb]']];
    text_handle=annotation('textbox','String',note_string,'Color','white',...
        'Units','pixels','Position',text_pos,'FontSize',fontsiz1);
    ylabel('\chi^2');
    xlabel('Function Evaluations');
    hleg2=legend(legtext,'Location','NorthEastOutside');
    set(hleg2,'TextColor','white','FontSize',fontsiz1);
    axis normal;
    pause(0.001);
    if (lmovie)
        frame=getframe(fig);
        writeVideo(mov,frame);
    end
    % To stop at a given frame
    %if i==1, return, end;
end
% Close Movie File
if (lmovie)
    close(mov);
    % Make a poster frame
    saveas(gca,[movname,'_poster.fig'],'fig');
    hold on
    text(0,1,'DONE','Color','yellow','FontSize',fontsiz1)
    hold off
end
warning('on');
end
