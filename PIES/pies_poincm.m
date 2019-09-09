function fig=pies_poincm(varargin)
%PIES_POINCM Plots the Poincare points files output by PIES.
%   This function plots the Poincaré points files output by PIES (fort.15,
%   fort.16, and fort.17).  If vessel data is passed to the function it
%   will plot the vessel data in addition to the Poincaré plots.  It
%   returns the handle to the current figure.
%
%   Options
%       pies_poincm('REG')          Plot all three cutplanes in one figure. (default)
%       pies_poincm('PHI_0')        Plot the phi=0 cutplane only.
%       pies_poincm('phi_pi2')      Plot the quarter field period cutplane.
%       pies_poincm('phi_pi')       Plot the half field period cutplane.
%       pies_poincm('PHI_0',gcf)    Plot in current figure.
%
%   See also read_vessel, plot_vessel, read_vmec, read_pies_netcdf.
%
%   Example
%       ves_data=read_vessel('test_vessel.dat');
%       pies_poincm('reg',ves_data);
%
%   Written by: Samuel Lazerson (lazerson@pppl.gov)
%   Version:    1.0
%   Date:       05/19/2011

% Set Defaults
i=1;
ves_data=[];
vmec_data=[];
pies_data=[];
pies_magco=[];
data0=[];
datapi2=[];
datapi3=[];
plottype='REG';
leg_text={'Fieldline'};
fig=[];
while i <= nargin
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch upper(varargin{i}.datatype)
                case 'VESSEL'
                    ves_data=varargin{i};
                case 'WOUT'
                    vmec_data=varargin{i};
                case 'PIES_OUT'
                    pies_data=varargin{i};
                case 'PIES_MAGCO'
                    pies_magco=varargin{i};
            end
        end
    elseif ischar(varargin{i})
        switch upper(varargin{i})
            case {'REG','PHI_0','PHI_PI2','PHI_PI'}
                plottype=strtrim(upper(varargin{i}));
        end
    elseif ishghandle(varargin{i})
        fig=varargin{i};
        hold on;
    end
    i=i+1;
end

    

% Plot the window
if isempty(fig)
    fig=figure('Position',[1 1 1920 1080]);
end

% Read the files
switch plottype
    case 'PHI_0'
        % Read first file
        fid=fopen('fort.15','r');
        data0=fscanf(fid,'%e %e orbit # = %d',[3 inf]);
        fclose(fid);
    case 'PHI_PI2'
        % Read second file
        fid=fopen('fort.16','r');
        datapi2=fscanf(fid,'%e %e',[2 inf]);
        fclose(fid);
    case 'PHI_PI'
        % Read third file
        fid=fopen('fort.17','r');
        datapi=fscanf(fid,'%e %e',[2 inf]);
        fclose(fid);
    case {'REG'}
        % Read first file
        fid=fopen('fort.15','r');
        data0=fscanf(fid,'%e %e orbit # = %d',[3 inf]);
        fclose(fid);
        % Read second file
        fid=fopen('fort.16','r');
        datapi2=fscanf(fid,'%e %e',[2 inf]);
        fclose(fid);
        % Read third file
        fid=fopen('fort.17','r');
        datapi=fscanf(fid,'%e %e',[2 inf]);
        fclose(fid);
    otherwise
        disp(['Error: Unknown plottype - ' plottype]);
        close(fig);
        fig=-1;
end

% Make the plots
switch plottype
    case 'PHI_0'
        plot(data0(1,:),data0(2,:),'.k','MarkerSize',1);
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincaré Plot (phi=0)');
        axis equal
    case 'PHI_PI2'
        plot(datapi2(1,:),datapi2(2,:),'.k','MarkerSize',1);
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincaré Plot (phi=\pi/2)');
        axis equal
    case 'PHI_PI'
        plot(datapi(1,:),datapi(2,:),'.k','MarkerSize',1);
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincaré Plot (phi=\pi)');
        axis equal
    case 'REG'
        subplot(1,3,1);
        plot(data0(1,:),data0(2,:),'.k','MarkerSize',1);
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincaré Plot (phi=0)');
        axis equal
        subplot(1,3,2);
        plot(datapi2(1,:),datapi2(2,:),'.k','MarkerSize',1);
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincaré Plot (phi=\pi/2)');
        axis equal
        subplot(1,3,3);
        plot(datapi(1,:),datapi(2,:),'.k','MarkerSize',1);
        xlabel('R [m]');
        ylabel('Z [m]');
        title('Poincaré Plot (phi=\pi)');
        axis equal
end

% Plot the vessel if passed
if ~isempty(ves_data)
    leg_text=[leg_text 'Vessel'];
    switch plottype
        case 'PHI_0'
            hold on
            hves=plot_vessel(ves_data,'phi',0);
            set(hves,'Color','k');
            hold off
        case 'PHI_PI2'
            hold on
            hves=plot_vessel(ves_data,'phi',3*pi/2/10);
            set(hves,'Color','k');
            hold off
        case 'PHI_PI'
            hold on
            hves=plot_vessel(ves_data,'phi',pi/10);
            set(hves,'Color','k');
            hold off
        case 'REG'
            phi=[0 3*pi/2 pi]./10;
            for i=1:3
                subplot(1,3,i);
                hold on
                hves=plot_vessel(ves_data,'phi',phi(i));
                set(hves,'Color','k');
                hold off
            end
    end
end

% Plot the VMEC data if passed
if ~isempty(vmec_data)
    leg_text=[leg_text 'VMEC Boundary'];
    theta=0:2*pi/359:2*pi;
    switch plottype
        case 'PHI_0'
            zeta=0;
            r=cfunct(theta,zeta,vmec_data.rmnc(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            z=sfunct(theta,zeta,vmec_data.zmns(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            hold on
            plot(r,z,'b');
            hold off
        case 'PHI_PI2'
            zeta=pi/2./vmec_data.nfp;
            r=cfunct(theta,zeta,vmec_data.rmnc(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            z=sfunct(theta,zeta,vmec_data.zmns(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            hold on
            plot(r,z,'b');
            hold off
        case 'PHI_PI'
            zeta=pi./vmec_data.nfp;
            r=cfunct(theta,zeta,vmec_data.rmnc(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            z=sfunct(theta,zeta,vmec_data.zmns(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            hold on
            plot(r,z,'b');
            hold off
        case 'REG'
            zeta=[0 pi/2 pi]./vmec_data.nfp;
            r=cfunct(theta,zeta,vmec_data.rmnc(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            z=sfunct(theta,zeta,vmec_data.zmns(:,vmec_data.ns),vmec_data.xm,vmec_data.xn);
            for i=1:3
                subplot(1,3,i);
                hold on
                plot(r(1,:,i),z(1,:,i),'b');
                hold off
            end
    end
end

% Plot the PIES data if passed
if ~isempty(pies_data)
    theta=0:2*pi/359:2*pi;
    leg_text=[leg_text 'PIES Boundary' 'PIES LCFS'];
    switch plottype
        case 'PHI_0'
            zeta=0;
            r=cfunct(-theta,zeta,pies_data.rmnc(:,pies_data.ns),pies_data.xm,pies_data.xn);
            z=sfunct(-theta,zeta,pies_data.zmns(:,pies_data.ns),pies_data.xm,pies_data.xn);
            dex=(data0(3,:)==(pies_data.hitsrf-1));
            hold on
            plot(r,z,'r');
            plot(data0(1,dex),data0(2,dex),'r.');
            hold off
        case 'PHI_PI2'
            zeta=pi/2;
            r=cfunct(-theta,zeta,pies_data.rmnc(:,pies_data.ns),pies_data.xm,pies_data.xn);
            z=sfunct(-theta,zeta,pies_data.zmns(:,pies_data.ns),pies_data.xm,pies_data.xn);
            hold on
            plot(r,z,'r');
            hold off
        case 'PHI_PI'
            zeta=pi;
            r=cfunct(-theta,zeta,pies_data.rmnc(:,pies_data.ns),pies_data.xm,pies_data.xn);
            z=sfunct(-theta,zeta,pies_data.zmns(:,pies_data.ns),pies_data.xm,pies_data.xn);
            hold on
            plot(r,z,'r');
            hold off
        case 'REG'
            zeta=[0 pi/2 pi];
            dex=(data0(3,:)==(pies_data.hitsrf-1));
            rp=cfunct(-theta,zeta,pies_data.rmnc(:,pies_data.ns),pies_data.xm,pies_data.xn);
            zp=sfunct(-theta,zeta,pies_data.zmns(:,pies_data.ns),pies_data.xm,pies_data.xn);
            for i=1:3
                subplot(1,3,i);
                hold on
                plot(rp(1,:,i),zp(1,:,i),'r');
                hold off
            end
            subplot(1,3,1)
            hold on
            plot(data0(1,dex),data0(2,dex),'r.');
            hold off
    end
end

% Plot the PIES magco if passed
if ~isempty(pies_magco)
    theta=0:2*pi/359:2*pi;
    leg_text=[leg_text 'MAG_AXIS' 'PIES GOOD SURFACES'];
    switch plottype
        case 'PHI_0'
            zeta=0;
            r=cfunct(-theta,zeta,pies_magco.rmnc,pies_magco.xm,pies_magco.xn);
            z=sfunct(-theta,zeta,pies_magco.zmns,pies_magco.xm,pies_magco.xn);
            hold on
            plot(r(1,1),z(1,1),'+g');
            if isfield(pies_magco,'hitsrf')
                k=pies_magco.hitsrf-1;
            else
                k=pies_magco.ns;
            end
            for i=2:k
                plot(r(i,:),z(i,:),'g');
            end
            hold off
        case 'PHI_PI2'
            zeta=pi/2;
            r=cfunct(-theta,zeta,pies_magco.rmnc,pies_magco.xm,pies_magco.xn);
            z=sfunct(-theta,zeta,pies_magco.zmns,pies_magco.xm,pies_magco.xn);
            hold on
            plot(r(1,1),z(1,1),'+g');
            if isfield(pies_magco,'hitsrf')
                k=pies_magco.hitsrf-1;
            else
                k=pies_magco.ns;
            end
            for i=2:k
                plot(r(i,:),z(i,:),'g');
            end
            hold off
        case 'PHI_PI'
            zeta=pi;
            r=cfunct(-theta,zeta,pies_magco.rmnc,pies_magco.xm,pies_magco.xn);
            z=sfunct(-theta,zeta,pies_magco.zmns,pies_magco.xm,pies_magco.xn);
            hold on
            plot(r(1,1),z(1,1),'+g');
            if isfield(pies_magco,'hitsrf')
                k=pies_magco.hitsrf-1;
            else
                k=pies_magco.ns;
            end
            for i=2:k
                plot(r(i,:),z(i,:),'g');
            end
            hold off
        case 'REG'
            zeta=[0 pi/2 pi];
            r=cfunct(-theta,zeta,pies_magco.rmnc,pies_magco.xm,pies_magco.xn);
            z=sfunct(-theta,zeta,pies_magco.zmns,pies_magco.xm,pies_magco.xn);
            for i=1:3
                subplot(1,3,i);
                hold on
                plot(r(1,1,i),z(1,1,i),'+g');
                if isfield(pies_magco,'hitsrf')
                    k=pies_magco.hitsrf-1;
                else
                    k=pies_magco.ns;
                end
                for j=2:k
                    plot(r(j,:,i),z(j,:,i),'g');
                end
                hold off
            end
            subplot(1,3,1)
            hold on
            plot(data0(1,dex),data0(2,dex),'r.');
            hold off
    end
end

% Add legend
legend(leg_text);

end

