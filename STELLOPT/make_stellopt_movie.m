function make_stellopt_movie(varargin)
%MAKE_STELLOPT_MOVIE(movietype[,vessel_data]) Create stellopt movie
%   MAKE_STELLOPT_MOVIE(movietype) Creates a movie from stellopt output
%   files.  The routine supports 3 resolutions (1080, 720, 480; 1080 is
%   default).  If vessel data is provided then all VMEC plots will include
%   that data.  The plot options are:
%       'opt':      (default) R-Z plot of the VMEC vaccum surfaces at phi=0
%       'pressure': Pressure fitting of pressure profile.
%       'recon':    Reconstruction showing R-Z at multiple planes, pressure
%                   profile fitting and chi minimization.
%       'ncsx':     Set views suitable for NCSX
%       'lhd':      Set views suitable for LHD (default)
%
%   Example:
%       ves_data=read_vessel('vessel.dat');
%       make_stellopt_movie('1080','recon',ves_data);
%
%   See also plot_stellopt, read_stellopt, and read_vessel.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/17/11

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
            end
        else
            switch varargin{i}
                case {'pressure','opt','recon','current','mse',...
                        'thomson_te','thomson_ne','ti'}
                    plottype=varargin{i};
                case {'1080','720','480'}
                    restype=varargin{i};
                case {'realspace','fluxspace'}
                    axistype=varargin{i};
                case {'old'}
                    lold_data=1;
                case {'ncsx'}
                    vmec_xmin=0.75; vmec_xmax=2.25;
                    vmec_ymin=-0.75; vmec_ymax=0.75;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                case {'D3D','DIIID'}
                    vmec_xmin=0.1; vmec_xmax=2.9;
                    vmec_ymin=-1.4; vmec_ymax=1.4;
                    dx=(vmec_xmax-vmec_xmin)/20;  % 5%
                    dy=(vmec_ymax-vmec_ymin)/100;  % 1%
                case {'isote'}
                    isote = 1;
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
if strcmp(restype,'720')
    win_res=[1 1 1280 720];
elseif strcmp(restype,'480')
    win_res=[1 1 640 480];
else
    win_res=[1 1 1920 1080];
    text_pos=[1020 100 100 200];
end
fig=figure('DoubleBuffer','on','Position',win_res,...
        'Color','black','BackingStore','on','MenuBar','none',...
        'Name',plottype,'InvertHardcopy','off');
set(gca,'nextplot','replacechildren','XColor','white',...
        'YColor','white','ZColor','white','Color','black','FontSize',20);
% Get File lists
o_filelist=dir('output.*');
if lold_data
    wout_filename='wout.*.*';
else
    wout_filename='wout_*.*.nc';
end
switch plottype
    case {'pressure','thomson_te','thomson_ne','ti'}
        filelist=dir(wout_filename);
        p_filelist=dir('p_prof.*.*');
        input_filelist=dir('input.*.*');
    case 'mse'
        filelist=dir(wout_filename);
        mse_filelist=dir('mse_prof.*');
        input_filelist=dir('input.*.*');
    case 'current'
        filelist=dir(wout_filename);
        j_filelist=dir('mag_diags.*.*');
        input_filelist=dir('input.*.*');
    case {'opt'}
        filelist=dir(wout_filename);
        input_filelist=dir('input.*.*');
    case 'recon'
        filelist=dir(wout_filename);
        p_filelist=dir('p_prof.*.*');
        input_filelist=dir('input.*.*');
        j_filelist=dir('mag_diags.*.*');
end

% Process File lists
o_filename={o_filelist(1).name};
filename={filelist(1).name};
input_filename={input_filelist(1).name};
nfiles=numel(filelist);
for i=2:nfiles
    filename=[filename; filelist(i).name];
    input_filename=[input_filename; input_filelist(i).name];
end
switch plottype
    case {'opt'}
    case {'pressure','thomson_te','thomson_ne','ti'}
        p_filename={p_filelist(1).name};
        for i=2:nfiles
            p_filename=[p_filename; p_filelist(i).name];
        end
    case {'mse'}
        mse_filename={mse_filelist(1).name};
        for i=2:nfiles
            mse_filename=[mse_filename; mse_filelist(i).name];
        end
    case {'current'}
        j_filename={j_filelist(1).name};
        for i=2:nfiles
            j_filename=[j_filename; j_filelist(i).name];
        end
    case {'recon'}
        p_filename={p_filelist(1).name};
        j_filename={j_filelist(1).name};
        for i=2:nfileso_filelist
            p_filename=[p_filename; p_filelist(i).name];
            j_filename=[j_filename; j_filelist(i).name];
        end        
end

% Creat movie filename
temp=o_filename{1};
dex1=strfind(temp,'.');
dex2=length(temp);
movname=strcat(plottype,temp(dex1:dex2));
if strcmp(axistype,'fluxspace')
    movname=strcat(movname,'_flux');
end
% Setup Movie File
mov=VideoWriter([movname '_movie.avi']);
mov.Quality= 100;
mov.FrameRate=1;
open(mov);

filedex=[];
% Get file numbers
for i=1:nfiles-1
    dex=strfind(input_filename{i},'.');
    dex=dex(numel(dex))+1;
    filedex=[filedex str2double(input_filename{i}(dex:numel(input_filename{i})))];
end
if (nfiles > 1)
    filedex=[filedex filedex(nfiles-1)];  % To Get the minimum element
end

% Read Output File
odata=read_stellopt(o_filename{1});
if ~isempty(filedex)
    for j=1:numel(filedex)
        filedex(j)=find(odata.func_eval==filedex(j),1);
    end
else
    filedex(1) = 1;
end

% Read Input File for nproc


% Plot Files
switch plottype
    case 'pressure'
        for i=1:nfiles
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            p=zeros(data.ns,361,3);
            for j=1:data.ns, p(j,:,:)=data.presf(j); end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',10,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','bottom');
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            pdata=read_stellopt(p_filename{i});
            %Find VMEC DOMAIN
            dex1=find(pdata.data(:,4)<1,1,'first');
            dex2=find(pdata.data(:,4)<1,1,'last');
            x1=pdata.data(dex1,1);
            x2=pdata.data(dex2,1);
            thom_data=pdata.data(:,7).*pdata.norm;    % Factor * p-data
            vmec_pres=pdata.data(:,8);
            sigmas=pdata.norm.*pdata.data(:,9);
            %sigmas=pdata.norm.*pdata.data(:,9)./sqrt(length(thom_data));
            thom_wegt=(vmec_pres-thom_data)./sigmas;
            y1=max(thom_data)*.01;
            switch axistype
                case 'realspace'
                    vmec_line=zeros(1,dex2-dex1+1);
                    hp1=plot(pdata.data(:,1),thom_data,'ow');
                    xlim([2.5 5]);
                    hold on
                    hpatch=fill([x1 x2 x2 x1],[0 0 y1 y1],'red');
                    set(hpatch,'FaceColor','red','EdgeColor','red');
                    hp2=plot(pdata.data(:,1),vmec_pres,'w');
                    text((x2-x1)/3+x1,y1*5,'VMEC Domain','Color','white');
                    hp3=plot([r(1,1,3) r(1,1,3)],[0 max([thom_data; vmec_pres])],'Color','yellow');
                    hold off
                    axis tight
                    ylabel('Pressure [Pa]');
                    ylim([0 max([max(thom_data) max(vmec_pres)])]);
                    hleg=legend([hp1,hp2],'Experiment','Simulation','boxoff');
                    p_xlim=xlim;
                    set(hleg,'TextColor','white','FontSize',12);
                    % Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    pdata.data(:,9)=pdata.data(:,9)*sqrt(length(pdata.data(:,9)));  % Renormalization factor
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(pdata.data(:,1),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(pdata.data(:,1),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(pdata.data(:,1),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('R[m]');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([p_xlim(1) p_xlim(2)]);
                    ylim([-5 5]);
                case 'fluxspace'
                    hp1=plot(pdata.data(:,4),thom_data,'ow');
                    if (max(xlim) > 2.0) xlim([0 2.0]); end
                    hold on
                    hp2=plot(data.phi./data.phi(data.ns),data.presf./max(data.presf),'w');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    % Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    pdata.data(:,9)=pdata.data(:,9)*sqrt(length(pdata.data(:,9)));  % Renormalization factor
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(1:size(pdata.data,1),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:size(pdata.data,1),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:size(pdata.data,1),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 size(pdata.data,1)]);
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(pdata.data(:,1),pdata.data(:,2),'+g');
            end
            % Plot Output
            subplot(4,3,[7 8 10 11]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            %colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot([odata.func_eval(filedex(i)) odata.func_eval(filedex(i))],[1e-7 max(odata.chi_total)],'yellow','LineWidth',4.0);
            legtext={''};
            hold on
            plot(odata.func_eval(filedex),odata.chi_total(filedex),'o','Color','white',...
                'MarkerSize',10,'LineWidth',2.0);
            legtext=[legtext 'Total'];
            note_string={['\chi^2_{Total}=' num2str(odata.chi_total(filedex(i)),'%7.4e')]};
            plot(odata.func_eval(filedex),odata.chisq_PressureProfile(filedex),'.r','MarkerSize',20,'LineWidth',2.0);
            legtext=[legtext 'Pressure Profile'];
            note_string=[note_string; ['\chi^2_{Pressure}=' num2str(odata.chisq_PressureProfile(filedex(i)),'%7.4e')]];
            %plot(odata.func_eval(filedex),odata.chisq_PressureProfile(filedex)./odata.nw_PressureProfile(filedex),'xr','MarkerSize',20,'LineWidth',2.0);
            plot(odata.func_eval(filedex),odata.chisq_EdgePressure(filedex),'.b','MarkerSize',20,'LineWidth',2.0);
            legtext=[legtext 'Edge Pressure'];
            note_string=[note_string; ['\chi^2_{P-Edge}=' num2str(odata.chisq_EdgePressure(filedex(i)),'%7.4e')]];
            hold off
            % Add P_Grad String
            if isfield(odata,'chisq_PressureGradient')
                note_string=[note_string; ['\chi^2_{<P-Grad>}=' num2str(odata.chisq_PressureGradient(filedex(i)),'%7.4e')]];
                legtext=[legtext 'P-Grad'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_PressureGradient(filedex),'+b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add CURTOR if available
            if isfield(odata,'chisq_Totaltor_current')
                note_string=[note_string; ['\chi^2_{I-tor}=' num2str(odata.chisq_Totaltor_current(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Toroidal Current'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_Totaltor_current(filedex),'.y','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add DIAGNO if available
            if isfield(odata,'chisq_DIAGNOmag_diag_')
                note_string=[note_string; ['\chi^2_{Mag. Diag.}=' num2str(odata.chisq_DIAGNOmag_diag_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Mag. Diags.'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_DIAGNOmag_diag_(filedex),'.m','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Beta String
            if isfield(odata,'chisq_VolAvgBeta')
                note_string=[note_string; ['\chi^2_{<Beta>}=' num2str(odata.chisq_VolAvgBeta(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Beta'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_VolAvgBeta(filedex),'.c','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add MSE String
            if isfield(odata,'chisq_MSEDiagnostic')
                note_string=[note_string; ['\chi^2_{<MSE>}=' num2str(odata.chisq_MSEDiagnostic(filedex(i)),'%7.4e')]];
                legtext=[legtext 'MSE'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_MSEDiagnostic(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Tack on some Beta and toroidal current info
            note_string=[note_string; strtrim(sprintf('<\\beta>=%3.2f%%,\\beta_0=%3.2f%%',odata.avgbeta(filedex(i))*100,data.betaxis*100))];
            note_string=[note_string; ['I=' num2str(data.Itor*1e-3,'%5.1f') '[kA]']];
            text_handle=annotation('textbox','String',note_string,'Color','white',...
                'Units','pixels','Position',text_pos,'FontSize',20);
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend(legtext,'Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',20);
            %ylim([0.001 max(odata.chi_total(filedex))]);
            ylim([1e-7 10000]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    case 'thomson_te'
        for i=1:nfiles
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            p=zeros(data.ns,361,3);
            for j=1:data.ns, p(j,:,:)=data.presf(j); end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',10,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','bottom');
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            pdata=read_stellopt(p_filename{i});
            %Find VMEC DOMAIN
            dex1=find(pdata.te_data(:,4)<1,1,'first');
            dex2=find(pdata.te_data(:,4)<1,1,'last');
            thom_r=pdata.te_data(:,1);
            thom_phi=pdata.te_data(:,3);
            thom_z=pdata.te_data(:,2);
            thom_s=pdata.te_data(:,4);
            thom_exp=pdata.te_data(:,5);
            thom_sim=pdata.te_data(:,6);
            thom_sigma=pdata.te_data(:,7);
            thom_wegt=pdata.te_data(:,8);
            x1=pdata.te_data(dex1,1);
            x2=pdata.te_data(dex2,1);
            switch axistype
                case 'realspace'
                case 'fluxspace'
                    errorbar(thom_s,thom_exp,thom_sigma,'ow')
                    xlim([0 1.2]);
                    hold on
                    plot(thom_s,thom_sim,'or');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 size(thom_r,1)]);
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(thom_r,thom_z,'+g');
            end
            % Plot Output
            subplot(4,3,[7 8 10 11]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            %colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot([odata.func_eval(filedex(i)) odata.func_eval(filedex(i))],[1e-7 max(odata.chi_total)],'yellow','LineWidth',4.0);
            legtext={''};
            hold on
            plot(odata.func_eval(filedex),odata.chi_total(filedex),'o','Color','white',...
                'MarkerSize',10,'LineWidth',2.0);
            legtext=[legtext 'Total'];
            note_string={['\chi^2_{Total}=' num2str(odata.chi_total(filedex(i)),'%7.4e')]};
            hold off
            % Add Ne String
            if isfield(odata,'chisq_ElectronDensity')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_ElectronDensity(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Ne'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_ElectronDensity(filedex),'.r','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Te String
            if isfield(odata,'chisq_ElectronTemp_')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_ElectronTemp_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Te'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_ElectronTemp_(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Ti String
            if isfield(odata,'chisq_IonTemperature')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_IonTemperature(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Ti'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_IonTemperature(filedex),'.b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add P_Grad String
            if isfield(odata,'chisq_PressureGradient')
                note_string=[note_string; ['\chi^2_{<P-Grad>}=' num2str(odata.chisq_PressureGradient(filedex(i)),'%7.4e')]];
                legtext=[legtext 'P-Grad'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_PressureGradient(filedex),'+b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add CURTOR if available
            if isfield(odata,'chisq_Totaltor_current')
                note_string=[note_string; ['\chi^2_{I-tor}=' num2str(odata.chisq_Totaltor_current(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Toroidal Current'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_Totaltor_current(filedex),'.y','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add DIAGNO if available
            if isfield(odata,'chisq_DIAGNOmag_diag_')
                note_string=[note_string; ['\chi^2_{Mag. Diag.}=' num2str(odata.chisq_DIAGNOmag_diag_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Mag. Diags.'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_DIAGNOmag_diag_(filedex),'.m','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Beta String
            if isfield(odata,'chisq_VolAvgBeta')
                note_string=[note_string; ['\chi^2_{<Beta>}=' num2str(odata.chisq_VolAvgBeta(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Beta'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_VolAvgBeta(filedex),'.c','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add MSE String
            if isfield(odata,'chisq_MSEDiagnostic')
                note_string=[note_string; ['\chi^2_{<MSE>}=' num2str(odata.chisq_MSEDiagnostic(filedex(i)),'%7.4e')]];
                legtext=[legtext 'MSE'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_MSEDiagnostic(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Tack on some Beta and toroidal current info
            note_string=[note_string; strtrim(sprintf('<\\beta>=%3.2f%%,\\beta_0=%3.2f%%',odata.avgbeta(filedex(i))*100,data.betaxis*100))];
            note_string=[note_string; ['I=' num2str(data.Itor*1e-3,'%5.1f') '[kA]']];
            text_handle=annotation('textbox','String',note_string,'Color','white',...
                'Units','pixels','Position',text_pos,'FontSize',20);
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend(legtext,'Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',20);
            %ylim([0.001 max(odata.chi_total(filedex))]);
            ylim([1e-7 10000]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    case 'thomson_ne'
        for i=1:nfiles
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            p=zeros(data.ns,361,3);
            for j=1:data.ns, p(j,:,:)=data.presf(j); end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',10,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','bottom');
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            pdata=read_stellopt(p_filename{i});
            %Find VMEC DOMAIN
            dex1=find(pdata.ne_data(:,4)<1,1,'first');
            dex2=find(pdata.ne_data(:,4)<1,1,'last');
            thom_r=pdata.ne_data(:,1);
            thom_phi=pdata.ne_data(:,3);
            thom_z=pdata.ne_data(:,2);
            thom_s=pdata.ne_data(:,4);
            thom_exp=pdata.ne_data(:,5);
            thom_sim=pdata.ne_data(:,6);
            thom_sigma=pdata.ne_data(:,7);
            thom_wegt=pdata.ne_data(:,8);
            x1=pdata.ne_data(dex1,1);
            x2=pdata.ne_data(dex2,1);
            switch axistype
                case 'realspace'
                case 'fluxspace'
                    errorbar(thom_s,thom_exp,thom_sigma,'ow')
                    xlim([0 1.2]);
                    hold on
                    plot(thom_s,thom_sim,'or');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 size(thom_r,1)]);
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(thom_r,thom_z,'+g');
            end
            % Plot Output
            subplot(4,3,[7 8 10 11]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            %colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot([odata.func_eval(filedex(i)) odata.func_eval(filedex(i))],[1e-7 max(odata.chi_total)],'yellow','LineWidth',4.0);
            legtext={''};
            hold on
            plot(odata.func_eval(filedex),odata.chi_total(filedex),'o','Color','white',...
                'MarkerSize',10,'LineWidth',2.0);
            legtext=[legtext 'Total'];
            note_string={['\chi^2_{Total}=' num2str(odata.chi_total(filedex(i)),'%7.4e')]};
            hold off
            % Add Ne String
            if isfield(odata,'chisq_ElectronDensity')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_ElectronDensity(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Ne'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_ElectronDensity(filedex),'.r','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Te String
            if isfield(odata,'chisq_ElectronTemp_')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_ElectronTemp_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Te'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_ElectronTemp_(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Ti String
            if isfield(odata,'chisq_IonTemperature')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_IonTemperature(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Ti'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_IonTemperature(filedex),'.b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add P_Grad String
            if isfield(odata,'chisq_PressureGradient')
                note_string=[note_string; ['\chi^2_{<P-Grad>}=' num2str(odata.chisq_PressureGradient(filedex(i)),'%7.4e')]];
                legtext=[legtext 'P-Grad'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_PressureGradient(filedex),'+b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add CURTOR if available
            if isfield(odata,'chisq_Totaltor_current')
                note_string=[note_string; ['\chi^2_{I-tor}=' num2str(odata.chisq_Totaltor_current(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Toroidal Current'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_Totaltor_current(filedex),'.y','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add DIAGNO if available
            if isfield(odata,'chisq_DIAGNOmag_diag_')
                note_string=[note_string; ['\chi^2_{Mag. Diag.}=' num2str(odata.chisq_DIAGNOmag_diag_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Mag. Diags.'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_DIAGNOmag_diag_(filedex),'.m','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Beta String
            if isfield(odata,'chisq_VolAvgBeta')
                note_string=[note_string; ['\chi^2_{<Beta>}=' num2str(odata.chisq_VolAvgBeta(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Beta'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_VolAvgBeta(filedex),'.c','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add MSE String
            if isfield(odata,'chisq_MSEDiagnostic')
                note_string=[note_string; ['\chi^2_{<MSE>}=' num2str(odata.chisq_MSEDiagnostic(filedex(i)),'%7.4e')]];
                legtext=[legtext 'MSE'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_MSEDiagnostic(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Tack on some Beta and toroidal current info
            note_string=[note_string; strtrim(sprintf('<\\beta>=%3.2f%%,\\beta_0=%3.2f%%',odata.avgbeta(filedex(i))*100,data.betaxis*100))];
            note_string=[note_string; ['I=' num2str(data.Itor*1e-3,'%5.1f') '[kA]']];
            text_handle=annotation('textbox','String',note_string,'Color','white',...
                'Units','pixels','Position',text_pos,'FontSize',20);
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend(legtext,'Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',20);
            %ylim([0.001 max(odata.chi_total(filedex))]);
            ylim([1e-7 10000]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    case 'mse'
        for i=1:nfiles
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            jcurv=cfunct(theta,zeta,data.currvmnc,data.xm,data.xn);
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                torocont(r,z,jcurv,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',10,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','bottom');
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                colormap('jet');
            end
            % Plot MSE profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            msedata=read_stellopt(mse_filename{i});
            mse_sigma=msedata.data(:,8);
            dex = (mse_sigma < 1e10);
            mse_pol=msedata.data(dex,5);
            mse_vmec=msedata.data(dex,7);
            mse_sigma=msedata.data(dex,8);
            mse_s=msedata.data(dex,4);
            if isfield(data,'acauxs')
                nels = sum(data.acauxs >0);
                acauxs=data.acauxs(1:nels+1);
                acauxf=data.acauxf(1:nels+1);
                %plot(acauxs,acauxs.*0,'+r');
                plot(acauxs,acauxf,'r');
                ylim([min(acauxf) max(acauxf)]);
            else
                hp1=plot(0:1/(data.ns-1):1,data.jcurv./1000,'w');
                hold on
                plot([1 1],[min(data.jcurv) max(data.jcurv)]./1000,'r');
                xlim([0 1]);
                ylim([min(data.jcurv)./1000 max(data.jcurv)./1000]);
            end
            %set(gca,'XScale','log');
            ylabel('Toridal Current profile [kA]');
            % Plot MSE Signal
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            hp2=errorbar(mse_s,180*mse_pol/pi,180*mse_sigma/pi,'wx');  % Measured
            %set(gca,'XTick',1:size(msedata.data,1));
            hold on
            plot(mse_s,180*mse_vmec/pi,'or');
            ylim([-10 10]);
            %plot([1 1],180.*[min(mse_pol-max(mse_sigma)) max(mse_pol+max(mse_sigma))]./pi,'r');
            %ylim(180.*[min(mse_pol-max(mse_sigma)) max(mse_pol+max(mse_sigma))]./pi);
            if isfield(data,'acauxs')
                acauxs=data.acauxs(data.acauxs >= 0);
                plot(acauxs,acauxs.*0,'+r');
            end
            %ylim([min([min(mse_pol) min(mse_vmec)]) max([max(mse_pol) max(mse_vmec)])]);
            xlim([0 1]);
            title('MSE Polarization');
            ylabel('Degrees');
            xlabel('Normalized Flux');
            %set(gca,'XScale','log');
            % Plot Output
            subplot(4,3,[7 8 10 11]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            %colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot([odata.func_eval(filedex(i)) odata.func_eval(filedex(i))],[1e-7 max(odata.chi_total)],'yellow','LineWidth',4.0);
            legtext={''};
            hold on
            plot(odata.func_eval(filedex),odata.chi_total(filedex),'o','Color','white',...
                'MarkerSize',10,'LineWidth',2.0);
            legtext=[legtext 'Total'];
            note_string={['\chi^2_{Total}=' num2str(odata.chi_total(filedex(i)),'%7.4e')]};
            % Add P String
            if isfield(odata,'chisq_PressureProfile')
                note_string=[note_string; ['\chi^2_{Pressure}=' num2str(odata.chisq_PressureProfile(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Pressure Profile'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_PressureProfile(filedex),'.r','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add P String
            if isfield(odata,'chisq_EdgePressure')
                note_string=[note_string; ['\chi^2_{P-Edge}=' num2str(odata.chisq_EdgePressure(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Edge Pressure'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_EdgePressure(filedex),'.b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add P_Grad String
            if isfield(odata,'chisq_PressureGradient')
                note_string=[note_string; ['\chi^2_{<P-Grad>}=' num2str(odata.chisq_PressureGradient(filedex(i)),'%7.4e')]];
                legtext=[legtext 'P-Grad'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_PressureGradient(filedex),'+b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add CURTOR if available
            if isfield(odata,'chisq_Totaltor_current')
                note_string=[note_string; ['\chi^2_{I-tor}=' num2str(odata.chisq_Totaltor_current(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Toroidal Current'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_Totaltor_current(filedex),'.y','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add DIAGNO if available
            if isfield(odata,'chisq_DIAGNOmag_diag_')
                note_string=[note_string; ['\chi^2_{Mag. Diag.}=' num2str(odata.chisq_DIAGNOmag_diag_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Mag. Diags.'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_DIAGNOmag_diag_(filedex),'.m','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Beta String
            if isfield(odata,'chisq_VolAvgBeta')
                note_string=[note_string; ['\chi^2_{<Beta>}=' num2str(odata.chisq_VolAvgBeta(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Beta'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_VolAvgBeta(filedex),'.c','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add MSE String
            if isfield(odata,'chisq_MSEDiagnostic')
                note_string=[note_string; ['\chi^2_{<MSE>}=' num2str(odata.chisq_MSEDiagnostic(filedex(i)),'%7.4e')]];
                legtext=[legtext 'MSE'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_MSEDiagnostic(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Tack on some Beta and toroidal current info
            note_string=[note_string; strtrim(sprintf('<\\beta>=%3.2f%%,\\beta_0=%3.2f%%',odata.avgbeta(filedex(i))*100,data.betaxis*100))];
            note_string=[note_string; ['I=' num2str(data.Itor*1e-3,'%5.1f') '[kA]']];
            text_handle=annotation('textbox','String',note_string,'Color','white',...
                'Units','pixels','Position',text_pos,'FontSize',20);
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend(legtext,'Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',20);
            %ylim([0.001 max(odata.chi_total(filedex))]);
            ylim([1e-7 10000]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    case 'ti'
        for i=1:nfiles
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            if (data.iasym == 1)
                r=r+sfunct(theta,zeta,data.rmns,data.xm,data.xn);
                z=z+cfunct(theta,zeta,data.zmnc,data.xm,data.xn);
            end
            % Make 3D pressure
            p=zeros(data.ns,361,3);
            for j=1:data.ns, p(j,:,:)=data.presf(j); end
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                torocont(r,z,p,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',10,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),vmec_ymin+dy,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),vmec_ymin+dy,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','bottom');
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                % Add title to second plot
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(vmec_xmin+dx,vmec_ymax,'0 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==2, text(vmec_xmin+dx,vmec_ymax,'1/4 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==3, text(vmec_xmin+dx,vmec_ymax,'1/2 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                colormap('hot');
            end
            % Plot Pressure profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            pdata=read_stellopt(p_filename{i});
            %Find VMEC DOMAIN
            dex1=find(pdata.ti_data(:,4)<1,1,'first');
            dex2=find(pdata.ti_data(:,4)<1,1,'last');
            thom_r=pdata.ti_data(:,1);
            thom_phi=pdata.ti_data(:,3);
            thom_z=pdata.ti_data(:,2);
            thom_s=pdata.ti_data(:,4);
            thom_exp=pdata.ti_data(:,5);
            thom_sim=pdata.ti_data(:,6);
            thom_sigma=pdata.ti_data(:,7);
            thom_wegt=pdata.ti_data(:,8);
            x1=pdata.ti_data(dex1,1);
            x2=pdata.ti_data(dex2,1);
            switch axistype
                case 'realspace'
                case 'fluxspace'
                    errorbar(thom_s,thom_exp,thom_sigma,'ow')
                    xlim([0 1.2]);
                    hold on
                    plot(thom_s,thom_sim,'or');
                    hold off
                    xlabel('Normalized Flux');
                    if (min(ylim) < 0.0), ylim([0 max(ylim)]); end
                    %Plot Pressure Residuals
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    reddex=abs(thom_wegt)>2.0;
                    yeldex=(abs(thom_wegt)>1.0).*(abs(thom_wegt)<=2.0);
                    gredex=abs(thom_wegt)<=1.0;
                    if (max(reddex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*reddex,'r','EdgeColor','red');
                        hold on
                    end
                    if (max(yeldex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*yeldex,'y','EdgeColor','yellow');
                        hold on
                    end
                    if (max(gredex) > 0)
                        bar(1:size(thom_r,1),thom_wegt.*gredex,'g','EdgeColor','green');
                        hold off
                    end
                    hold off
                    title('Weighted Deviation','Color','white')
                    xlabel('Index');
                    ylabel('(p_{sim}-p_{exp})/\sigma');
                    xlim([0 size(thom_r,1)]);
                    ylim([-5 5]);
                    subplot(4,3,[3 6]);
                    hold on
                    plot(thom_r,thom_z,'+g');
            end
            % Plot Output
            subplot(4,3,[7 8 10 11]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            %colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot([odata.func_eval(filedex(i)) odata.func_eval(filedex(i))],[1e-7 max(odata.chi_total)],'yellow','LineWidth',4.0);
            legtext={''};
            hold on
            plot(odata.func_eval(filedex),odata.chi_total(filedex),'o','Color','white',...
                'MarkerSize',10,'LineWidth',2.0);
            legtext=[legtext 'Total'];
            note_string={['\chi^2_{Total}=' num2str(odata.chi_total(filedex(i)),'%7.4e')]};
            hold off
            % Add Ne String
            if isfield(odata,'chisq_ElectronDensity')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_ElectronDensity(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Ne'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_ElectronDensity(filedex),'.r','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Te String
            if isfield(odata,'chisq_ElectronTemp_')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_ElectronTemp_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Te'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_ElectronTemp_(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Ti String
            if isfield(odata,'chisq_IonTemperature')
                note_string=[note_string; ['\chi^2_{Ne}=' num2str(odata.chisq_IonTemperature(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Ti'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_IonTemperature(filedex),'.b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add P_Grad String
            if isfield(odata,'chisq_PressureGradient')
                note_string=[note_string; ['\chi^2_{<P-Grad>}=' num2str(odata.chisq_PressureGradient(filedex(i)),'%7.4e')]];
                legtext=[legtext 'P-Grad'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_PressureGradient(filedex),'+b','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add CURTOR if available
            if isfield(odata,'chisq_Totaltor_current')
                note_string=[note_string; ['\chi^2_{I-tor}=' num2str(odata.chisq_Totaltor_current(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Toroidal Current'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_Totaltor_current(filedex),'.y','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add DIAGNO if available
            if isfield(odata,'chisq_DIAGNOmag_diag_')
                note_string=[note_string; ['\chi^2_{Mag. Diag.}=' num2str(odata.chisq_DIAGNOmag_diag_(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Mag. Diags.'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_DIAGNOmag_diag_(filedex),'.m','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add Beta String
            if isfield(odata,'chisq_VolAvgBeta')
                note_string=[note_string; ['\chi^2_{<Beta>}=' num2str(odata.chisq_VolAvgBeta(filedex(i)),'%7.4e')]];
                legtext=[legtext 'Beta'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_VolAvgBeta(filedex),'.c','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Add MSE String
            if isfield(odata,'chisq_MSEDiagnostic')
                note_string=[note_string; ['\chi^2_{<MSE>}=' num2str(odata.chisq_MSEDiagnostic(filedex(i)),'%7.4e')]];
                legtext=[legtext 'MSE'];
                hold on
                plot(odata.func_eval(filedex),odata.chisq_MSEDiagnostic(filedex),'.g','MarkerSize',20,'LineWidth',2.0);
                hold off
            end
            % Tack on some Beta and toroidal current info
            note_string=[note_string; strtrim(sprintf('<\\beta>=%3.2f%%,\\beta_0=%3.2f%%',odata.avgbeta(filedex(i))*100,data.betaxis*100))];
            note_string=[note_string; ['I=' num2str(data.Itor*1e-3,'%5.1f') '[kA]']];
            text_handle=annotation('textbox','String',note_string,'Color','white',...
                'Units','pixels','Position',text_pos,'FontSize',20);
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend(legtext,'Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',20);
            %ylim([0.001 max(odata.chi_total(filedex))]);
            ylim([1e-7 10000]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    case 'current'
        for i=1:nfiles
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Process VMEC data
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            % Make toroidal current
            currv=cfunct(theta,zeta,data.currvmnc,data.nfp,data.xm,data.xn);
            jphi=r.*currv;
            % Plot VMEC data
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                torocont(r,z,jphi,j);
                hold on
                plot(r(1,1,j),z(1,1,j),'+k','MarkerSize',10,'LineWidth',2);
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),-1.45,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),-1.45,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','bottom');
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([-1.5 1.5]);
                xlim([2.5 5]);
                % Add title to second plot
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
                % Add text to each plot
                if j==1, text(2.55,1.5,'0 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==2, text(2.55,1.5,'1/4 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                if j==3, text(2.55,1.5,'1/2 Field Period','Color','white',...
                        'FontSize',20,'VerticalAlignment','top'); end
                colormap('jet');
            end
            % Plot Current profile
            subplot(4,3,9);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            hp1=plot(0:1/(data.ns-1):1,data.jcurv,'w');
            xlabel('Normalized Flux');
            axis tight
            ylabel('Toroidal Current [A]');
            % Plot Magnetic Diagnostics
            jdata=read_stellopt(j_filename{i});
            subplot(4,3,12);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            jdata.data(:,4)=jdata.data(:,4).*sqrt(length(jdata.data(:,4)));  % Renormalization factor
            reddex=abs(jdata.data(:,4))>2.0;
            yeldex=(abs(jdata.data(:,4))>1.0).*(abs(jdata.data(:,4))<=2.0);
            gredex=abs(jdata.data(:,4))<=1.0;
            if (max(reddex) > 0)
                bar(jdata.data(:,4).*reddex,'r','EdgeColor','red');
                hold on
            end
            if (max(yeldex) > 0)
                bar(jdata.data(:,4).*yeldex,'y','EdgeColor','yellow');
                hold on
            end
            if (max(gredex) > 0)
                bar(jdata.data(:,4).*gredex,'g','EdgeColor','green');
                hold off
            end
            hold off
            title('Weighted Deviation','Color','white')
            xlabel('Magnetic Diagnostic');
            ylim([-5 5]);
            % Plot Output
            subplot(4,3,[7 8 10 11]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            %colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot([odata.func_eval(filedex(i)) odata.func_eval(filedex(i))],[0.01 max(odata.chi_total)],'yellow','LineWidth',4.0);
            hold on
            plot(odata.func_eval(filedex),odata.chi_total(filedex),'o','Color','white',...
                'MarkerSize',10,'LineWidth',2.0);
            plot(odata.func_eval(filedex),odata.chisq_DIAGNOmag_diag_(filedex),'.r','MarkerSize',10,'LineWidth',2.0);
            plot(odata.func_eval(filedex),odata.chisq_Totaltor_current(filedex),'.b','MarkerSize',10,'LineWidth',2.0);
            hold off
            % Create Chi String
            note_string={['\chi^2_{Total}=' num2str(odata.chi_total(filedex(i)),'%7.4e')]};
            note_string=[note_string; ['\chi^2_{Mag. Diags.}=' num2str(odata.chisq_DIAGNOmag_diag_(filedex(i)),'%7.4e')]];
            note_string=[note_string; ['\chi^2_{Cur. Tor.}=' num2str(odata.chisq_Totaltor_current(filedex(i)),'%7.4e')]];
            note_string=[note_string; ['<\beta>=' num2str(odata.avgbeta(filedex(i)),'%6.4f')]];
            text_handle=annotation('textbox','String',note_string,'Color','white',...
                'Units','pixels','Position',text_pos,'FontSize',20);
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend('','Total','Mag. Diags.','Total Toroidal Current','Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',20);
            ylim([0.01 max(odata.chi_total(filedex))]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    case 'opt'
        for i=1:nfiles
            % Remove the text_handle if it exists
            if ishandle(text_handle), delete(text_handle); end
            % Plot VMEC Surfaces
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                plot(r(1,1,j),z(1,1,j),'+w','MarkerSize',10);
                hold on
                for k=2:data.ns-1
                    plot(r(k,:,j),z(k,:,j),'w');
                end
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                % Add Magnetic Axis indicator
                plot(r(1,1,j),-1.45,'v','Color','yellow','MarkerSize',10,'MarkerFaceColor','yellow');
                text(r(1,1,j),-1.45,num2str(r(1,1,j),'%4.3f'),'Color','white','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','bottom');
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([vmec_ymin vmec_ymax]);
                xlim([vmec_xmin vmec_xmax]);
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
                % Ghost Limiter
                %hold on
                %plot(3.664+1.07.*cos(theta),1.07.*sin(theta),'--w','Linewidth',4);
                %hold off
                % Add text to each plot
                if j==1, text(2.55,1.5,'0 Field Period','Color','white',...
                        'FontSize',20); end
                if j==2, text(2.55,1.5,'1/4 Field Period','Color','white',...
                        'FontSize',20); end
                if j==3, text(2.55,1.5,'1/2 Field Period','Color','white',...
                        'FontSize',20); end
            end
            % Plot Pressure and current
            input_data=read_vmec_input(input_filename{i});
            switch axistype
                case 'realspace'
                    %pressure=polyval(fliplr(input_data.am),0:1./(data.ns-1):1);
                    %current=polyval(fliplr(input_data.ac),0:1./(data.ns-1):1);
                    pressure=data.presf;
                    current=data.jcurv;
                    pressure=[fliplr(pressure) pressure(2:numel(pressure))];
                    current=[fliplr(current) current(2:numel(current))];
                    r_real=[r(data.ns:-1:2,181,3)' r(:,1,3)'];
                    % Pressure
                    subplot(4,3,9);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    plot(r_real,pressure,'Color','white')
                    hold on
                    plot([r_real(data.ns) r_real(data.ns)],[min(pressure) max(pressure)],...
                        'Color','yellow')
                    hold off
                    ylabel('Pressure [Pa]');
                    axis tight
                    xlim([vmec_xmin vmec_xmax]);
                    % Current
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    plot(r_real,current/1000.,'Color','white')
                    hold on
                    plot([r_real(data.ns) r_real(data.ns)],[min(current/1000.) max(current/1000.)],...
                        'Color','yellow')
                    hold off
                    ylabel('Itor [kA]');
                    xlabel(['R [m] (phi=18°)' ]);
                    axis tight
                    xlim([vmec_xmin vmec_xmax]);
                case 'fluxspace'
                    %pressure=polyval(fliplr(input_data.am),0:.02:1);
                    %current=polyval(fliplr(input_data.ac),0:.02:1);
                    pressure=data.presf;
                    current=data.jcurv;
                    % Pressure
                    subplot(4,3,9);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    plot(0:1/(data.ns-1):1,pressure,'Color','white')
                    ylabel('Pressure [Pa]');
                    axis tight
                    % Current
                    subplot(4,3,12);
                    set(gca,'nextplot','replacechildren','XColor','white',...
                        'YColor','white','ZColor','white','Color','black','FontSize',20);
                    plot(0:1/(data.ns-1):1,current/1000.,'Color','white')
                    ylabel('Itor [kA]');
                    xlabel('Normalized Flux');
                    axis tight
            end
            % Plot Output
            odata=read_stellopt(o_filename{1});
            subplot(4,3,[7 8 10 11]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot(odata.func_eval(filedex),odata.chi_total(filedex),'o','Color','white',...
                'MarkerSize',10,'LineWidth',2.0);
            hold on
            for j=1:odata.numels
                name=strcat('chisq_',odata.struc_names{j});
                plot(odata.func_eval(filedex),odata.(name)(filedex),'.','Color',colors(j,:),'MarkerSize',10,'LineWidth',2.0);
            end
            plot([odata.func_eval(filedex(i)) odata.func_eval(filedex(i))],[0.1 max(odata.chi_total)],'yellow','LineWidth',2.0);
            hold off
            % Create Chi String
            note_string={['\chi^2_{Total}=' num2str(odata.chi_total(filedex(i)),'%7.4e')]};
            for j=1:odata.numels
                note_string=[note_string; [strcat('\chi^2_{',odata.struc_names{j},'}') num2str(odata.(strcat('chisq_',odata.struc_names{j}))(filedex(j)),'%7.4e')]];
            end
            text_handle=annotation('textbox','String',note_string,'Color','white',...
                'Units','pixels','Position',text_pos,'FontSize',20);
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend(names,'Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',12);
            ylim([0.1 max(odata.chi_total(filedex))]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    case 'recon'
        for i=1:nfiles
            % Get File Number
            dex=strfind(filename{i},'.');
            dex=dex(numel(dex))+1;
            if ~strcmp(filename{i}(dex:numel(filename{i})),'min')
                filedex=str2num(filename{i}(dex:numel(filename{i})));
            %else
            %    filedex=-1;
            end
            % Plot VMEC Surfaces
            theta=0:2*pi/360:2*pi;
            data=read_vmec(filename{i});
            zeta=[0 2*pi/data.nfp/4. 2*pi/data.nfp/2.];
            r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
            z=sfunct(theta,zeta,data.zmns,data.xm,data.xn);
            for j=1:3
                subplot(4,3,[j j+3]);
                set(gca,'nextplot','replacechildren','XColor','white',...
                    'YColor','white','ZColor','white','Color','black','FontSize',20);
                plot(r(1,1,j),z(1,1,j),'+w');
                hold on
                for k=2:data.ns-1
                    plot(r(k,:,j),z(k,:,j),'w');
                end
                plot(r(data.ns,:,j),z(data.ns,:,j),'r','LineWidth',4);
                if plotves
                    hves=plot_vessel(ves_data,'phi',zeta(j));
                    set(hves,'Color','white','LineWidth',4);
                end
                hold off
                ylabel('Z [m]');
                xlabel('R [m]');
                ylim([-1.5 1.5]);
                xlim([2.5 5]);
                if j==2, title(filename{i},'Color','white','Interpreter','none');end
                axis square;
            end
            % Plot Pressure profile
            subplot(4,3,7);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            pdata=read_stellopt(p_filename{i});
            %Find VMEC DOMAIN
            dex1=find(pdata.data(:,4)<1,1,'first');
            dex2=find(pdata.data(:,4)<1,1,'last');
            x1=pdata.data(dex1,1);
            x2=pdata.data(dex2,1);
            y1=max(pdata.data(:,5))*.01;
            switch axistype
                case 'realspace'
                    vmec_line=zeros(1,dex2-dex1+1);
                    hp1=plot(pdata.data(:,1),pdata.data(:,5),'ow');
                    hold on
                    hpatch=fill([x1 x2 x2 x1],[0 0 y1 y1],'red');
                    set(hpatch,'FaceColor','red','EdgeColor','red');
                    hp2=plot(pdata.data(:,1),pdata.data(:,6),'w');
                    text((x2-x1)/3+x1,y1*5,'VMEC Domain','Color','white');
                    hp3=plot([r(1,1,3) r(1,1,3)],[y1 y1].*100,'--w');
                    hold off
                case 'fluxspace'
                    hp1=plot(pdata.data(:,4),pdata(:,5),'ow');
                    hold on
                    hp2=plot(pdata.data(:,4),pdata.data(:,6),'w');
                    hold off
            end
            %title('Pressure Profile','Color','white');
            axis tight
            ylabel('Pressure [Pa]');
            ylim([0 max([max(pdata.data(:,5)) max(pdata.data(:,6))])]);
            hleg=legend([hp1,hp2],'Experiment','Simulation','boxoff');
            set(hleg,'TextColor','white','FontSize',12);
            % Plot Current Density
            subplot(4,3,10);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            input_data=makemovie_read_input(input_filename{i});
            current=polyval(fliplr(input_data.ac),pdata.data(:,4)).*data.Itor;
            current(pdata.data(:,4)>1.0)=0.0;
            switch axistype
                case 'realspace'
                    plot(pdata.data(:,1),current/1000.,'Color','white')
                    xlabel('R[m]');
                case 'fluxspace'
                    plot(pdata.data(:,4),current/1000.,'Color','white')
                    xlable('Normalized Flux');
            end
            title('VMEC Current Density','Color','white')
            ylabel('Itor [kA]');
            axis tight
            % Plot Output
            odata=read_stellopt(o_filename{1});
            subplot(4,3,[8 9 11 12]);
            set(gca,'nextplot','replacechildren','XColor','white',...
                'YColor','white','ZColor','white','Color','black','FontSize',20);
            colormap('lines');
            colors=get(gca,'ColorOrder');
            names=['Total';odata.names(:)];
            plot(odata.func_eval,odata.chi_total,'o','Color','white',...
                'MarkerSize',10);
            hold on
            plot(odata.func_eval,odata.chisq_PressureProfile,'.r');
            plot(odata.func_eval,odata.chisq_PressureProfile./odata.nw_PressureProfile,'xr');
            plot(odata.func_eval,odata.chisq_EdgePressure,'.b');
            plot(odata.func_eval,odata.chisq_PressureGradient,'.g');
            %for j=1:odata.numels
            %    name=strcat('chisq_',odata.struc_names{j});
            %    plot(odata.func_eval,odata.(name),'.','Color',colors(j,:));
            %end
            plot(odata.func_eval,odata.chisq_PressureProfile./odata.nw_PressureProfile,'x','Color','white');
            if (filedex==-1 && isfield(odata,'min_func'))
                filedex=odata.min_func;
            end
            plot([filedex filedex],[0.1 max(odata.chi_total)],'yellow');
            hold off
            set(gca,'Yscale','log');
            ylabel('\chi^2');
            xlabel('Function Evaluations');
            hleg2=legend(names,'Location','NorthEastOutside');
            set(hleg2,'TextColor','white','FontSize',12);
            ylim([0.1 max(odata.chi_total)]);
            pause(0.001);
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
end
% Close Movie File
close(mov);
% Make a poster frame
saveas(gca,[movname,'_poster.fig'],'fig');
hold on
text(0,1,'DONE','Color','yellow','FontSize',50)
hold off
warning('on');
end
