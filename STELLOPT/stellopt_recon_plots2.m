function stellopt_recon_plots2
%stellopt_recon_plots2 Makes reconstruction plots based on STELLOPT runs
%   Just call from command line, it takes care of the rest.
%
% Example usage
%      stellopt_recon_plots2;
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

% Defaults
ntheta = 128;

% Load file
files=dir('stellopt.*');
data = read_stellopt(files.name);
ext = files.name(10:end);
files=dir('wout_*.*.nc');
vmec_data=read_vmec(files(end).name);
files=dir('jacobian.*');
if isempty(files)
    jac_data=[];
else
    jac_data=read_stellopt(files(end).name);
end
files=[];
tprof=[];
dprof=[];

% Process Jacobian
if isempty(jac_data)
    sigma_y = [];
    y_fit   = [];
    ljac=0;
else
    Npar = jac_data.n;
    Npnt = jac_data.m;
    Nfit = jac_data.n;
    y_dat = data.TARGETS(end,:);
    y_fit = data.VALS(end,:);
    y_sig = data.SIGMAS(end,:);
    jac   = jac_data.jac';
    chisq = ((data.TARGETS-data.VALS)./data.SIGMAS).^2;
    chisq_tot = sum(chisq,2);
    delta_y=(y_dat-y_fit)./y_sig;
    weights_sq = (Npnt-Nfit+1)./((delta_y'*delta_y) * ones(Npnt,1));
    weights_sq(delta_y' == 0) = 1.0E-18;
    JtWJ       = jac' * ( jac.* ( weights_sq * ones(1,Npar)));
    covar      = inv(JtWJ);
    sigma_p    = abs(sqrt(diag(covar)));
    sigma_y    = zeros(Npnt,1);
    for i = 1:Npnt
        sigma_y(i) = jac(i,:)*covar*jac(i,:)';
    end
    sigma_y    = abs(sqrt(sigma_y).*y_sig');
    sigma_yp   = abs(sqrt(weights_sq+sigma_y).*y_sig');
    ljac=1;
end


% Cycle through options
if isfield(data,'TE_target')
    if isempty(tprof)
        files = dir('tprof.*.*');
        tprof=importdata(files(end).name);
        tprof=tprof.data;
    end
    zeta = mean(data.TE_PHI(1,:));
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    f = pchip(tprof(:,1),tprof(:,3)./1000,vmec_data.phi./vmec_data.phi(end));
    f = repmat(f',[1 ntheta]);
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(1,2,2);
    torocont(r,z,f,1);
    hold on;
    plot(data.TE_R(1,:),data.TE_Z(1,:),'ow','MarkerSize',8,'LineWidth',2);
    hold off;
    colormap jet;
    set(gca,'FontSize',24);
    xlabel('R [m]');
    ylabel('Z [m]');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'T_e [keV]');
    subplot(1,2,1);
    plot(sqrt(data.TE_S(1,:)),data.TE_target(end,:)./1000,'ok','MarkerSize',8,'LineWidth',2);
    hold on;
    plot(sqrt(tprof(:,1)),tprof(:,3)./1000,'b','LineWidth',4);
    %%%%%%%%% Red lines
    if ljac
        IndexC=strfind(jac_data.target_name,'Electron Temperature');
        Index = find(not(cellfun('isempty',IndexC)));
        yup=[]; ydn=[]; n=1;
        for i=Index
            yup(n) = y_fit(i) + 1.96.*real(sigma_y(i));
            ydn(n) = y_fit(i) - 1.96.*real(sigma_y(i));
            n = n +1;
        end
        s=sqrt(data.TE_S(1,:));
        [~,dex]=sort(s,'ascend');
        fill([s(dex) fliplr(s(dex))],[yup(dex) fliplr(ydn(dex))]./1E3,'blue','EdgeColor','none','FaceAlpha',0.33);
    end
    %%%%%%%%%
    hold off;
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('T_e [keV]');
    ylim([0 1.05.*max(tprof(:,3))/1000]);
    xlim([0 1.2]);
    legend('Exp.','Recon.');
    set(gca,'Position',[0.162 0.237 0.303 0.576]);
    annotation('textbox',[0.1 0.85 0.8 0.1],'string','Electron Temperature Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_te_' ext '.fig']);
    saveas(fig,['recon_te_' ext '.png']);
    close(fig);
end


if isfield(data,'NE_target')
    if isempty(tprof)
        files = dir('tprof.*.*');
        tprof=importdata(files(end).name);
        tprof=tprof.data;
    end
    zeta = mean(data.TE_PHI(1,:));
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    f = pchip(tprof(:,1),tprof(:,2)./1E19,vmec_data.phi./vmec_data.phi(end));
    f = repmat(f',[1 ntheta]);
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(1,2,2);
    torocont(r,z,f,1);
    hold on;
    plot(data.NE_R(1,:),data.NE_Z(1,:),'ow','MarkerSize',8,'LineWidth',2);
    hold off;
    colormap jet;
    set(gca,'FontSize',24);
    xlabel('R [m]');
    ylabel('Z [m]');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'n_e x10^{19} [m^{-3}]');
    subplot(1,2,1);
    plot(sqrt(data.NE_S(1,:)),max(tprof(:,2)).*data.NE_target(end,:)./1E19,'ok','MarkerSize',8,'LineWidth',2);
    hold on;
    plot(sqrt(tprof(:,1)),tprof(:,2)./1E19,'b','LineWidth',4);
    %%%%%%%%% Red lines
    if ljac
        IndexC=strfind(jac_data.target_name,'Electron Density');
        Index = find(not(cellfun('isempty',IndexC)));
        yup=[]; ydn=[]; n=1;
        for i=Index
            if ~isempty(strfind(jac_data.target_name{i},'Line Integrated')), continue; end;
            yup(n) = y_fit(i) + 1.96.*real(sigma_y(i));
            ydn(n) = y_fit(i) - 1.96.*real(sigma_y(i));
            n = n +1;
        end
        s=sqrt(data.NE_S(1,:));
        [~,dex]=sort(s,'ascend');
        fill([s(dex) fliplr(s(dex))],max(tprof(:,2)).*[yup(dex) fliplr(ydn(dex))]./1E19,'blue','EdgeColor','none','FaceAlpha',0.33);
    end
    %%%%%%%%%
    hold off;
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('n_e 10^{19} [m^{-3}]');
    ylim([0 1.05.*max(tprof(:,2))/1E19]);
    xlim([0 1.2]);
    legend('Exp.','Recon.');
    set(gca,'Position',[0.162 0.237 0.303 0.576]);
    annotation('textbox',[0.1 0.85 0.8 0.1],'string','Electron Density Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_ne_' ext '.fig']);
    saveas(fig,['recon_ne_' ext '.png']);
    close(fig);
end

if isfield(data,'NELINE_target')
    if isempty(tprof)
        files = dir('tprof.*.*');
        tprof=importdata(files(end).name);
        tprof=tprof.data;
    end
    zeta = mean(data.NELINE_PHI0(1,:));
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    f = pchip(tprof(:,1),tprof(:,2)./1E19,vmec_data.phi./vmec_data.phi(end));
    f = repmat(f',[1 ntheta]);
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(1,2,2);
    torocont(r,z,f,1);
    hold on;
    plot([data.NELINE_R0(1,:); data.NELINE_R1(1,:)],[data.NELINE_Z0(1,:); data.NELINE_Z1(1,:)],'k','LineWidth',2);
    hold off;
    colormap jet;
    set(gca,'FontSize',24);
    xlabel('R [m]');
    ylabel('Z [m]');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'n_e x10^{19} [m^{-3}]');
    subplot(1,2,1);
    plot(1:length(data.NELINE_target(end,:)),data.NELINE_target(end,:)./1E19,'ok','MarkerSize',8,'LineWidth',2);
    hold on;
    plot(1:length(data.NELINE_target(end,:)),data.NELINE_equil(end,:)./1E19,'+b','MarkerSize',8,'LineWidth',2);
    hold off;
    axis tight;
    ylim([0 2.0*max(ylim)]);
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('Line Int. Density [m^{-2}]');
    legend('Exp.','Recon.');
    set(gca,'Position',[0.162 0.237 0.303 0.576]);
    annotation('textbox',[0.1 0.85 0.8 0.1],'string','Line Int Density Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_neline_' ext '.fig']);
    saveas(fig,['recon_neline_' ext '.png']);
    close(fig);
end

if isfield(data,'XICS_BRIGHT_target')
    if isempty(dprof)
        files = dir('dprof.*.*');
        dprof=importdata(files(end).name);
        dprof=dprof.data;
    end
    zeta = mean(data.XICS_BRIGHT_PHI0(1,:));
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    f = pchip(dprof(:,1),dprof(:,2),vmec_data.phi./vmec_data.phi(end));
    f = repmat(f',[1 ntheta]);
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(1,2,2);
    torocont(r,z,f,1);
    hold on;
    plot([data.XICS_BRIGHT_R0(1,:); data.XICS_BRIGHT_R1(1,:)],[data.XICS_BRIGHT_Z0(1,:); data.XICS_BRIGHT_Z1(1,:)],'k','LineWidth',1);
    text(5.15,1,'1');
    text(5.15,0,num2str(length(data.XICS_BRIGHT_R0(1,:))));
    hold off;
    colormap jet;
    set(gca,'FontSize',24);
    xlabel('R [m]');
    ylabel('Z [m]');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'Eff. Emiss. [arb]');
    subplot(1,2,1);
    hold on;
    plot(1:length(data.XICS_BRIGHT_target(end,:)),data.XICS_BRIGHT_target(end,:),'ok','MarkerSize',8,'LineWidth',2);
    plot(1:length(data.XICS_BRIGHT_target(end,:)),data.XICS_BRIGHT_equil(end,:),'+b','MarkerSize',8,'LineWidth',2);
    %%%%%%%%% Red lines
    if ljac
        IndexC=strfind(jac_data.target_name,'XICS Brightness');
        Index = find(not(cellfun('isempty',IndexC)));
        yup=[]; ydn=[]; n=1;
        for i=Index
            yup(n) = y_fit(i) + 1.96.*real(sigma_y(i));
            ydn(n) = y_fit(i) - 1.96.*real(sigma_y(i));
            n = n +1;
        end
        s=1:length(data.XICS_BRIGHT_target(end,:));
        fill([s fliplr(s)],[yup fliplr(ydn)],'blue','EdgeColor','none','FaceAlpha',0.33);
    end
    %%%%%%%%%
    hold off;
    ylim([0 1.2*max(ylim)]);
    xlim([0 length(data.XICS_BRIGHT_target(end,:))+1]);
    set(gca,'FontSize',24);
    xlabel('Channel');
    ylabel('\int\epsilon_{XICS}dl   [arb]');
    legend('Exp.','Recon.');
    set(gca,'Position',[0.162 0.237 0.303 0.576]);
    annotation('textbox',[0.1 0.85 0.8 0.1],'string','XICS Emissivity Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_xics_bright_' ext '.fig']);
    saveas(fig,['recon_xics_bright_' ext '.png']);
    close(fig);
end

if isfield(data,'XICS_target')
    if isempty(tprof)
        files = dir('tprof.*.*');
        tprof=importdata(files(end).name);
        tprof=tprof.data;
    end
    zeta = mean(data.XICS_PHI0(1,:));
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    f = pchip(tprof(:,1),tprof(:,4)./1E3,vmec_data.phi./vmec_data.phi(end));
    f = repmat(f',[1 ntheta]);
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(1,2,2);
    torocont(r,z,f,1);
    hold on;
    plot([data.XICS_R0(1,:); data.XICS_R1(1,:)],[data.XICS_Z0(1,:); data.XICS_Z1(1,:)],'k','LineWidth',1);
    text(5.15,1,'1');
    text(5.15,0,num2str(length(data.XICS_R0(1,:))));
    hold off;
    colormap jet;
    set(gca,'FontSize',24);
    xlabel('R [m]');
    ylabel('Z [m]');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'T_i [keV]');
    subplot(1,2,1);
    hold on;
    plot(1:length(data.XICS_target(end,:)),data.XICS_target(end,:),'ok','MarkerSize',8,'LineWidth',2);
    plot(1:length(data.XICS_target(end,:)),data.XICS_equil(end,:),'+b','MarkerSize',8,'LineWidth',2);
    %%%%%%%%% Red lines
    if ljac
        IndexC=strfind(jac_data.target_name,'XICS Signal');
        Index = find(not(cellfun('isempty',IndexC)));
        yup=[]; ydn=[]; n=1;
        for i=Index
            yup(n) = y_fit(i) + 1.96.*real(sigma_y(i));
            ydn(n) = y_fit(i) - 1.96.*real(sigma_y(i));
            n = n +1;
        end
        s=1:length(data.XICS_target(end,:));
        fill([s fliplr(s)],[yup fliplr(ydn)],'blue','EdgeColor','none','FaceAlpha',0.33);
    end
    %%%%%%%%%
    hold off;
    ylim([0 1.2*max(ylim)]);
    xlim([0 length(data.XICS_target(end,:))+1]);
    set(gca,'FontSize',24);
    xlabel('Channel');
    ylabel('\int\epsilon_{XICS}T_i dl   [arb]');
    legend('Exp.','Recon.');
    set(gca,'Position',[0.162 0.237 0.303 0.576]);
    annotation('textbox',[0.1 0.85 0.8 0.1],'string','XICS Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_xics_ti_' ext '.fig']);
    saveas(fig,['recon_xics_ti_' ext '.png']);
    close(fig);
end

if isfield(data,'XICS_W3_target')
    if isempty(tprof)
        files = dir('tprof.*.*');
        tprof=importdata(files(end).name);
        tprof=tprof.data;
    end
    zeta = mean(data.XICS_W3_PHI0(1,:));
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    f = pchip(tprof(:,1),tprof(:,3)./1E3,vmec_data.phi./vmec_data.phi(end));
    f = repmat(f',[1 ntheta]);
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(1,2,2);
    torocont(r,z,f,1);
    hold on;
    plot([data.XICS_R0(1,:); data.XICS_R1(1,:)],[data.XICS_Z0(1,:); data.XICS_Z1(1,:)],'k','LineWidth',1);
    text(5.15,1,'1');
    text(5.15,0,num2str(length(data.XICS_R0(1,:))));
    hold off;
    colormap jet;
    set(gca,'FontSize',24);
    xlabel('R [m]');
    ylabel('Z [m]');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'T_e [keV]');
    subplot(1,2,1);
    hold on;
    plot(1:length(data.XICS_W3_target(end,:)),data.XICS_W3_target(end,:),'ok','MarkerSize',8,'LineWidth',2);
    plot(1:length(data.XICS_W3_target(end,:)),data.XICS_W3_equil(end,:),'+b','MarkerSize',8,'LineWidth',2);
    %%%%%%%%% Red lines
    if ljac
        IndexC=strfind(jac_data.target_name,'XICS W3 Factor');
        Index = find(not(cellfun('isempty',IndexC)));
        yup=[]; ydn=[]; n=1;
        for i=Index
            yup(n) = y_fit(i) + 1.96.*real(sigma_y(i));
            ydn(n) = y_fit(i) - 1.96.*real(sigma_y(i));
            n = n +1;
        end
        s=1:length(data.XICS_W3_target(end,:));
        fill([s fliplr(s)],[yup fliplr(ydn)],'blue','EdgeColor','none','FaceAlpha',0.33);
    end
    %%%%%%%%%
    hold off;
    ylim([0 1.2*max(ylim)]);
    xlim([0 length(data.XICS_W3_target(end,:))+1]);
    set(gca,'FontSize',24);
    xlabel('Channel');
    ylabel('\int\epsilon_{W3}dl   [arb]');
    legend('Exp.','Recon.');
    set(gca,'Position',[0.162 0.237 0.303 0.576]);
    annotation('textbox',[0.1 0.85 0.8 0.1],'string','XICS W3 Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_xics_te_' ext '.fig']);
    saveas(fig,['recon_xics_te_' ext '.png']);
    close(fig);
end



if isfield(data,'XICS_V_target')
    if isempty(dprof)
        files = dir('dprof.*.*');
        dprof=importdata(files(end).name);
        dprof=dprof.data;
    end
    zeta = mean(data.XICS_V_PHI0(1,:));
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    f = pchip(dprof(:,1),dprof(:,3)./1E3,vmec_data.phi./vmec_data.phi(end));
    f = repmat(f',[1 ntheta]);
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(1,2,2);
    torocont(r,z,f,1);
    hold on;
    plot([data.XICS_R0(1,:); data.XICS_R1(1,:)],[data.XICS_Z0(1,:); data.XICS_Z1(1,:)],'k','LineWidth',1);
    text(5.15,1,'1');
    text(5.15,0,num2str(length(data.XICS_R0(1,:))));
    hold off;
    colormap jet;
    set(gca,'FontSize',24);
    xlabel('R [m]');
    ylabel('Z [m]');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'\phi [kV]');
    subplot(1,2,1);
    hold on;
    plot(1:length(data.XICS_V_target(end,:)),data.XICS_V_target(end,:)./1E3,'ok','MarkerSize',8,'LineWidth',2);
    plot(1:length(data.XICS_V_target(end,:)),data.XICS_V_equil(end,:)./1E3,'+b','MarkerSize',8,'LineWidth',2);
    %%%%%%%%% Red lines
    if ljac
        IndexC=strfind(jac_data.target_name,'XICS Perp. Velocity');
        Index = find(not(cellfun('isempty',IndexC)));
        yup=[]; ydn=[]; n=1;
        for i=Index
            yup(n) = y_fit(i) + 1.96.*real(sigma_y(i));
            ydn(n) = y_fit(i) - 1.96.*real(sigma_y(i));
            n = n +1;
        end
        s=1:length(data.XICS_V_target(end,:));
        fill([s fliplr(s)],[yup fliplr(ydn)]./1E3,'blue','EdgeColor','none','FaceAlpha',0.33);
    end
    %%%%%%%%%
    hold off;
    ylim([1.2*min(ylim) 1.2*max(ylim)]);
    xlim([0 length(data.XICS_V_target(end,:))+1]);
    set(gca,'FontSize',24);
    xlabel('Channel');
    ylabel('\int v dl   [k.Rad]');
    legend('Exp.','Recon.');
    set(gca,'Position',[0.162 0.237 0.303 0.576]);
    annotation('textbox',[0.1 0.85 0.8 0.1],'string','XICS V Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_xics_v_' ext '.fig']);
    saveas(fig,['recon_xics_v_' ext '.png']);
    close(fig);
end

if isfield(data,'ECEREFLECT_target')
    fig = figure('Position',[1 1 1024 768],'Color','white');
    plot(data.ECEREFLECT_freq(end,:),data.ECEREFLECT_target(end,:)./1000,'ok','MarkerSize',8,'LineWidth',2);
    hold on;
    plot(data.ECEREFLECT_freq(end,:),data.ECEREFLECT_equil(end,:)./1000,'+b','MarkerSize',8,'LineWidth',2);
    plot(data.ECEREFLECT_freq(end,:),data.ECEREFLECT_tradx(end,:)./1000,'xb','MarkerSize',8,'LineWidth',2);
    plot(data.ECEREFLECT_freq(end,:),data.ECEREFLECT_trado(end,:)./1000,'^r','MarkerSize',8,'LineWidth',2);
    %%%%%%%%% Red lines
    if ljac
        IndexC=strfind(jac_data.target_name,'ECE Reflectometry Diagnostic');
        Index = find(not(cellfun('isempty',IndexC)));
        yup=[]; ydn=[]; n=1;
        for i=Index
            yup(n) = y_fit(i) + 1.96.*real(sigma_y(i));
            ydn(n) = y_fit(i) - 1.96.*real(sigma_y(i));
            n = n +1;
        end
        s=data.ECEREFLECT_freq(end,:);
        [~,dex]=sort(s,'ascend');
        fill([s(dex) fliplr(s(dex))],[yup(dex) fliplr(ydn(dex))]./1E3,'blue','EdgeColor','none','FaceAlpha',0.33);
    end
    %%%%%%%%%
    hold off;
    set(gca,'FontSize',24);
    xlabel('Freq. [GHz]');
    ylabel('T_{Rad} [keV]');
    axis tight;
    legend('Exp.','Recon.','X-Mode','O-Mode');
    title('ECE Signal Reconstruction');
    saveas(fig,['recon_ece_' ext '.fig']);
    saveas(fig,['recon_ece_' ext '.png']);
    close(fig);
end

if isfield(data,'FLUXLOOPS_target')
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(2,1,1);
    [ax,h1,h2]=plotyy(sqrt(vmec_data.phi./vmec_data.phi(end)),vmec_data.presf./1E3,...
        sqrt(vmec_data.phi./vmec_data.phi(end)),vmec_data.jcurv./1E3);
    set(h1,'Color','b','LineWidth',2);
    set(h2,'Color','r','LineWidth',2);
    set(ax(1),'FontSize',12,'YColor','b');
    set(ax(2),'FontSize',12,'YColor','r');
    xlabel(ax(1),'Rho');
    ylabel(ax(1),'Press. [kPa]');
    ylabel(ax(2),'Current [kA/m^{-2}]');
    subplot(2,1,2);
    %f=data.FLUXLOOPS_equil(end,:)./data.FLUXLOOPS_target(end,:);
    f=(data.FLUXLOOPS_equil(end,:)-data.FLUXLOOPS_target(end,:))./data.FLUXLOOPS_target(end,:);
    f(data.FLUXLOOPS_target(end,:)==0) = 0;
    bar(1:length(data.FLUXLOOPS_target(end,:)),f.*100,'k');
    %plot(1:length(data.FLUXLOOPS_target(end,:)),data.FLUXLOOPS_target(end,:).*1000,'ok','MarkerSize',8,'LineWidth',2);
    %hold on;
    %plot(1:length(data.FLUXLOOPS_equil(end,:)),data.FLUXLOOPS_equil(end,:).*1000,'+b','MarkerSize',8,'LineWidth',2);
    %hold off;
    set(gca,'FontSize',24);
    xlabel('Flux Loop Index');
    %ylabel('Signal [mWb]');
    ylabel('Signal [% Diff.]');
    axis tight;
    %ylim([-20 max(ylim)]);
    %legend('Exp.','Recon.');
    annotation('textbox',[0.1 0.9 0.8 0.1],'string','Flux Loop Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_fluxloop_' ext '.fig']);
    saveas(fig,['recon_fluxloop_' ext '.png']);
    close(fig);
end

if isfield(data,'SEGROG_target')
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(2,1,1);
    [ax,h1,h2]=plotyy(sqrt(vmec_data.phi./vmec_data.phi(end)),vmec_data.presf./1E3,...
        sqrt(vmec_data.phi./vmec_data.phi(end)),vmec_data.jcurv./1E3);
    set(h1,'Color','b','LineWidth',2);
    set(h2,'Color','r','LineWidth',2);
    set(ax(1),'FontSize',12,'YColor','b');
    set(ax(2),'FontSize',12,'YColor','r');
    xlabel(ax(1),'Rho');
    ylabel(ax(1),'Press. [kPa]');
    ylabel(ax(2),'Current [kA/m^{-2}]');
    subplot(2,1,2);
    f=(data.SEGROG_equil(end,:)-data.SEGROG_target(end,:))./data.SEGROG_target(end,:);
    f(data.SEGROG_target(end,:)==0) = 0;
    bar(1:length(data.SEGROG_target(end,:)),f.*100,'k');
    %plot(1:length(data.SEGROG_target(end,:)),data.SEGROG_target(end,:).*1000,'ok','MarkerSize',8,'LineWidth',2);
    %hold on;
    %plot(1:length(data.SEGROG_equil(end,:)),data.SEGROG_equil(end,:).*1000,'+b','MarkerSize',8,'LineWidth',2);
    %hold off;
    set(gca,'FontSize',24);
    xlabel('Rogowski Index');
    %ylabel('Signal [mT-m]');
    ylabel('Signal [% Diff.]');
    axis tight;
    ylim([-100 100]);
    legend('Exp.','Recon.');
    annotation('textbox',[0.1 0.9 0.8 0.1],'string','Rogowski Reconstruction','FontSize',24,'LineStyle','none','HorizontalAlignment','center');
    saveas(fig,['recon_segrog_' ext '.fig']);
    saveas(fig,['recon_segrog_' ext '.png']);
    close(fig);
end

% Plot Profiles
if ~isempty(tprof)
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(2,2,1);
    plot(sqrt(tprof(:,1)),tprof(:,3)./1E3,'b','LineWidth',2);
    hold on;
    plot(sqrt(tprof(:,1)),tprof(:,4)./1E3,'r','LineWidth',2);
    hold off;
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('T [keV]');
    legend('Electron','Ion');
    title('Temperatures');
    subplot(2,2,2);
    plot(sqrt(tprof(:,1)),tprof(:,2)./1E19,'b','LineWidth',2);
    hold on;
    plot(sqrt(tprof(:,1)),tprof(:,5).*tprof(:,2)./1E19,'r','LineWidth',2);
    hold off;
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('n x10^{19} [m^{-3}]');
    legend('Electron','Ion');
    title('Densities');
    subplot(2,2,3);
    plot(sqrt(tprof(:,1)),tprof(:,5),'k','LineWidth',2);
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('Z_{eff}');
    title('Effective Ion Charge');
    subplot(2,2,4);
    plot(sqrt(tprof(:,1)),tprof(:,6)./1E3,'k','LineWidth',2);
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('P [kPA]');
    title('Pressure');
    saveas(fig,['recon_tprof_' ext '.fig']);
    saveas(fig,['recon_tprof_' ext '.png']);
    close(fig);
end


if ~isempty(dprof)
    fig = figure('Position',[1 1 1024 768],'Color','white');
    subplot(2,2,1);
    plot(sqrt(dprof(:,1)),dprof(:,2),'k','LineWidth',2);
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('\epsilon_{XICS} [arb]');
    title('XICS Effective Emissivity');
    ylim([0 1.05*max(ylim)]);
    subplot(2,2,3);
    plot(sqrt(dprof(:,1)),dprof(:,3)./1E3,'k','LineWidth',2);
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('\Phi [kV]');
    title('E-Static Potential');
    subplot(2,2,4);
    plot(sqrt(dprof(:,1)),2.*sqrt(dprof(:,1)).*gradient(dprof(:,3))./(1E3.*vmec_data.Aminor),'k','LineWidth',2);
    set(gca,'FontSize',24);
    xlabel('Rho');
    ylabel('E_\rho [kV/m]');
    title('Rad. Electric Field');
    axis tight;
    saveas(fig,['recon_dprof_' ext '.fig']);
    saveas(fig,['recon_dprof_' ext '.png']);
    close(fig);
end


if ~isempty(vmec_data)
    fig = figure('Position',[1 1 1024 768],'Color','white');
    zeta = [0 0.25 0.5].*2*pi/vmec_data.nfp;
    theta = 0:2*pi./(ntheta-1):2*pi;
    r = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    b0=sum(vmec_data.bmnc);
    b0=b0(1);
    subplot(2,2,1);
    plot(r(end,:,3),z(end,:,3),'r','LineWidth',4);
    hold on;
    plot(r(end,:,2),z(end,:,2),'g','LineWidth',4);
    plot(r(end,:,1),z(end,:,1),'b','LineWidth',4);
    plot(r(1,1,3),z(1,1,3),'+r','MarkerSize',12,'LineWidth',4);
    plot(r(1,1,2),z(1,1,2),'+g','MarkerSize',12,'LineWidth',4);
    plot(r(1,1,1),z(1,1,1),'+b','MarkerSize',12,'LineWidth',4);
    hold off;
    axis equal;
    text(0.02*diff(xlim)+min(xlim),0.29*diff(ylim)+min(ylim),['V=' num2str(vmec_data.Volume,'%5.1f') ' m^3'],'FontSize',18);
    text(0.02*diff(xlim)+min(xlim),0.17*diff(ylim)+min(ylim),['\Phi=' num2str(vmec_data.phi(end),'%6.3f') ' Wb'],'FontSize',18);
    text(0.02*diff(xlim)+min(xlim),0.05*diff(ylim)+min(ylim),['B_0=' num2str(b0,'%5.2f') ' T'],'FontSize',18);
    set(gca,'FontSize',18);
    xlabel('R [m]');
    ylabel('Z [m]');
    title('VMEC Equilibrium');
    subplot(2,2,2);
    plot(vmec_data.phi./vmec_data.phi(end),vmec_data.presf./1E3,'k','LineWidth',2);
    text(0.02*diff(xlim)+min(xlim),0.07*diff(ylim)+min(ylim),['\beta=' num2str(vmec_data.betatot*100,'%4.2f') '%'],'FontSize',18);
    set(gca,'FontSize',18);
    xlabel('Norm. Flux (s)');
    ylabel('P [kPa]');
    title('Pressure Profile');
    subplot(2,2,3);
    plot(vmec_data.phi./vmec_data.phi(end),vmec_data.iotaf,'k','LineWidth',2);
    set(gca,'FontSize',18);
    xlabel('Norm. Flux (s)');
    ylabel('\iota');
    title('Rot. Trans. Profile');
    subplot(2,2,4);
    plot(vmec_data.phi./vmec_data.phi(end),vmec_data.jcurv./1E3,'k','LineWidth',2);
    text(0.02*diff(xlim)+min(xlim),0.07*diff(ylim)+min(ylim),['I_{tor}=' num2str(vmec_data.Itor./1E3,'%6.2f') ' kA'],'FontSize',18);
    set(gca,'FontSize',18);
    xlabel('Norm. Flux (s)');
    ylabel('dI/ds [kA/m^{-2}]');
    title('Current Profile');
    
    saveas(fig,['recon_vmec_' ext '.fig']);
    saveas(fig,['recon_vmec_' ext '.png']);
    close(fig);
end

end

