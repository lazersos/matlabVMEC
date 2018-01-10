function [ output_args ] = stellopt_recon_plots( type,extension )
%STELLOPT_RECON_PLOTS Make a stellopt reconstruction plot
%   Detailed explanation goes here

dex = -1;  % Set to -1 to use last iteration


% Get Data
stel_data=read_stellopt(['stellopt.' strtrim(extension)]);
vmec_data=read_vmec('wout_reset_file.nc');
if (dex==-1)
    opt_data=read_namelist(['input.' strtrim(extension) '_min'],'OPTIMUM');
else
    opt_data=read_namelist(['input.' strtrim(extension)],'OPTIMUM');
end

% Make Plots
switch type
    case{'TE'}
        subplot(2,1,1);
        if (dex <=0), dex=size(stel_data.TE_S,1); end
        s=stel_data.TE_S(dex,:);
        te_target=stel_data.TE_target(dex,:);
        te_sigma=stel_data.TE_sigma(dex,:);
        phi=0:.01:1;
        te_equil=pchip(opt_data.te_aux_s,opt_data.te_aux_f,phi);
        errorbar(s,te_target./1000.,te_sigma./1000.,'ok');
        hold on
        plot(phi,te_equil./1000.,'r');
        xlim([0 1.1]);
        title('Electron Temperature');
        ylabel('T_e [keV]');
        xlabel('Normalized Toroidal Flux');
        subplot(2,1,2);
        error=abs(te_target-stel_data.TE_equil(dex,:))./te_sigma;
        error=error(s<=1.0);
        s_temp=s(s<=1.0);
        dex_r=error > 3.;
        dex_y=(error<=3).*(error > 1.);
        dex_g=error <= 1.0;
        if (max(dex_r) > 0)
            hr=bar(s_temp,error.*dex_r,'hist');
            set(hr,'FaceColor','red');
            hold on
        end
        if (max(dex_y) > 0)
            hy=bar(s_temp,error.*dex_y,'hist');
            set(hy,'FaceColor','yellow');
            hold on
        end
        if (max(dex_g) > 0)
            hg=bar(s_temp,error.*dex_g,'hist');
            set(hg,'FaceColor','green');
            hold on
        end
        ylabel('Error');
        xlabel('Index');
        xlim([0 1.1]);
        ylim([0 5]);
        hold off;
    case{'NE'}
        subplot(2,1,1);
        if (dex <=0), dex=size(stel_data.NE_S,1); end
        s=stel_data.NE_S(dex,:);
        ne_target=stel_data.NE_target(dex,:);
        ne_sigma=stel_data.NE_sigma(dex,:);
        phi=0:.01:1;
        ne_equil=pchip(opt_data.ne_aux_s,opt_data.ne_aux_f,phi);
        errorbar(s,ne_target./1E19,ne_sigma./1E19,'ok');
        hold on
        plot(phi,ne_equil./1E19,'r');
        xlim([0 1.1]);
        title('Electron Density');
        ylabel('N_e x10^{19}[m^{-3}]');
        xlabel('Normalized Toroidal Flux');
        subplot(2,1,2);
        error=abs(ne_target-stel_data.NE_equil(dex,:))./ne_sigma;
        error=error(s<=1.0);
        s_temp=s(s<=1.0);
        dex_r=error > 3.;
        dex_y=(error<=3).*(error > 1.);
        dex_g=error <= 1.0;
        if (max(dex_r) > 0)
            hr=bar(s_temp,error.*dex_r,'hist');
            set(hr,'FaceColor','red');
            hold on
        end
        if (max(dex_y) > 0)
            hy=bar(s_temp,error.*dex_y,'hist');
            set(hy,'FaceColor','yellow');
            hold on
        end
        if (max(dex_g) > 0)
            hg=bar(s_temp,error.*dex_g,'hist');
            set(hg,'FaceColor','green');
            hold on
        end
        ylabel('Error');
        xlabel('Index');
        xlim([0 1.1]);
        ylim([0 5]);
        hold off;
    case{'TI'}
        subplot(2,1,1);
        if (dex <=0), dex=size(stel_data.TI_S,1); end
        s=stel_data.TI_S(dex,:);
        ti_target=stel_data.TI_target(dex,:);
        ti_sigma=stel_data.TI_sigma(dex,:);
        phi=0:.01:1;
        ti_equil=pchip(opt_data.ti_aux_s,opt_data.ti_aux_f,phi);
        errorbar(s,ti_target./1000.,ti_sigma./1000.,'ok');
        hold on
        plot(phi,ti_equil./1000.,'r');
        xlim([0 1.1]);
        title('Ion Temperature');
        ylabel('T_i [keV]');
        xlabel('Normalized Toroidal Flux');
        subplot(2,1,2);
        error=abs(ti_target-stel_data.TI_equil(dex,:))./ti_sigma;
        error=error(s<=1.0);
        s_temp=s(s<=1.0);
        dex_r=error > 3.;
        dex_y=(error<=3).*(error > 1.);
        dex_g=error <= 1.0;
        if (max(dex_r) > 0)
            hr=bar(s_temp,error.*dex_r,'hist');
            set(hr,'FaceColor','red');
            hold on
        end
        if (max(dex_y) > 0)
            hy=bar(s_temp,error.*dex_y,'hist');
            set(hy,'FaceColor','yellow');
            hold on
        end
        if (max(dex_g) > 0)
            hg=bar(s_temp,error.*dex_g,'hist');
            set(hg,'FaceColor','green');
            hold on
        end
        ylabel('Error');
        xlabel('Index');
        xlim([0 1.1]);
        ylim([0 5]);
        hold off;
    case{'T_total'}
        subplot(2,2,1);
        if (dex <=0), dex=size(stel_data.NE_S,1); end
        s=stel_data.NE_S(dex,:);
        ne_target=stel_data.NE_target(dex,:);
        ne_sigma=stel_data.NE_sigma(dex,:);
        phi=0:.01:1;
        ne_equil=pchip(opt_data.ne_aux_s,opt_data.ne_aux_f,phi);
        errorbar(s,ne_target./1E19,ne_sigma./1E19,'ok');
        hold on
        plot(phi,ne_equil./1E19,'r');
        xlim([0 1.1]);
        title('Electron Density');
        ylabel('N_e x10^{19}[m^{-3}]');
        xlabel('Normalized Toroidal Flux');
        subplot(2,2,2);
        dex=size(stel_data.TE_S,1);
        s=stel_data.TE_S(dex,:);
        te_target=stel_data.TE_target(dex,:);
        te_sigma=stel_data.TE_sigma(dex,:);
        phi=0:.01:1;
        te_equil=pchip(opt_data.te_aux_s,opt_data.te_aux_f,phi);
        errorbar(s,te_target./1000.,te_sigma./1000.,'ok');
        hold on
        plot(phi,te_equil./1000.,'r');
        xlim([0 1.1]);
        title('Electron Temperature');
        ylabel('T_e [keV]');
        xlabel('Normalized Toroidal Flux');
        subplot(2,2,3);
        dex=size(stel_data.TI_S,1);
        s=stel_data.TI_S(dex,:);
        ti_target=stel_data.TI_target(dex,:);
        ti_sigma=stel_data.TI_sigma(dex,:);
        phi=0:.01:1;
        ti_equil=pchip(opt_data.ti_aux_s,opt_data.ti_aux_f,phi);
        errorbar(s,ti_target./1000.,ti_sigma./1000.,'ok');
        hold on
        plot(phi,ti_equil./1000.,'r');
        xlim([0 1.1]);
        title('Ion Temperature');
        ylabel('T_i [keV]');
        xlabel('Normalized Toroidal Flux');
        subplot(2,2,4);
        phi=vmec_data.phi./vmec_data.phi(vmec_data.ns);
        plot(phi,vmec_data.presf./1000.,'r');
        xlim([0 1.1]);
        temp=ylim;
        ylim([0 temp(2)]);
        title('Plasma Pressure');
        ylabel('P [kPa]');
        xlabel('Normalized Toroidal Flux');
    case{'MSE','MSE_iota'}
        subplot(1,2,1);
        phi=vmec_data.phi./vmec_data.phi(vmec_data.ns);
        q=1./vmec_data.iotaf;
        plotyy(phi,vmec_data.iotaf,phi,vmec_data.jcurv);
        subplot(1,2,2);
        if (dex <=0), dex=size(stel_data.MSE_S,1); end
        s=stel_data.MSE_S(dex,:);
        mse_target=stel_data.MSE_target(dex,:);
        mse_sigma=stel_data.MSE_sigma(dex,:);
        mse_equil=stel_data.MSE_equil(dex,:);
        errorbar(s,mse_target,mse_sigma,'ok');
        xlim([0 1.1]);
        hold on
        plot(s,mse_equil,'+r');
    case{'MSE_q'}
        subplot(1,2,1);
        phi=vmec_data.phi./vmec_data.phi(vmec_data.ns);
        q=1./vmec_data.iotaf;
        plot(phi,q);
        title('q-Profile');
        ylabel('q');
        xlabel('Normalized Toroidal Flux');
        subplot(1,2,2);
        if (dex <=0), dex=size(stel_data.MSE_S,1); end
        s=stel_data.MSE_S(dex,:);
        mse_target=stel_data.MSE_target(dex,:);
        mse_sigma=stel_data.MSE_sigma(dex,:);
        mse_equil=stel_data.MSE_equil(dex,:);
        errorbar(s,mse_target,mse_sigma,'ok');
        xlim([0 1.1]);
        hold on
        plot(s,mse_equil,'+r');
        title('Motional Stark Effect');
        xlabel('Normalized Toroidal Flux');
        ylabel('Signal (radians)');
    case{'FLUXLOOPS'}
        subplot(2,1,1);
        if (dex <=0), dex=size(stel_data.FLUXLOOPS_target,1); end
        fluxloop_target=stel_data.FLUXLOOPS_target(dex,:);
        fluxloop_sigma=stel_data.FLUXLOOPS_sigma(dex,:);
        fluxloop_equil=stel_data.FLUXLOOPS_equil(dex,:);
        dex2=fluxloop_sigma<1E10;
        errorbar(fluxloop_target.*dex2,fluxloop_sigma.*dex2,'ok');
        hold on
        plot(fluxloop_equil.*dex2,'or');
        hold off
        title('FluxLoop Reconstruction');
        ylabel('Flux [Wb]');
        subplot(2,1,2);
        error = abs(fluxloop_target-fluxloop_equil)./fluxloop_sigma;
        dex_r=error > 3.;
        dex_y=(error<=3).*(error > 1.);
        dex_g=error <= 1.0;
        if (max(dex_r) > 0)
            hr=bar(error.*dex_r,'hist');
            set(hr,'FaceColor','red');
            hold on
        end
        if (max(dex_y) > 0)
            hy=bar(error.*dex_y,'hist');
            set(hy,'FaceColor','yellow');
            hold on
        end
        if (max(dex_g) > 0)
            hg=bar(error.*dex_g,'hist');
            set(hg,'FaceColor','green');
            hold on
        end
        ylabel('Error');
        xlabel('Index');
        ylim([0 5]);
        hold off;
        
        
        
        
        
end


end

