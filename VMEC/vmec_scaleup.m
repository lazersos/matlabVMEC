function vmec_out = vmec_scaleup(vmec_data,scale,varargin)
%VMEC_SCALEUP Scales up VMEC geometry
%   The VMEC_SCALEUP routine creates a photographic enlargement of a VMEC
%   equilibria while preserving the major radius of the device.  It takes a
%   VMEC data structure as returned by READ_VMEC and a scale factor as
%   input and returns the VMEC data structure with enlarged surfaces.  Note
%   that no other arrays are touched so the GMNC/S, LMNC/S, and magnetic 
%   field arrays will no longer be correct. Optional arguments are 'plots'
%   to make example plots. 'aspect' Hold aspect ratio fixed.
%
%   Exmaple Usage (data assumed to have at least 10 flux surfaces)
%       data=read_vmec('wout_test.nc');
%       data_big = vmec_scaleup(data,1.5);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version: 1.0

% Defaults
lplot=0;
lfixmajor=1; % Fix

% Handle varargin
if ~isempty(varargin)
    j=1;
    while j<=length(varargin)
        switch varargin{j}
            case{'plots'}
                lplot = 1;
            case{'aspect'}
                lfixmajor=0;
        end
        j=j+1;
    end
end

vmec_out=vmec_data;
dexa = vmec_data.xm==0;
%dexa = abs(vmec_data.rmnc(:,1))>0;

% Remove geometric axis
if lfixmajor
    for i = 1:vmec_out.ns
        vmec_out.rmnc(dexa,i) = vmec_out.rmnc(dexa,i)-vmec_data.rmnc(dexa,1);
        vmec_out.zmns(dexa,i) = vmec_out.zmns(dexa,i)-vmec_data.zmns(dexa,1);
    end
    
    if vmec_out.iasym==1
        for i = 1:vmec_out.ns
            vmec_out.rmns(dexa,i) = vmec_out.rmns(dexa,i)-vmec_data.rmns(dexa,1);
            vmec_out.zmnc(dexa,i) = vmec_out.zmnc(dexa,i)-vmec_data.zmnc(dexa,1);
        end
    end
end

% Scale up
vmec_out.rmnc = vmec_out.rmnc.*scale;
vmec_out.zmns = vmec_out.zmns.*scale;

if vmec_out.iasym==1
    vmec_out.rmns = vmec_out.rmns.*scale;
    vmec_out.zmnc = vmec_out.zmnc.*scale;
end

if lfixmajor
    % Add Geometric axis back in
    for i = 1:vmec_out.ns
        vmec_out.rmnc(dexa,i) = vmec_out.rmnc(dexa,i)+vmec_data.rmnc(dexa,1);
        vmec_out.zmns(dexa,i) = vmec_out.zmns(dexa,i)+vmec_data.zmns(dexa,1);
    end
    
    if vmec_out.iasym==1
        for i = 1:vmec_out.ns
            vmec_out.rmns(dexa,i) = vmec_out.rmns(dexa,i)+vmec_data.rmns(dexa,1);
            vmec_out.zmnc(dexa,i) = vmec_out.zmnc(dexa,i)+vmec_data.zmnc(dexa,1);
        end
    end
end

% Make a plot

if lplot
    ns = vmec_data.ns;
    ns2 = round(ns/2);
    theta = deg2rad(0:(360));
    zeta  = deg2rad([0 0.25 0.5 0.75].*round(360./vmec_data.nfp));
    rv    = cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    zv    = sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    if vmec_data.iasym==1
        rv    = rv+sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
        zv    = zv+cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
    end
        
    ro    = cfunct(theta,zeta,vmec_out.rmnc,vmec_out.xm,vmec_out.xn);
    zo    = sfunct(theta,zeta,vmec_out.zmns,vmec_out.xm,vmec_out.xn);
    if vmec_out.iasym==1
        ro    = ro+sfunct(theta,zeta,vmec_out.rmns,vmec_out.xm,vmec_out.xn);
        zo    = zo+cfunct(theta,zeta,vmec_out.zmnc,vmec_out.xm,vmec_out.xn);
    end
    figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    for i=1:4
        subplot(1,4,i);
        plot(rv(1,1,i),zv(1,1,i),'+k'); hold on;
        plot(squeeze(rv(ns2,:,i)),squeeze(zv(ns2,:,i)),'--k','LineWidth',1);
        plot(squeeze(rv(ns,:,i)),squeeze(zv(ns,:,i)),'k','LineWidth',4);
        plot(ro(1,1,i),zo(1,1,i),'xr'); 
        plot(squeeze(ro(ns2,:,i)),squeeze(zo(ns2,:,i)),'--r','LineWidth',1);
        plot(squeeze(ro(ns,:,i)),squeeze(zo(ns,:,i)),'r','LineWidth',4);
        set(gca,'FontSize',24);
        xlabel('R [m]');
        ylabel('Z [m]');
        axis equal;
        title([' \phi = ' num2str(round(rad2deg(zeta(i))),'%i')]);
    end
    annotation('textbox',[0.05 0.025 0.1 0.1],...
        'String',['Scale = ' num2str(scale,'%5.2f')],...
        'FontSize',24,'LineStyle','none');
end

return;


end

