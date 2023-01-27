function output = VMECcomp( filename,datatype )
%VMECcomp(filename,datatype) Creates plots comparing VMEC equilibria.
%   The VMECcomp function plots allows the user to create comparrision
%   plots of specific quantities between VMEC equilibria.  This is
%   accomplished by passing either a string with wildcards as the first
%   parameter (or a cell array of file names) and a string inicating the
%   type of value to plot as the second parameter.  Available options are:
%       'curtor'        Toroidal Current
%       'extcur'        Vacuum Field Coil Currents
%       'iota'          Rotational Transform profile
%       'ac_aux'        Current spline
%       'am_aux'        Pressure spline
%       'q'             Safety Factor
%       'pressure'      Pressure profile
%       'current'       Toroidal current profile ('jcurv')
%       'jdotb'         <J*B>
%       'omega'         Rotation
%       'iota_press'    Pressure as a function of iota
%       'iota_pprime'   dp/ds as a function of iota
%       'flux'          Flux surfaces at phi=0
%       'fluxpi2'       Flux surfaces at quarter field period
%       'fluxpi'        Flux surfaces at half field period
%       'flux_edge'     VMEC edge and axis only phi=0
%       'magaxis'       3D plot of magnetic axis trajectory
%
% Example usage
%      haxis=VMEComp('wout*','iota'); % Comparrision plot of iota
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.1

if ~iscell(filename)
    file_struct=dir(filename);
    nfiles=length(file_struct);
    filename=cell(1,nfiles);
    vmec_data=cell(1,nfiles);
    for i=1:nfiles
        filename{i}=file_struct(i).name;
        try
            vmec_data{i}=read_vmec(filename{i});
        catch
            vmec_data{i}=[];
        end
    end
else
    if isstr(filename{1})
        nfiles=max(size(filename));
        vmec_data=cell(1,nfiles);
        for i=1:nfiles
            try
                vmec_data{i}=read_vmec(filename{i});
            catch
                vmec_data{i}=[];
            end
        end
    elseif isstruct(filename{1})
        nfiles=max(size(filename));
        vmec_data=cell(1,nfiles);
        for i=1:nfiles
            vmec_data{i}=filename{i};
            filename{i} = strtrim(vmec_data{i}.input_extension);
        end
    end
    
end

% Fix underscores in filesname
for i=1:length(filename)
    filename{i} = strrep(filename{i},'_','\_');
end

switch datatype
    case 'curtor'
        hold on
        curtor = [];
        for i=1:nfiles
            if ~isfield(vmec_data{i},'ctor'), continue; end
            curtor=[curtor; vmec_data{i}.ctor];
        end
        plot(1:nfiles,curtor,'o');
        set(gca,'XTick',1:nfiles);
        set(gca,'XTickLabel',filename);
        try; rotateXLabels(gca,90);end;
        ylabel('Net Toroidal Current');
        output=gca;
    case 'extcur'
        hold on
        extcur = [];
        for i=1:nfiles
            if ~isfield(vmec_data{i},'extcur'), continue; end
            ns=vmec_data{i}.ns;
            extcur=[extcur; vmec_data{i}.extcur];
        end
        bar3(extcur);
        xlabel('Current Group');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Vaccum Field Currents');
        view(3);
        output=gca;
    case 'iota'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            ns=vmec_data{i}.ns;
            iota=vmec_data{i}.iotaf;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),iota,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),iota,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),iota,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Rotational Transform');
        view(3);
        output=gca;
    case 'iota_flip'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            ns=vmec_data{i}.ns;
            iota=-vmec_data{i}.iotaf;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),iota,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),iota,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),iota,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Rotational Transform');
        view(3);
        output=gca;
    case 'omega'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'omega'), continue; end
            ns=vmec_data{i}.ns;
            omega=vmec_data{i}.omega;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),omega,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),omega,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),omega,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Rotational Transform');
        view(3);
        output=gca;
    case 'q'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'itoaf'), continue; end
            ns=vmec_data{i}.ns;
            iota=vmec_data{i}.iotaf;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),1./iota,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),1./iota,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),1./iota,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Safety Factor (q)');
        view(3);
        output=gca;
    case 'pressure'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'presf'), continue; end
            ns=vmec_data{i}.ns;
            presf=vmec_data{i}.presf;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),presf,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),presf,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),presf,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Pressure');
        view(3);
        output=gca;
    case 'jdotb'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jdotb'), continue; end
            ns=vmec_data{i}.ns;
            jdotb=vmec_data{i}.jdotb;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jdotb,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jdotb,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jdotb,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('<J\cdotB>');
        view(3);
        output=gca;
    case 'buco'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'buco'), continue; end
            ns=vmec_data{i}.ns;
            buco=vmec_data{i}.buco;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),buco,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),buco,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),buco,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('<B_u>');
        view(3);
        output=gca;
    case 'bvco'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bvco'), continue; end
            ns=vmec_data{i}.ns;
            bvco=vmec_data{i}.bvco;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),bvco,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),bvco,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),bvco,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('<B_V>');
        view(3);
        output=gca;
    case 'jcurv'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcurv'), continue; end
            ns=vmec_data{i}.ns;
            jcurv=vmec_data{i}.jcurv;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcurv,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcurv,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcurv,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('<J^v>');
        view(3);
        output=gca;
    case 'jcurv_phi'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcurv'), continue; end
            ns=vmec_data{i}.ns;
            jcurv=vmec_data{i}.jcurv;
            phi = vmec_data{i}.phi;
            if (i == 1)
                plot3(phi,i.*ones(1,ns),jcurv,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(phi,i.*ones(1,ns),jcurv,'r','LineWidth',2.0);
            else
                plot3(phi,i.*ones(1,ns),jcurv,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('<J^v>');
        view(3);
        output=gca;
    case 'jcuru'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcuru'), continue; end
            ns=vmec_data{i}.ns;
            jcuru=vmec_data{i}.jcuru;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcuru,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcuru,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcuru,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('<J^u>');
        view(3);
        output=gca;
    case 'current'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'jcurv'), continue; end
            ns=vmec_data{i}.ns;
            jcurv=vmec_data{i}.jcurv.*2.*pi;
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcurv,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcurv,'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),jcurv,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('dI/ds [A]');
        view(3);
        output=gca;
    case 'iota_press'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            ns=vmec_data{i}.ns;
            iotaf=vmec_data{i}.iotaf;
            presf=vmec_data{i}.presf;
            if (i == 1)
                plot3(iotaf,i.*ones(1,ns),presf,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(iotaf,i.*ones(1,ns),presf,'r','LineWidth',2.0);
            else
                plot3(iotaf,i.*ones(1,ns),presf,'k');
            end
        end
        xlabel('Rotational Transform');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Pressure');
        view(3);
        output=gca;
    case 'iota_pprime'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'iotaf'), continue; end
            ns=vmec_data{i}.ns;
            iotaf=vmec_data{i}.iotaf;
            presf=vmec_data{i}.presf;
            presf=gradient(presf,vmec_data{i}.phipf(1));
            if (i == 1)
                plot3(iotaf,i.*ones(1,ns),presf,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(iotaf,i.*ones(1,ns),presf,'r','LineWidth',2.0);
            else
                plot3(iotaf,i.*ones(1,ns),presf,'k');
            end
        end
        xlabel('Rotational Transform');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('dp/dpsi');
        view(3);
        output=gca;
    case 'ac_aux'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'acauxf'), continue; end
            ns = find(vmec_data{i}.acauxs > 0.0,1,'last');
            ac_aux_s = vmec_data{i}.acauxs(1:ns);
            ac_aux_f = vmec_data{i}.acauxf(1:ns);
            if (i == 1)
                plot3(ac_aux_s(1:ns),i.*ones(1,ns),ac_aux_f(1:ns),'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(ac_aux_s(1:ns),i.*ones(1,ns),ac_aux_f(1:ns),'r','LineWidth',2.0);
            else
                plot3(ac_aux_s(1:ns),i.*ones(1,ns),ac_aux_f(1:ns),'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Current Spline Coefficients');
        view(3);
        output=gca;
    case 'am_aux'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'amauxf'), continue; end
            ns = find(vmec_data{i}.amauxs > 0.0,1,'last');
            am_aux_s = vmec_data{i}.amauxs(1:ns);
            am_aux_f = vmec_data{i}.amauxf(1:ns);
            if (i == 1)
                plot3(am_aux_s(1:ns),i.*ones(1,ns),am_aux_f(1:ns),'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(am_aux_s(1:ns),i.*ones(1,ns),am_aux_f(1:ns),'r','LineWidth',2.0);
            else
                plot3(am_aux_s(1:ns),i.*ones(1,ns),am_aux_f(1:ns),'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Pressure Spline Coefficients');
        view(3);
        output=gca;
    case {'flux0','flux'}
        ntheta=90;
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = 'k';
            if i==1, color='b';end
            if i==nfiles, color='r';end
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,0,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,0,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,0,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,0,zmnc,xm,xn);
            end
            for j=10:10:ns
                plot3(r(j,:),i.*ones(1,ntheta),z(j,:),color);
            end
            plot3(r(1,1),i,z(1,1),'+','Color',color);
            plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
        end
        xlabel('R [m]');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
        axis equal
        output=gca;
    case {'flux_edge'}
        ntheta=90;
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = 'k';
            if i==1, color='b';end
            if i==nfiles, color='r';end
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,0,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,0,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,0,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,0,zmnc,xm,xn);
            end
            plot3(r(1,1),i,z(1,1),'+','Color',color);
            plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
        end
        xlabel('R [m]');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
        axis equal
        output=gca;
    case {'flux_edge3'}
        ntheta=90;
        subplot(1,3,1);
        hold on
        zeta = 0;
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = 'k';
            if i==1, color='b';end
            if i==nfiles, color='r';end
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmnc,xm,xn);
            end
            plot3(r(1,1),i,z(1,1),'+','Color',color);
            plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
        end
        xlabel('R [m]');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
        axis equal
        subplot(1,3,2);
        hold on
        zeta = 2*pi/4/vmec_data{i}.nfp;
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = 'k';
            if i==1, color='b';end
            if i==nfiles, color='r';end
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmnc,xm,xn);
            end
            plot3(r(1,1),i,z(1,1),'+','Color',color);
            plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
        end
        xlabel('R [m]');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
        axis equal
        subplot(1,3,3);
        hold on
        zeta = pi/vmec_data{i}.nfp;
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = 'k';
            if i==1, color='b';end
            if i==nfiles, color='r';end
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,zeta,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,zeta,zmnc,xm,xn);
            end
            plot3(r(1,1),i,z(1,1),'+','Color',color);
            plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
        end
        xlabel('R [m]');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
        axis equal
        output=gcf;
    case 'fluxpi2'
        ntheta=90;
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = 'k';
            if i==1, color='b';end
            if i==nfiles, color='r';end
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            nfp=vmec_data{i}.nfp;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,pi/2/nfp,zmnc,xm,xn);
            end
            for j=10:10:ns
                plot3(r(j,:),i.*ones(1,ntheta),z(j,:),color);
            end
            plot3(r(1,1),i,z(1,1),'+','Color',color);
            plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
        end
        xlabel('R [m]');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
        axis equal
        output=gca;
    case 'fluxpi'
        ntheta=90;
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            color = 'k';
            if i==1, color='b';end
            if i==nfiles, color='r';end
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            nfp=vmec_data{i}.nfp;
            r=cfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,rmnc,xm,xn);
            z=sfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,zmns,xm,xn);
            if (vmec_data{i}.iasym)
                rmns=vmec_data{i}.rmns;
                zmnc=vmec_data{i}.zmnc;
                r=r+sfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,rmns,xm,xn);
                z=z+cfunct(0:2*pi/(ntheta-1):2*pi,pi/nfp,zmnc,xm,xn);
            end
            for j=10:10:ns
                plot3(r(j,:),i.*ones(1,ntheta),z(j,:),color);
            end
            plot3(r(1,1),i,z(1,1),'+','Color',color);
            plot3(r(ns,:),i.*ones(1,ntheta),z(ns,:),color);
        end
        xlabel('R [m]');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('Z [m]');
        view(3);
        axis equal
        output=gca;
    case 'magaxis'
        nzeta=360;
        cosph=cos(0:2*pi/(nzeta-1):2*pi);
        sinph=sin(0:2*pi/(nzeta-1):2*pi);
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'rmnc'), continue; end
            nfp=vmec_data{i}.nfp;
            ns=vmec_data{i}.ns;
            rmnc=vmec_data{i}.rmnc;
            zmns=vmec_data{i}.zmns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            r=cfunct(0,0:2*pi/(nzeta-1):2*pi,rmnc,xm,xn);
            z=sfunct(0,0:2*pi/(nzeta-1):2*pi,zmns,xm,xn);
            plot3(squeeze(r(1,1,:)).*cosph',squeeze(r(1,1,:)).*sinph',squeeze(z(1,1,:)),'k');
        end
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Z [m]');
        view(3);
        axis equal
        output=gca;
    case 'g'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'gmnc'), continue; end
            nfp=vmec_data{i}.nfp;
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.gmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.gmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('g');
        view(3);
        output=gca;
    case 'modb'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bmnc'), continue; end
            nfp=vmec_data{i}.nfp;
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.bmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.bmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('|B|');
        view(3);
        output=gca;
    case 'bsupu'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bsupumnc'), continue; end
            nfp=vmec_data{i}.nfp;
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.bsupumnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.bsupumns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('B^U');
        view(3);
        output=gca;
    case 'bsupv'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'bsupvmnc'), continue; end
            nfp=vmec_data{i}.nfp;
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.bsupvmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.bsupvmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            if (i == 1)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'r','LineWidth',2.0);
            else
                plot3(0:1/(ns-1):1,i.*ones(1,ns),abs(f),'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('B^V');
        view(3);
        output=gca;
    case 'special'
        hold on
        for i=1:nfiles
            if ~isfield(vmec_data{i},'nfp'), continue; end
            mu0=pi*4E-7;
            nfp=vmec_data{i}.nfp;
            ns=vmec_data{i}.ns;
            xm=vmec_data{i}.xm;
            xn=vmec_data{i}.xn;
            fmnc = vmec_data{i}.currvmnc;
            f2mnc = vmec_data{i}.gmnc;
            f3mnc = vmec_data{i}.rmnc;
            f = cfunct(0,0,fmnc,xm,xn);
            f2 = cfunct(0,0,f2mnc,xm,xn);
            f3 = cfunct(0,0,f3mnc,xm,xn);
            f3 = (0:1/(vmec_data{i}.ns-1):1)';
            if (vmec_data{i}.iasym)
                fmns = vmec_data{i}.currvmns;
                f = f+sfunct(0,0,fmns,xm,xn);
            end
            if (i == 1)
                plot3(f3,i.*ones(1,ns),abs(f).*mu0,'b','LineWidth',2.0);
            elseif (i == nfiles)
                plot3(f3,i.*ones(1,ns),abs(f).*mu0,'r','LineWidth',2.0);
            else
                plot3(f3,i.*ones(1,ns),abs(f).*mu0,'k');
            end
        end
        xlabel('Normalized Flux');
        set(gca,'YTick',1:nfiles);
        set(gca,'YTickLabel',filename);
        zlabel('currvmnc');
        view(3);
        output=gca;
    otherwise
        disp(['Error: Datatype ' strtrim(datatype) ' is not supported']);
        output = -1;
end

return

end

