function [ figs, n_fida ] = plot_fidasim(filename,varargin)
%PLOT_FIDASIM Makes plots of FIDASIM data generated with the BEAMS3D
%Interface.
%The PLOT_FIDASIM function creates various canned plot of the data.  The
%
% Example usage
%     filename_b3d = fidasim_38581_4350_b3d_distribution.h5
%     filename_transp = 38581A03-04350_distribution.h5
%     [figs, n_b3d] = plot_fidasim(filename_b3d(1:end-16),'energy','pitch', 'profiles');
%     [figs, n_transp] = plot_fidasim(filename_transp(1:end-16),'energy','pitch','profiles','figs',figs);
%     fac = n_transp/n_b3d;
%     plot_fidasim(filename_b3d(1:end-16),'energy','pitch', 'profiles', 'fac', fac, 'save','figs',figs);
%     plot_fidasim(filename_b3d(1:end-16),'denf2d','ep2d');
%     plot_fidasim(filename_transp(1:end-16),'denf2d','ep2d');
%
% Maintained by: David Kulla (david.kulla@ipp.mpg.de)
% Version:       1.00

ec=1.6021773300E-19; % Charge of an electron (leave alone)

if (strcmp(filename(end-1:end),'h5'))
    disp('ERROR: Only give runid (dist and eq files are loaded automatically!');
    disp(['       Filename: ' filename]);
    data=[];
end

dist_name = [filename, '_distribution.h5'];
eq_name = [filename,'_equilibrium.h5'];
neut_name = [filename, '_neutrals.h5'];
spec_name = [filename,'_spectra.h5'];
geom_name = [filename,'_geometry.h5'];

lsave = 0;
plot_type = {};
figs = {};
fac = 1;
name = filename;
leq = 0;
ldist = 0;
lspec = 0;
lneut = 0;
lgeom =0;
if nargin > 1
    i = 1;
    while i < nargin
        switch varargin{i}
            case {'overview','profiles',...
                    'ba','br2d','bt2d','bz2d',...
                    'brtor','bttor','bztor'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                leq = 1;
            case {'fslice','denf2d','denf','fdenf','fdenf2d',...
                    'pitch','energy','ep2d',...
                    'fdenftor','denftor'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                ldist = 1;
            case{'ndensvert', 'ndenshorz', 'ndenscross'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lneut = 1;
            case{'spectrum'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lspec = 1;
                lgeom = 1;
                i=i+1;
                channel = varargin{i};
            case {'fida', 'bes', 'fidabes'}
                i=i+1;
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lspec = 1;
                lgeom = 1;
            case 'save'
                lsave = 1;
            case 'figs'
                i=i+1;
                figs = varargin{i};
            case 'fac'
                i = i+1;
                fac = varargin{i};
            case 'name'
                i = i+1;
                name = varargin{i};
            otherwise
                disp(['ERROR: Option ', varargin{i}, ' not found!']);
        end
        i=i+1;
    end
end

if ldist
    dist = read_hdf5(dist_name);
    if ~isstruct(dist)
        disp('ERROR: Distribution file not found, check filename!');
        disp(['       Filename: ' filename]);
        data=[];
    end

    dr = dist.r(2)-dist.r(1);
    dz = dist.z(2)-dist.z(1);
    nr = double(dist.nr);
    nz = double(dist.nz);
    if ndims(dist.f) == 5
        dphi = dist.phi(2) - dist.phi(1);
        % Area
        %area=dr*100.*dz*100;

        nphi=double(dist.nphi);
        area=dr.*dz;
        % Volume (function of R)
        vol = dist.r.*dphi.*area;
        vol2d=repmat(vol,[1 dist.nz dist.nphi]);
        disp(['Total vol2d  : ', num2str(sum(vol2d,'all'))])

        tmp = pi*((dist.r(end)+dr).^2-dist.r(1).^2)*(dist.z(end)-dist.z(1));
        disp(['Total cyl vol: ', num2str(sum(tmp,'all'))]);
        %Correct for trapz(ones(4,1)=3 (edges not correctly accounted for)
        z = trapz(dz*nz/(nz-1),ones(nz,1));
        rdphi = repmat(dist.r.*z,1,dist.nphi);
        tmp = squeeze(trapz(dr*nr/(nr-1),trapz(dphi*nphi/(nphi-1),rdphi,2)));

        disp(['Total int vol: ', num2str(sum(tmp,'all'))]);
        n_fida = sum(dist.denf.*vol2d,'all');
    else
        n_fida = 2*pi*dr*dz*sum(dist.r.*sum(squeeze(dist.denf(:,:,1)),2)); %Axisymmetric only.
    end

    [~,z0_ind]=min(abs(dist.z));
end
if leq
    eq = read_hdf5(eq_name);
    if ~isstruct(eq)
        disp('ERROR: Equilbirium file not found, check filename!');
        disp(['       Filename: ' filename]);
        data=[];
    end
end

if lneut
    neut = read_hdf5(neut_name);
    if ~isstruct(neut)
        disp('ERROR: Neutrals file not found, check filename!');
        disp(['       Filename: ' filename]);
        data=[];
    end
end

if lspec
    spec = read_hdf5(spec_name);
    if ~isstruct(spec)
        disp('ERROR: Spectra file not found, check filename!');
        disp(['       Filename: ' filename]);
        data=[];
    end
end
if lgeom
    geom = read_hdf5(geom_name);
    if ~isstruct(geom)
        disp('ERROR: Geometry file not found, check filename!');
        disp(['       Filename: ' filename]);
        data=[];
    end
end



for i = 1:size(plot_type,2)
    if i>numel(figs)
        figs{i}=figure;
        ax = gca;
        hold on;
    else
        allAxesInFigure = findall(figs{i},'type','axes');
        ax = allAxesInFigure(~ismember(get(allAxesInFigure,'Tag'),{'legend','Colobar'}));
    end
    %figure('Color','white','Position',[1 -100 1024 768])
    switch lower(plot_type{i})
        case 'energy'
            if ndims(dist.f) == 5
                tmp = squeeze(trapz(dz*nz/(nz-1),trapz(dist.pitch,dist.f,2),4));
                tmp =squeeze(trapz(dr*nr/(nr-1),trapz(dphi*nphi/(nphi-1),repmat(dist.r',size(dist.f,1),1,size(dist.f,5)).*tmp,3),2));
                %                                h = repmat(vol2d,1,1,1,size(dist.f,1));
                %                 h = permute(h,[4, 1, 2,3]);
                %                 tmp  = sum(squeeze(trapz(dist.pitch,dist.f,2)).*h,[2,3,4]);

            else
                tmp = squeeze(trapz(dz*nz/(nz-1),trapz(dist.pitch,dist.f,2),4))*2*pi;
                tmp=squeeze(trapz(dr*nr/(nr-1),repmat(dist.r',size(dist.f,1),1).*tmp,2));
            end
            if fac == 1
                plot(ax,dist.energy, tmp,'DisplayName',['Energy - ' name] );
            else
                plot(ax,dist.energy, fac*tmp,'DisplayName',['Energy - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('Energy [keV]')
            ylabel('Fast Ion Distribution [1/keV]')
        case 'pitch'
            if ndims(dist.f) == 5
                tmp = squeeze(trapz(dphi*nphi/(nphi-1),trapz(dz*nz/(nz-1),trapz(dist.energy,dist.f,1),4),5));
                tmp = squeeze(trapz(dr*nr/(nr-1),repmat(dist.r',size(dist.f,2),1).*tmp,2));
                %                 h = repmat(vol2d,1,1,1,size(dist.f,2));
                %                 h = permute(h,[4, 1, 2,3]);
                %                 tmp  = sum(squeeze(trapz(dist.energy,dist.f,1)).*h,[2,3,4]);
            else
                tmp = squeeze(trapz(dz*nz/(nz-1),trapz(dist.energy,dist.f,1),4))*2*pi;
                tmp = squeeze(trapz(dr*nr/(nr-1),repmat(dist.r',size(dist.f,2),1).*tmp,2));
            end
            fprintf('Total fast ions in %s: %3.2e\n',filename,n_fida);
            if fac == 1
                plot(ax,dist.pitch, tmp,'DisplayName',['Pitch - ' name] );
                fprintf('Total from pitch: %3.2e\n',trapz(dist.pitch, tmp));
            else
                plot(ax,dist.pitch, fac*tmp,'DisplayName',['Pitch - ' name ', scaling factor: ' num2str(fac)]);
            end
            %plot(dist.pitch, squeeze(trapz(dist.phi,trapz(dist.z,trapz(dist.r,trapz(dist.energy,dist.f,1),3),4),5)),'DisplayName','Pitch');
            xlabel('Pitch [-]')
            ylabel('Fast Ion Distribution [-]')
        case 'ep2d'
            if ndims(dist.f) == 5
                tmp = squeeze(trapz(dist.phi,trapz(dist.z,dist.f,4),5));
                rtmp = permute(repmat(dist.r,1,size(dist.f,1),size(dist.f,2),1),[2,3,1]);
            else
                tmp = squeeze(trapz(dist.z,dist.f,4))*2*pi;
                rtmp = permute(repmat(dist.r,1,size(dist.f,1),size(dist.f,2),1),[2,3,1]);
            end
            tmp = squeeze(trapz(dist.r,rtmp.*tmp,3));
            pixplot(dist.energy,dist.pitch,tmp)
            ylabel('Pitch [-]')
            xlabel('Energy [keV]')
            cstring='Fast Ion Distribution [1/keV]';
            c = colorbar;
            c.Label.String = cstring;
            if lsave
                sname = [filename, '_', plot_type{i}];
                savefig(sname)
                exportgraphics(gcf,[sname,'.png'],'Resolution',300);
            end
            return
        case 'profiles'
            yyaxis left
            plot(eq.plasma.r, squeeze(eq.plasma.te(:,z0_ind,1)), 'DisplayName','T_e');
            hold on
            plot(eq.plasma.r, squeeze(eq.plasma.ti(:,z0_ind,1)), 'DisplayName','T_i');
            plot(eq.plasma.r, squeeze(eq.plasma.zeff(:,z0_ind,1)), 'DisplayName','Zeff [-]');
            ylabel('T [keV]')
            legend('Location','best')
            yyaxis right
            plot(eq.plasma.r, squeeze(eq.plasma.dene(:,z0_ind,1)), 'DisplayName','n_e');
            xlabel('R [cm]')
            ylabel('n_e [cm^{-3}]')
        case 'ba'
            plot(ax,eq.plasma.r, squeeze(eq.fields.br(:,z0_ind,1)), 'DisplayName','B_r');
            hold on
            plot(ax,eq.plasma.r, squeeze(eq.fields.bt(:,z0_ind,1)), 'DisplayName','B_t');
            plot(ax,eq.plasma.r, squeeze(eq.fields.bz(:,z0_ind,1)), 'DisplayName','B_z');
            xlabel('R [cm]')
            ylabel('Magnetic Field [T]')
            legend()
        case 'fslice'
            [~,e_ind]=min(abs(dist.energy-20));
            [~,p_ind]=min(abs(dist.pitch));
            if fac ==1
                plot(ax,dist.r, squeeze(dist.f(e_ind,p_ind,:,z0_ind,1)),'DisplayName',['f slice - ' name] );
            else
                plot(ax,dist.r, fac*squeeze(dist.f(e_ind,p_ind,:,:,1)),'DisplayName',['f slice - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('R [cm]')
            ylabel('Fast ion distribution slice [1/cm^3/keV/dp]')
            legend()
        case 'denf'
            if fac == 1
                plot(ax,dist.r, squeeze(dist.denf(:,z0_ind,1)), 'DisplayName',['Denf - ' name] );
            else
                plot(ax,dist.r, fac*squeeze(dist.denf(:,z0_ind,1)), 'DisplayName',['Denf - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('R [m]')
            ylabel('Fast ion density [m^{-3}]')
            title('Fast ion density profile at z=0')
        case 'fdenf'
            tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            if fac == 1
                plot(ax,dist.r, tmp(:,z0_ind,1), 'DisplayName',['Denf from f - ' name] );
            else
                plot(ax,dist.r, fac*tmp(:,z0_ind,1), 'DisplayName',['Denf from f - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('R [m]')
            ylabel('Fast ion density [m^{-3}]')
            title('Fast ion density profile at z=0')

        case 'denf2d'
            tmp = dist.denf;
            cstring = 'Fast ion density [m^{-3}]';
        case 'fdenf2d'
            tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            cstring = 'Fast ion density [m^{-3}]';
        case 'br2d'
            tmp = eq.fields.br;
            cstring = 'Magnetic Field B_r [T]';
        case 'bt2d'
            tmp = eq.fields.bt;
            cstring = 'Magnetic Field B_t [T]';
        case 'bz2d'
            tmp = eq.fields.bz;
            cstring = 'Magnetic Field B_z [T]';
        case 'denftor'
            tmp = dist.denf;
            cstring = 'Fast ion density [m^{-3}]';
        case 'fdenftor'
            tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            cstring = 'Fast ion density [m^{-3}]';
        case 'brtor'
            tmp = eq.fields.br;
            cstring = 'Magnetic Field B_r [T]';
        case 'bttor'
            tmp = eq.fields.bt;
            cstring = 'Magnetic Field B_t [T]';
        case 'bztor'
            tmp = eq.fields.bz;
            cstring = 'Magnetic Field B_z [T]';
        case 'ndensvert'
            pixplot(neut.grid.x, neut.grid.z, squeeze(sum(neut.tdens(:,:,20,:) + neut.hdens(:,:,20,:) + neut.fdens(:,:,20,:), 1)))
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
        case 'ndenshorz'
            pixplot(neut.grid.x, neut.grid.y, squeeze(sum(neut.tdens(:,:,:,20) + neut.hdens(:,:,:,20)+ neut.fdens(:,:,:,20), 1)))
            xlabel('Beam Grid X [cm]')
            ylabel('Beam Grid Y [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;

        case 'ndenscross'
            pixplot(neut.grid.y, neut.grid.z, squeeze(sum(neut.tdens(:,40,:,:) + neut.hdens(:,40,:,:)+ neut.fdens(:,40,:,:), 1)))
            xlabel('Beam Grid Y [cm]')
            ylabel('Beam Grid Z [cm]')
            cstring = 'Beam neutral density [1/cm^3]';
            c = colorbar;
            c.Label.String = cstring;
        case 'spectrum'
            %[R_data, bes_data, fida_data,spec_data, lambda_data, names,bes_err,fida_err,fida_bes_err,dispersion_data] = get_bes_fida_aug_data(filename_cfr,time(tim),bes_range,fida_range);
            %[R, bes, fida, spec_sim, lambda_sim] = get_bes_fida(spec_name, bes_range(:,:),fida_range,dispersion_data,lambda_data);

            specr = spec.full + spec.half + spec.third + spec.halo + spec.dcx + spec.fida;% + spec.brems;
            plot(spec.lambda,specr(:,channel), 'DisplayName', 'Spectrum')
            hold on
            plot(spec.lambda, spec.full(:,channel), 'DisplayName','Full Energy')
            plot(spec.lambda, spec.half(:,channel), 'DisplayName','Half Energy')
            plot(spec.lambda, spec.third(:,channel), 'DisplayName','Third Energy')
            plot(spec.lambda, spec.halo(:,channel)+spec.dcx(:,channel), 'DisplayName','Halo+DCX') %+spec.brems(:,channel)
            plot(spec.lambda, spec.fida(:,channel), 'DisplayName','FIDA')
            xlabel('Wavelenght [nm]')
            ylabel(' Intensity [Ph/(s nm m^2 sr)]')
            set(gca,'YScale','log')
            xlim([650 663]);
            ylim([5e15, 1e19]);
            title(['Channel: ' geom.spec.id(channel)])
            legend();
        case 'fidabes'


    end
    if strcmp(plot_type{i}(end-1:end),'2d')
        pixplot(dist.r,dist.z,tmp(:,:,1))
        xlabel('R [cm]')
        ylabel('Z [cm]')
        c = colorbar;
        c.Label.String = cstring;
    elseif strcmp(plot_type{i}(end-2:end),'tor')
        if ndims(dist.f) == 5
            pixplot(dist.r,dist.phi,squeeze(tmp(:,z0_ind,:)))
            xlabel('R [cm]')
            ylabel('Z [cm]')
            c = colorbar;
            c.Label.String = cstring;
        else
            disp('4D Distribution has no toroidal information')
        end
    end
    if lsave
        legend(ax);
        sname = [filename, '_', plot_type{i}];
        savefig(figs{i},sname)
        exportgraphics(figs{i},[sname,'.png'],'Resolution',300);
    end

end

end

function [R, bes, fida, spectmp, lambda] = get_bes_fida(filename,bes_range, fida_range,dispersion,lambda_dat)
%h5info(filename);
R = h5read(filename,'/radius');
full = h5read(filename,'/full');
half= h5read(filename,'/half');
third= h5read(filename,'/third');
halo= h5read(filename,'/halo');
dcx= h5read(filename,'/dcx');
fida= h5read(filename,'/fida');
brems= h5read(filename,'/brems');
lambda = h5read(filename,'/lambda');
spec = full + half + third + halo + dcx + fida;% + brems;



%spect = interp1(lambda, (spec').', (lambda_dat').'); %Vectorized interpolation
for i = 1:size(lambda_dat,2)
    spectmp(:,i) = interp1(lambda, spec(:,i), lambda_dat(:,i),'spline');
end


bes_dex = (lambda_dat > repmat(bes_range(:,1)',size(lambda_dat,1),1)) & (lambda_dat < repmat(bes_range(:,2)',size(lambda_dat,1),1));
fida_dex = (lambda_dat > fida_range(1)) & (lambda_dat < fida_range(2));
spectmp = movmean(spectmp,5);
disp('Applying moving mean to FIDASIM data');
% for i = 1:size(lambda_dat,2)
% dispersion_tmp(:,i) = interp1(lambda_dat(:,i),dispersion(:,i),lambda);
% end
%
% bes = sum(spec.*dispersion_tmp.*bes_dex,1);
% fida = sum(spec.*dispersion_tmp.*fida_dex,1);


bes = sum(spectmp.*dispersion.*bes_dex,1,'omitnan');
fida = sum(spectmp.*dispersion.*fida_dex,1,'omitnan');

end


