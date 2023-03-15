function [ ax, n_fida ] = plot_fidasim(filename,varargin)
%PLOT_FIDASIM Makes plots of FIDASIM data generated with the BEAMS3D
%Interface.
%The PLOT_FIDASIM function creates various canned plots of the files used
%for running FIDASIM and the outputs. The function reads the necessary
%files automatically when supplied a FIDASIM runid. The function is quite
%flexible and can also be used to compare data from multiple runs and
%across codes. For plotting FIDA/BES profiles, see plot_fidasim_profiles.
%
% Example usage
%      plot_fidasim(runid);
%      plot_fidasim(runid,'overview'); %NOT IMPLEMENTED
%      plot_fidasim(runid,'profiles'); %Kinetic profiles
%      plot_fidasim(runid,'ba'); %All magnetic field components on midplane
%      plot_fidasim(runid,'bx_'); %Magnetic field on RZ plane (x=r,p,z)
%      plot_fidasim(runid,'bx_'); %Magnetic field on midplane (x=r,p,z)
%      plot_fidasim(runid,'spectrum', channel_no); %Spectral components
%      plot_fidasim(runid,'fslice'); %Outboard midplane profile of dist.(f)
%      plot_fidasim(runid,'denf_'); %FI density on RZ plane from denf
%      plot_fidasim(runid,'fdenf_'); %FI density on RZ plane from f
%      plot_fidasim(runid,'energy'); %Energy dist, int. over real space
%      plot_fidasim(runid,'pitch'); %Pitch dist, int. over real space
%      plot_fidasim(runid,'ep_'); %E-p dist.,  int. over real space
%      !!! '_' can be '2d' (RZ plane) or 'tor' (midplane)
%      plot_fidasim(runid,'ndensvert'); %Neutral density, vertical
%      plot_fidasim(runid,'ndenshorz'); %Neutral density, horizontal
%      plot_fidasim(runid,'ax', ax); %Figure handles for sharing plots
%
% Miscellaneous Arguments
%      plot_fidasim(runid,'mean'); %Apply moving mean to spectrum
%      plot_fidasim(runid,'sim_data',sim_data); %Used for dispersion
%      plot_fidasim(runid,'save'); %Export figures (.fig and .png)
%      plot_fidasim(runid,'name', 'test'); %ID Name for legend
%      plot_fidasim(runid,'fac', 1.0); %Scaling factor
%
% Example usage of comparing BEAMS3D and TRANSP results:
%     filename_b3d = fidasim_b3d
%     filename_transp = transp
%     [ax, n_b3d] = plot_fidasim(filename_b3d,'energy','pitch', 'profiles');
%     [ax, n_transp] = plot_fidasim(filename_transp,'energy','pitch','profiles','ax',ax);
%
% Maintained by: David Kulla (david.kulla@ipp.mpg.de)
% Version:       1.00

ec=1.6021773300E-19; % Charge of an electron (leave alone)

if (strcmp(filename(end-1:end),'h5'))
    disp('ERROR: Only give runid (dist and eq files are loaded automatically!');
    disp(['       Filename: ' filename]);
end

dist_name = [filename, '_distribution.h5'];
eq_name = [filename,'_equilibrium.h5'];
neut_name = [filename, '_neutrals.h5'];
spec_name = [filename,'_spectra.h5'];
geom_name = [filename,'_geometry.h5'];

lsave = 0;
plot_type = {};
ax = {};
fac = 1;
name = filename;
leq = 0;
ldist = 0;
lspec = 0;
lneut = 0;
lmean=0;
lgeom =0;
linput=0;
n_fida=-1;
sim_data = {};
channel = 0;
linestyle = '-';
color='k';
if nargin > 1
    i = 1;
    while i < nargin
        switch varargin{i}
            case {'overview','profiles',...
                    'ba','br2d','bt2d','bz2d',...
                    'brtor','bttor','bztor','q2d'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                leq = 1;
            case {'fslice','denf2d','denf','fdenf','fdenf2d',...
                    'pitch','energy','ep2d',...
                    'fdenftor','denftor'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                ldist = 1;
                leq=1;
            case{'ndensvert', 'ndenshorz', 'ndenscross'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lneut = 1;
                lgeom=1;
                linput=1;
            case{'ndens2d', 'ndenstor','ndens'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lneut = 1;
                lgeom=1;
                linput=1;
                leq=1;
            case{'spectrum'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lspec = 1;
                lgeom = 1;
                i=i+1;
                channel = varargin{i};
            case{'los3d', 'lostor','los2d'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                lspec = 1;
                lgeom = 1;
                i=i+1;
                channel = varargin{i};
                i=i+1;
                length=varargin{i};
                %                 i=i+1;
                %                 color=varargin{i};
            case {'fida', 'bes', 'fidabes'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                disp(['ERROR: Option ', varargin{i}, ' not implemented here. Use plot_fidasim_profiles instead.']);
                lspec = 1;
                lgeom = 1;
            case 'mean'
                lmean = 1;
            case 'sim_data'
                i=i+1;
                sim_data = varargin{i};
            case 'save'
                lsave = 1;
            case 'ax'
                i=i+1;
                ax = varargin{i};
            case 'fac'
                i = i+1;
                fac = varargin{i};
            case 'name'
                i = i+1;
                name = varargin{i};
            case 'style'
                i = i+1;
                linestyle = varargin{i};
            otherwise
                disp(['ERROR: Option ', varargin{i}, ' not found!']);
        end
        i=i+1;
    end
end

if linput
    input=read_namelist([filename,'_inputs.dat'],'fidasim_inputs');
end

if ldist
    dist = read_hdf5(dist_name);
    if ~isstruct(dist)
        disp('ERROR: Distribution file not found, check filename!');
        disp(['       Filename: ' filename]);
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
        %disp(['Total vol2d  : ', num2str(sum(vol2d,'all'))])

        tmp = pi*((dist.r(end)+dr).^2-dist.r(1).^2)*(dist.z(end)-dist.z(1));
        %disp(['Total cyl vol: ', num2str(sum(tmp,'all'))]);
        %Correct for trapz(ones(4,1)=3 (edges not correctly accounted for)
        z = trapz(dz*nz/(nz-1),ones(nz,1));
        rdphi = repmat(dist.r.*z,1,dist.nphi);
        tmp = squeeze(trapz(dr*nr/(nr-1),trapz(dphi*nphi/(nphi-1),rdphi,2)));

        %disp(['Total int vol: ', num2str(sum(tmp,'all'))]);
        n_fida = sum(dist.denf.*vol2d,'all');
    else
        n_fida = 2*pi*dr*dz*sum(dist.r.*sum(squeeze(dist.denf(:,:,1)),2)); %Axisymmetric only.
    end

    [~,z0_ind]=min(abs(dist.z-2));
end
if leq
    eq = read_hdf5(eq_name);
    if ~isstruct(eq)
        disp('ERROR: Equilbirium file not found, check filename!');
        disp(['       Filename: ' filename]);
    end
    [~,z0_ind]=min(abs(eq.fields.z-2));
end

if lneut
    neut = read_hdf5(neut_name);
    if ~isstruct(neut)
        disp('ERROR: Neutrals file not found, check filename!');
        disp(['       Filename: ' filename]);
    end
    neut.dens = neut.dcxdens+neut.fdens+neut.hdens+neut.tdens+neut.halodens;
    neut.dens=squeeze(sum(neut.dens,1));%Sum over all levels
    neut.grid.vol = repmat(mean(diff(neut.grid.x))*mean(diff(neut.grid.y))*mean(diff(neut.grid.z)),size(neut.grid.x_grid));
    neut.nparts=neut.dens.*neut.grid.vol;
      nneutrals=1.d6*input.pinj/ (1.d3*input.einj*ec...
         *( input.current_fractions(1)      ...
         +  input.current_fractions(2)/2.d0 ...
         +  input.current_fractions(3)/3.d0 ) );
end

if lspec
    spec = read_hdf5(spec_name);
    if ~isstruct(spec)
        disp('ERROR: Spectra file not found, check filename!');
        disp(['       Filename: ' filename]);
    end
    %[~,I] = sort(spec.radius);

end
if lgeom
    geom = read_hdf5(geom_name);
    if ~isstruct(geom)
        disp('ERROR: Geometry file not found, check filename!');
        disp(['       Filename: ' filename]);
        lgeom=0;
    end
end


if ischar(channel)
    chan_description = channel;
    if ~isempty(sim_data)
        channel = find(strcmp(channel,cellstr(deblank(sim_data.names'))));
        %channel = I(channel);%-1;
    elseif strcmp(channel,'all')
        channel = true(size(deblank(geom.spec.id)));
    else
        channel = contains(deblank(geom.spec.id),channel);
    end
elseif iscell(channel)
    chan_description=channel;
    if ~isempty(sim_data)
        channel_tmp = [];
        for i=1:numel(channel)
            channel_tmp =[channel_tmp, find(strcmp(channel{i},cellstr(deblank(sim_data.names'))))];
        end
        channel=sort(channel_tmp);
        %channel = I(channel);%-1;
    else
        channel_tmp = false(geom.spec.nchan,1);
        for i=1:numel(channel)
            %channel_tmp(:,i) = or(channel_tmp,contains(deblank(geom.spec.id),channel{i}));
            channel_tmp(:,i) = contains(deblank(geom.spec.id),channel{i});
        end
        channel = channel_tmp;
    end

end

for i = 1:size(plot_type,2)
    if i>numel(ax)
        %figs{i}=figure;
        figure;
        ax{i} = gca;
        hold(ax{i},'on');
    else
        %allAxesInFigure = findall(figs{i},'type','axes');
        %ax{i} = allAxesInFigure(~ismember(get(allAxesInFigure,'Tag'),{'legend','Colobar'}));
        hold(ax{i},'on');
    end
    %figure('Color','white','Position',[1 -100 1024 768])
    switch lower(plot_type{i})
        case 'energy'
            if ndims(dist.f) == 5
                tmp = squeeze(trapz(dz*nz/(nz-1),trapz(dist.pitch,dist.f,2),4));
                tmp =squeeze(trapz(dr*nr/(nr-1),trapz(dphi*nphi/(nphi-1),repmat(dist.r',size(dist.f,1),1,size(dist.f,5)).*tmp,3),2));
            else
                tmp = squeeze(trapz(dz*nz/(nz-1),trapz(dist.pitch,dist.f,2),4))*2*pi;
                tmp=squeeze(trapz(dr*nr/(nr-1),repmat(dist.r',size(dist.f,1),1).*tmp,2));
            end
            if fac == 1
                % plot(ax,dist.energy(2:end), tmp(1:end-1),'DisplayName',['Energy - ' name] );
                plot(ax{i},dist.energy, tmp,'DisplayName',['Energy - ' name] );
            else
                plot(ax{i},dist.energy, fac*tmp,'DisplayName',['Energy - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('Energy [keV]')
            ylabel('Fast Ion Distribution [1/keV]')
        case 'pitch'
            if ndims(dist.f) == 5
                tmp = squeeze(trapz(dphi*nphi/(nphi-1),trapz(dz*nz/(nz-1),trapz(dist.energy,dist.f,1),4),5));
                tmp = squeeze(trapz(dr*nr/(nr-1),repmat(dist.r',size(dist.f,2),1).*tmp,2));
            else
                tmp = squeeze(trapz(dz*nz/(nz-1),trapz(dist.energy,dist.f,1),4))*2*pi;
                tmp = squeeze(trapz(dr*nr/(nr-1),repmat(dist.r',size(dist.f,2),1).*tmp,2));
            end
            fprintf('Total fast ions in %s: %3.2e\n',filename,n_fida);
            if fac == 1
                plot(ax{i},dist.pitch, tmp,'DisplayName',['Pitch - ' name] );
                fprintf('Total from pitch: %3.2e\n',trapz(dist.pitch, tmp));
            else
                plot(ax{i},dist.pitch, fac*tmp,'DisplayName',['Pitch - ' name ', scaling factor: ' num2str(fac)]);
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
            xlim([dist.energy(1) dist.energy(end)])
            ylim([dist.pitch(1) dist.pitch(end)])
            if lsave
                sname = [filename, '_', plot_type{i}];
                savefig(sname)
                exportgraphics(gcf,[sname,'.png'],'Resolution',300);
            end
            return
        case 'profiles'
            yyaxis(ax{i},'left')
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.te(:,z0_ind,1)), 'DisplayName',['T_e - ' name] );
            hold on
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.ti(:,z0_ind,1)), 'DisplayName',['T_i - ' name] );
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.zeff(:,z0_ind,1)), 'DisplayName',['Zeff [-] - ' name] );
            ylabel(ax{i},'T [keV]')
            legend(ax{i},'Location','best')
            yyaxis(ax{i},'right')
            plot(ax{i},eq.plasma.r, squeeze(eq.plasma.dene(:,z0_ind,1)), 'DisplayName',['n_e - ' name] );
            xlabel(ax{i},'R [cm]')
            ylabel(ax{i},'n_e [cm^{-3}]')
        case 'ba'
            plot(ax{i},eq.plasma.r, squeeze(eq.fields.br(:,z0_ind,1)),linestyle, 'DisplayName','B_r');
            plot(ax{i},eq.plasma.r, squeeze(eq.fields.bt(:,z0_ind,1)),linestyle, 'DisplayName','B_t');
            plot(ax{i},eq.plasma.r, squeeze(eq.fields.bz(:,z0_ind,1)),linestyle, 'DisplayName','B_z');
            xlabel(ax{i},'R [cm]')
            ylabel(ax{i},'Magnetic Field [T]')
            legend()
        case 'fslice'
            [~,e_ind]=min(abs(dist.energy-20));
            [~,p_ind]=min(abs(dist.pitch));
            if fac ==1
                plot(ax{i},dist.r, squeeze(dist.f(e_ind,p_ind,:,z0_ind,1)),linestyle,'DisplayName',['f slice - ' name] );
            else
                plot(ax{i},dist.r, fac*squeeze(dist.f(e_ind,p_ind,:,:,1)),linestyle,'DisplayName',['f slice - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('R [cm]')
            ylabel('Fast ion distribution slice [1/cm^3/keV/dp]')
            legend()
        case 'denf'
            if fac == 1
                plot(ax{i},dist.r, squeeze(dist.denf(:,z0_ind,1)), 'DisplayName',['Denf - ' name] );
            else
                plot(ax{i},dist.r, fac*squeeze(dist.denf(:,z0_ind,1)), 'DisplayName',['Denf - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('R [m]')
            ylabel('Fast ion density [m^{-3}]')
            %title(sprintf('FI density profile at z= %.2f',dist.z(z0_ind)))
        case 'fdenf'
            tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            if fac == 1
                plot(ax{i},dist.r, tmp(:,z0_ind,1), 'DisplayName',['Denf from f - ' name] );
            else
                plot(ax{i},dist.r, fac*tmp(:,z0_ind,1), 'DisplayName',['Denf from f - ' name ', scaling factor: ' num2str(fac)]);
            end
            xlabel('R [m]')
            ylabel('Fast ion density [m^{-3}]')
            title('Fast ion density profile at z=0')
        case 'denf2d'
            r = dist.r;
            z = dist.z;
            tmp = dist.denf;
            cstring = 'Fast ion density [m^{-3}]';
        case 'fdenf2d'
            r = dist.r;
            z = dist.z;
            tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            cstring = 'Fast ion density [m^{-3}]';
        case 'br2d'
            r = eq.fields.r;
            z = eq.fields.z;
            tmp = eq.fields.br;
            cstring = 'Magnetic Field B_r [T]';
        case 'bt2d'
            r = eq.fields.r;
            z = eq.fields.z;
            tmp = eq.fields.bt;
            cstring = 'Magnetic Field B_t [T]';
        case 'bz2d'
            r = eq.fields.r;
            z = eq.fields.z;
            tmp = eq.fields.bz;
            cstring = 'Magnetic Field B_z [T]';
        case 'q2d'
            r = eq.fields.r;
            z = eq.fields.z;
            tmp = (eq.fields.bt(:,:,1) .* eq.fields.r2d) ./ sqrt(eq.fields.bz(:,:,1).^2 + eq.fields.br(:,:,1).^2);
            cstring = 'Magnetic Field B_z [T]';
        case 'denftor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            tmp = dist.denf;
            cstring = 'Fast ion density [m^{-3}]';
        case 'fdenftor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            cstring = 'Fast ion density [m^{-3}]';
        case 'brtor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            tmp = eq.fields.br;
            cstring = 'Magnetic Field B_r [T]';
        case 'bttor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            tmp = eq.fields.bt;
            cstring = 'Magnetic Field B_t [T]';
        case 'bztor'
            r = eq.fields.r;
            phi = eq.fields.phi;
            tmp = eq.fields.bz;
            cstring = 'Magnetic Field B_z [T]';
        case 'ndens'
            neut_r = sqrt(neut.grid.x_grid(:).^2+neut.grid.y_grid(:).^2);
            [discr,redges]=discretize(neut_r,128);
            tmp = accumarray(discr,neut.nparts(:));
            plot(redges(1:end-1),tmp./diff(redges));
            xlabel('R [m]')
            ylabel('Deposition [particles/s]')
            title('Radial Deposition profile')
        case 'ndens2d'
            %eq.fields.phi=linspace(0,1,32);
            %[r,z,phi] = ndgrid(eq.fields.r,eq.fields.z,eq.fields.phi);
            %uvw = [r(:).*cos(phi(:)), -r(:).*sin(phi(:)),z(:)];
            %uvw = [u(:),v(:),w(:)];
            %[xyz] = uvw_to_xyz(input.alpha, input.beta, input.gamma, uvw, input.origin);
            %ngrid = [neut.grid.x_grid(:),neut.grid.y_grid(:),neut.grid.z_grid(:)];
            neut_r = sqrt(neut.grid.x_grid(:).^2+neut.grid.y_grid(:).^2);
            neut_phi= atan2(neut.grid.y_grid(:),neut.grid.x_grid(:));
            neut_z = neut.grid.z_grid(:);
            [discphi,phiedges]=discretize(neut_phi,8);
            [discr,redges]=discretize(neut_r,64);
            [discz,zedges]=discretize(neut_z,32);
            r = redges(1:end-1)+mean(diff(redges))/2;
            z = zedges(1:end-1)+mean(diff(zedges))/2;
            tmp = accumarray([discr,discz,discphi],neut_r.*neut.dens(:));
            tmp=sum(tmp,3)*(phiedges(2)-phiedges(1));%Sum over phi
            %             ndens_F = scatteredInterpolant(ngrid,neut.dens(:));
            %             ndens_F.Method = 'linear';
            %             ndens_F.ExtrapolationMethod = 'none';
            %             ndens = ndens_F(uvw);
            %             tmp=reshape(ndens,size(r));
            cstring='Neutral Density [neutrals/cm^3]';

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
        case 'fida'

        case 'spectrum'
            specr = spec.full + spec.half + spec.third + spec.halo + spec.dcx + spec.fida;% + spec.brems;
            specr = specr.*fac;
            if lmean
                k = 12;
                specr = movmean(specr,k);
                disp(['Applying moving mean with length ', num2str(k), ' to FIDASIM data: ', filename]);
            end
            if ~isempty(sim_data)
                %                 for j = 1:size(sim_data.lambda,2)
                %                     spectmp(:,j) = interp1(spec.lambda, specr(:,j), sim_data.lambda(:,j),'pchip');
                %                 end
                %                 disp('Interpolated wavelength to match data');
                if fac~=1.0
                    name = [name,', scale=',num2str(fac)];
                end
                cwav_mid=mean(spec.lambda);
                instfu = box_gauss_funct(spec.lambda,0.,1.,cwav_mid,sim_data.instfu_gamma,sim_data.instfu_box_nm);
                plot(spec.lambda,conv(specr(:,channel),instfu(:,channel),'same'), 'DisplayName', ['Spectrum - ' name] );
                %                 plot(spec.lambda, conv(spec.full(:,channel),instfu(:,channel),'same'), 'DisplayName',['Full - ' name] );
                %                 plot(spec.lambda, conv(spec.half(:,channel),instfu(:,channel),'same'),  'DisplayName',['Half - ' name] );
                %                 plot(spec.lambda, conv(spec.third(:,channel),instfu(:,channel),'same'),  'DisplayName',['Third - ' name] );
                %                 plot(spec.lambda, conv(spec.halo(:,channel)+spec.dcx(:,channel),instfu(:,channel),'same'),  'DisplayName',['Halo+DCX - ' name] ); %+spec.brems(:,channel)
                %                 fprintf('Halo Centered at %.3f nm\n', sum(spec.lambda.*conv(spec.halo(:,channel)+spec.dcx(:,channel),instfu(:,channel),'same'))./sum(conv(spec.halo(:,channel)+spec.dcx(:,channel),instfu(:,channel),'same')));
                %                 plot(spec.lambda, conv(spec.fida(:,channel),instfu(:,channel),'same'),  'DisplayName',['FIDA - ' name] );
            else
                plot(spec.lambda,specr(:,channel), 'DisplayName', ['Spectrum - ' name] );
                disp('Supply in_data from e.g. get_bes_fida_aug_data for more plots!')
            end
            hold on
            xlabel('Wavelength [nm]')
            ylabel('Intensity [Ph/(s nm m^2 sr)]')
            set(gca,'YScale','log')
            xlim([650 663]);
            ylim([1e15, 3e19]);
            if lgeom
                %title(['Channel: ' geom.spec.id(channel)])
                disp(['Channel: ', char(geom.spec.id(channel))])
                disp(['R= ', num2str(geom.spec.radius(channel))])
            end
            legend(ax{i},'Location','best');
        case 'los3d'
            vec = [0, 0, -1];
            rotation = 0;% 67.5;
            lens = rotate_points(geom.spec.lens,vec,deg2rad(rotation));
            axi = rotate_points(geom.spec.axis,vec,deg2rad(rotation));
            los = [lens, lens + axi.*max(geom.spec.radius)*length];
            los = reshape(los,geom.spec.nchan,3,2);
            h=plot3(ax{i},squeeze(los(channel,1,:))'*fac,squeeze(los(channel,2,:))'*fac,squeeze(los(channel,3,:))'*fac,linestyle);
            src = rotate_points(geom.nbi.src,vec,deg2rad(rotation));
            axi_nbi = rotate_points(geom.nbi.axis,vec,deg2rad(rotation));
            los_nbi = [src, src + axi_nbi.*length.*sqrt(sum(src.^2)/3)];
            los_nbi = reshape(los_nbi,3,2);
            plot3(ax{i},squeeze(los_nbi(1,:))'*fac,squeeze(los_nbi(2,:))'*fac,squeeze(los_nbi(3,:))'*fac,'-r');
            plot3(ax{i},squeeze(los_nbi(1,1))'*fac,squeeze(los_nbi(2,1))'*fac,squeeze(los_nbi(3,1))'*fac,'+k');
            sname = [filename, '_', plot_type{i}];
            %writematrix(los_nbi,sname,linestyle);
            % set(h, {'DisplayName'}, cellstr(deblank(geom.spec.id(channel))))
            %legend(h,'Location','bestoutside');
            axis equal
            rotate3d on
        case 'lostor'
            vec = [0, 0, -1];
            %lens =geom.spec.lens';
            %axi = (geom.spec.lens + geom.spec.axis.*max(geom.spec.radius)*1.1)';
            lens = rotate_points(geom.spec.lens',vec,deg2rad(rotation))'; %AUG: 67.5
            axi = rotate_points((geom.spec.lens + geom.spec.axis.*max(geom.spec.radius)*length)',vec,deg2rad(rotation))';
            los = [lens, axi];
            los = reshape(los,3,geom.spec.nchan,2);
            for k = 1:size(channel,2)
                %,'Color',ax{i}.ColorOrder(k,:)
                if iscell(chan_description)
                    displ=chan_description{k};
                else
                    displ=chan_description;
                end
                h=plot(ax{i},squeeze(los(1,channel(:,k),:))'*fac,squeeze(los(2,channel(:,k),:))'*fac ,linestyle, 'DisplayName', displ);
                %set(h, 'DisplayName', chan_description{k});
            end
            %h=plot(ax,squeeze(los(1,channel,:))'*fac,squeeze(los(2,channel,:))'*fac, 'k');
            %set(h, {'DisplayName'}, cellstr(deblank(geom.spec.id(channel))));
            %legend(h,'Location','bestoutside');
            xlabel('X [cm]')
            ylabel('Y [cm]')
        case 'los2d'
            vec = [0, 0, -1];
            %lens = rotate_points(geom.spec.lens',vec,deg2rad(67.5))';
            %axi = rotate_points(geom.spec.axis',vec,deg2rad(67.5))';
            lens = rotate_points(geom.spec.lens',vec,deg2rad(rotation))';
            axi = rotate_points(geom.spec.axis',vec,deg2rad(rotation))';
            los = [lens, lens + axi.*max(geom.spec.radius)*length];
            los = reshape(los,3,geom.spec.nchan,2);
            lost = permute(los,[3,2,1]);
            los2 = reshape(lost,2,geom.spec.nchan*3);
            losre=interp1([0,1],los2,linspace(0,1,200));
            los= reshape(losre,[],size(lost,2),size(lost,3));
            for k = 1:size(channel,2)
                r = sqrt(los(:,channel(:,k),1).^2 + los(:,channel(:,k),2).^2);
                %,'Color',ax{i}.ColorOrder(k,:)
                if iscell(chan_description)
                    displ=chan_description{k};
                else
                    displ=chan_description;
                end
                h=plot(ax{i},r*fac,squeeze(los(:,channel(:,k),3))*fac,linestyle, 'DisplayName', displ);
                %set(h, 'DisplayName', chan_description{k});
            end
            %legend(h,'Location','bestoutside');

    end
    %disp(plot_type{i});
    %if numel(plot_type{i}) > 2
    if strcmp(plot_type{i}(end-1:end),'2d')
        %if ldist
        pixplot(r,z,tmp(:,:,1));
        c = colorbar;
        c.Label.String = cstring;
        xlim([r(1) r(end)])
        ylim([z(1) z(end)])
        %end
        xlabel('R [cm]')
        ylabel('Z [cm]')
    elseif strcmp(plot_type{i}(end-2:end),'tor') && ldist
        if ndims(dist.f) < 5
            disp('4D Distribution has no toroidal information')
            return;
        end
        pixplot(r,phi,squeeze(tmp(:,z0_ind,:)))
        xlabel('R [cm]')
        ylabel('Phi [rad]')
        c = colorbar;
        c.Label.String = cstring;
        xlim([r(1) r(end)])
    elseif strcmp(plot_type{i}(end-2:end),'tor')
        pixplot(r,phi,squeeze(tmp(:,z0_ind,:)))
        xlabel('R [cm]')
        ylabel('Phi [rad]')
        c = colorbar;
        c.Label.String = cstring;
        xlim([r(1) r(end)])
    end
    %end
    if lsave
        %caxis([0 3e11])
        legend(ax,'Location','best');
        sname = [filename, '_', plot_type{i}];
        %         savefig(figs{i},sname)
        %         exportgraphics(figs{i},[sname,'.png'],'Resolution',300);
        savefig(gcf,sname)
        exportgraphics(gcf,[sname,'.png'],'Resolution',300);
    end

end

end




function F = box_gauss_funct(X,A,B,C,D,E) % From /afs/ipp/home/s/sprd/XXX_DIAG/LIB
gam   = double(D);
width = double(E);
rl    = abs(0.5d0*width./gam);
Z     = abs((double(X)-double(C))./gam);
F     = double(B)*(0.5d0./width.*(erf((Z+rl)) - erf((Z-rl))))+double(A);

% Normalization and cutoff
F = F./sum(F,1);
F(F<1e-5) = 0;
end


function [xyz] = uvw_to_xyz(alpha, beta, gamma, uvw, origin) % From uvw_to_xyz.pro (D3D FIDASIM)
sa = sin(alpha); ca = cos(alpha);
sb = sin(beta) ; cb = cos(beta);
sg = sin(gamma) ; cg = cos(gamma);

R = zeros(3);
R(1,1) = ca*cb ; R(2,1) = ca*sb*sg - cg*sa ; R(3,1) = sa*sg + ca*cg*sb;
R(1,2) = cb*sa ; R(2,2) = ca*cg + sa*sb*sg ; R(3,2) = cg*sa*sb - ca*sg;
R(1,3)= -sb   ; R(2,3) = cb*sg            ; R(3,3)= cb*cg;

uvw_shifted = uvw-repmat(origin,size(uvw,1),1);
xyz = uvw_shifted*R;
end
