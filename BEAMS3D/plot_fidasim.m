function [ figs ] = plot_fidasim(filename,varargin)
%PLOT_FIDASIM Makes plots of FIDASIM data generated with the BEAMS3D
%Interface.
%The PLOT_FIDASIM function creates various canned plot of the data.  The
%
% Example usage
%      plot_beams('output','denf');

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

dist = read_hdf5(dist_name);
if ~isstruct(dist)
    disp('ERROR: Distribution file not found, check filename!');
    disp(['       Filename: ' filename]);
    data=[];
end
eq = read_hdf5(eq_name);
if ~isstruct(eq)
    disp('ERROR: Equilbirium file not found, check filename!');
    disp(['       Filename: ' filename]);
    data=[];
end

lsave = 0;
plot_type = {};
figs = {};
if nargin > 1
    i = 1;
    while i < nargin
        switch varargin{i}
            case {'overview','profiles',...
                    'denf2d','denf','fdenf','fdenf2d',...
                    'pitch','energy','ep2d',...
                    'ba','br2d','bt2d','bz2d'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case 'save'
                lsave = 1;
            case 'figs'
                i=i+1;
                figs = varargin{i};
            otherwise
                disp(['ERROR: Option ', varargin{i}, ' not found!']);
        end
        i=i+1;
    end
end

[~,z0_ind]=min(abs(dist.z));
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

            plot(ax,dist.energy, squeeze(trapz(dist.phi,trapz(dist.z,trapz(dist.r,trapz(dist.pitch,dist.f,2),3),4),5)),'DisplayName','Pitch');
            xlabel('Energy [keV]')
            ylabel('Fast Ion Distribution [1/keV]')
        case 'pitch'
            dr = dist.r(2)-dist.r(1);
            dz = dist.z(2)-dist.z(1);
            dphi = dist.phi(2) - dist.phi(1);
            % Area
            %area=dr*100.*dz*100;
            area=dr.*dz;
            % Volume (function of R)
            vol = dist.r.*dphi.*area;
            vol2d=repmat(vol,[1 dist.nz dist.nphi]);

            n_fida = sum(dist.denf.*vol2d,'all');

            fprintf('Total fast ions in %s: %3.2e\n',filename,n_fida);
            tmp = squeeze(trapz(dist.phi,trapz(dist.z,trapz(dist.energy,dist.f,1),4),5));
            plot(ax,dist.pitch, squeeze(trapz(dist.r,repmat(dist.r',size(dist.f,2),1).*tmp,2)),'DisplayName','Pitch');

            %plot(dist.pitch, squeeze(trapz(dist.phi,trapz(dist.z,trapz(dist.r,trapz(dist.energy,dist.f,1),3),4),5)),'DisplayName','Pitch');
            xlabel('Pitch [-]')
            ylabel('Fast Ion Distribution [-]')
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
        case 'denf'
            plot(ax,dist.r, squeeze(dist.denf(:,z0_ind,1)), 'DisplayName','Denf ');
            xlabel('R [m]')
            ylabel('Fast ion density [m^{-3}]')
            title('Fast ion density profile at z=0')
        case 'fdenf'
            tmp = squeeze(trapz(dist.pitch,trapz(dist.energy,dist.f,1),2));
            plot(ax,dist.r, tmp(:,z0_ind,1), 'DisplayName',['Denf ',dist_name]);
            xlabel('R [m]')
            ylabel('Fast ion density [m^{-3}]')
            title('Fast ion density profile at z=0')
        case 'ep2d'
            pixplot(dist.energy,dist.pitch,squeeze(trapz(dist.phi,trapz(dist.z,trapz(dist.r,dist.f,3),4),5)))
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
    end
    if strcmp(plot_type{i}(end-1:end),'2d')
        pixplot(dist.r,dist.z,tmp(:,:,1))
        xlabel('R [m]')
        ylabel('Z [m]')
        c = colorbar;
        c.Label.String = cstring;
    end
    if lsave
        sname = [filename, '_', plot_type{i}];
        savefig(figs{i},sname)
        exportgraphics(figs{i},[sname,'.png'],'Resolution',300);
    end

end

end



