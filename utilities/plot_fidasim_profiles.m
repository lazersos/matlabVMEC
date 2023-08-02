function plot_data = plot_fidasim_profiles(filename,sim_data,varargin)
%PLOT_FIDASIM_PROFILES plots radial profiles from the FIDASIM code. As
%these are used to compare to experimental data, the integration range and
%dispersion is set by variables in the sim_data structure generated e.g. by
%the get_bes_fida_aug_data function. Plots are shared between the
%two functions.
%
% Example usage
%      [~,sim_data] = get_bes_fida_aug_data(filename,'t_point',3.5,'fidabes');
%      plot_fidasim_profiles(filename,sim_data,'X);
%      !!! 'X' can be 'fida', 'bes', or 'fidabes'
%
% Miscellaneous Arguments
%      plot_fidasim(runid,'mean'); %Apply moving mean to spectrum
%      plot_fidasim(runid,'sim_data',sim_data); %Used for dispersion
%      plot_fidasim(runid,'save'); %Export figures (.fig and .png)
%      plot_fidasim(runid,'name', 'test'); %ID Name for legend
%      plot_fidasim(runid,'fac', 1.0); %Scaling factor
%
bg_range=sim_data.bg_range;

bes_range = sim_data.bes_range;
fida_range = sim_data.fida_range;
dispersion = sim_data.dispersion;
lambda_dat = sim_data.lambda;
if isfield(sim_data,'dex')
    dex = sim_data.dex;
else
    dex = [];
end

if isfield(sim_data,'ax')
    ax=sim_data.ax;
else
    ax = {};
end

lsave = 0;
lmean = 0;
plot_type = {};
linestyle = '+';
fac = 1;
name = filename(1:end-2);
if nargin > 2
    i = 1;
    while i < nargin-1
        switch varargin{i}
            case {'FIDA','BES','FIDABES','fida','bes','fidabes','bck'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case 'mean'
                lmean =1;
            case 'save'
                lsave = 1;
            case 'avg_frames'
                i = i+1;
                avg_frames = varargin{i};
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


data=read_hdf5(filename);
% R = h5read(filename,'/radius');
% full = h5read(filename,'/full');
% half= h5read(filename,'/half');
% third= h5read(filename,'/third');
% halo= h5read(filename,'/halo');
% dcx= h5read(filename,'/dcx');
% fida= h5read(filename,'/fida');
% pfida= h5read(filename,'/pfida');
% brems= h5read(filename,'/brems');
% lambda = h5read(filename,'/lambda');
R = data.radius;
full =data.full;
half= data.half;
third= data.third;
halo= data.halo;
dcx= data.dcx;
fida= data.fida;
brems= data.brems;
lambda =data.lambda;
if isfield(data,'pfida')
    pfida=data.pfida;
else
    pfida=zeros(size(fida));
end

spec = full + half + third + halo + dcx + fida+pfida;% + brems;

cwav_mid=mean(lambda);
instfu = box_gauss_funct(lambda,0.,1.,cwav_mid,sim_data.instfu_gamma,sim_data.instfu_box_nm);
disp(['Applying Instrument function to FIDASIM data: ', filename]);
for i = 1:size(spec,2)
    spec(:,i) = conv(spec(:,i),instfu(:,i),'same');    
end

if lmean == 1
    spec = movmean(spec,15);
    disp(['Applying moving mean to FIDASIM data: ', filename]);
end

% dispersion_tmp = diff(lambda_dat,1);
dispersion_tmp = diff(lambda,1);
dispersion_tmp = [dispersion_tmp; dispersion_tmp(end,:)];
dispersion_tmp = repmat(dispersion_tmp,1,size(spec,2));

% bes_dex = (lambda_dat > repmat(bes_range(:,1)',size(lambda_dat,1),1)) & (lambda_dat < repmat(bes_range(:,2)',size(lambda_dat,1),1));
% fida_dex = (lambda_dat > fida_range(1)) & (lambda_dat < fida_range(2));
bg_dex = (lambda > bg_range(1)) & (lambda < bg_range(2));
bes_dex = (lambda > repmat(bes_range(:,1)',size(lambda,1),1)) & (lambda < repmat(bes_range(:,2)',size(lambda,1),1));
fida_dex = (lambda > fida_range(1)) & (lambda < fida_range(2));

% for i = 1:size(lambda_dat,2)
% dispersion_tmp(:,i) = interp1(lambda_dat(:,i),dispersion(:,i),lambda);
% end
%
% bes = sum(spec.*dispersion_tmp.*bes_dex,1);
% fida = sum(spec.*dispersion_tmp.*fida_dex,1);

bg = sum(brems.*dispersion_tmp.*bg_dex,1,'omitnan');
bes = sum(spec.*dispersion_tmp.*bes_dex,1,'omitnan');
fida = sum(spec.*dispersion_tmp.*fida_dex,1,'omitnan');

for i = 1:size(plot_type,2)
    if i>numel(ax)
        figure;
        ax{i} = gca;
        hold on;
    end
    legend(ax{i},'Location','best');
    switch lower(plot_type{i})
        case 'bck'
            tmp = bg(dex);
            ystr = 'BG';
        case 'bes'
            tmp = bes(dex);
            ystr = 'BES';
        case 'fida'
            tmp = fida(dex);
            ystr = 'FIDA';
        case 'fidabes'
            tmp = fida(dex)./bes(dex);
            ystr = 'FIDA/BES';

    end

    if fac~=1
        tmp = tmp.*fac;
        dispname = ['', name, ', scaling factor: ' num2str(fac)];
    else
        dispname = ['', name];
    end
    plot(ax{i},R(dex), tmp,linestyle,'DisplayName',dispname, 'LineWidth',2.0);
    xlabel(ax{i},'R [cm]')
    ylabel(ax{i},ystr)
    if lsave
        legend(ax{i},'Location','best');
        sname = [name, '_', plot_type{i}];
        savefig(gcf,sname)
        exportgraphics(gcf,[sname,'.eps'],'Resolution',300);
    end
end
plot_data.R = R;
plot_data.bes = bes;
plot_data.fida = fida;
plot_data.spec = spec;
plot_data.lambda=lambda;
plot_data.ax=ax;
plot_data.dex = dex;

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