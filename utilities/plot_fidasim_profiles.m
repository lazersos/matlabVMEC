function plot_data = plot_fidasim_profiles(filename,in_data,varargin)
%PLOT_FIDASIM_PROFILES plots radial profiles from the FIDASIM code. As
%these are used to compare to experimental data, the integration range and
%dispersion is set by variables in the in_data structure generated e.g. by
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
bes_range = in_data.bes_range;
fida_range = in_data.fida_range;
dispersion = in_data.dispersion;
lambda_dat = in_data.lambda;
if isfield(in_data,'dex')
    dex = in_data.dex;
else
    dex = [];
end

if isfield(in_data,'figs')
    figs=in_data.figs;
else
    figs = {};
end

lsave = 0;
lmean = 0;
plot_type = {};

fac = 1;
name = filename(1:end-2);
if nargin > 2
    i = 1;
    while i < nargin-1
        switch varargin{i}
            case {'FIDA','BES','FIDABES','fida','bes','fidabes'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case 'mean'
                lmean =1;
            case 'save'
                lsave = 1;
            case 'avg_frames'
                i = i+1;
                avg_frames = varargin{i};
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
% [R,I] = sort(R);
% full = full(:,I);
% half = half(:,I);
% third = third(:,I);
% dcx = dcx(:,I);
% fida = fida(:,I);
% brems = brems(:,I);

spec = full + half + third + halo + dcx + fida;% + brems;

if lmean == 1
    spec = movmean(spec,15);
    disp(['Applying moving mean to FIDASIM data: ', filename]);
end

%spect = interp1(lambda, (spec').', (lambda_dat').'); %Vectorized interpolation
for i = 1:size(lambda_dat,2)
    spectmp(:,i) = interp1(lambda, spec(:,i), lambda_dat(:,i),'spline');
end

for i = 1:size(spec,2)
    spectmp(:,i) = conv(spectmp(:,i),in_data.instfu(:,i),'same');
end

dispersion_tmp = diff(lambda_dat,1);
dispersion_tmp = [dispersion_tmp; dispersion_tmp(end,:)];

bes_dex = (lambda_dat > repmat(bes_range(:,1)',size(lambda_dat,1),1)) & (lambda_dat < repmat(bes_range(:,2)',size(lambda_dat,1),1));
fida_dex = (lambda_dat > fida_range(1)) & (lambda_dat < fida_range(2));


% for i = 1:size(lambda_dat,2)
% dispersion_tmp(:,i) = interp1(lambda_dat(:,i),dispersion(:,i),lambda);
% end
%
% bes = sum(spec.*dispersion_tmp.*bes_dex,1);
% fida = sum(spec.*dispersion_tmp.*fida_dex,1);


bes = sum(spectmp.*dispersion_tmp.*bes_dex,1,'omitnan');
fida = sum(spectmp.*dispersion_tmp.*fida_dex,1,'omitnan');

for i = 1:size(plot_type,2)
    if i>numel(figs)
        figs{i}=figure;
        ax = gca;
        hold on;
    else
        allAxesInFigure = findall(figs{i},'type','axes');
        ax = allAxesInFigure(~ismember(get(allAxesInFigure,'Tag'),{'legend','Colobar'}));
    end
    legend(ax,'Location','best');
    switch lower(plot_type{i})
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
        dispname = ['FIDASIM ', name, ', scaling factor: ' num2str(fac)];
    else
        dispname = ['FIDASIM ', name];
    end
    plot(ax,R(dex), tmp,'+','DisplayName',dispname, 'LineWidth',2.0);
    xlabel(ax,'R [cm]')
    ylabel(ax,ystr)
    if lsave
        legend(ax,'Location','best');
        sname = [name, '_', plot_type{i}];
        savefig(figs{i},sname)
        exportgraphics(figs{i},[sname,'.png'],'Resolution',300);
    end
end
plot_data.R = R;
plot_data.bes = bes;
plot_data.fida = fida;
plot_data.spec = spectmp;
plot_data.lambda=lambda;
plot_data.figs=figs;
plot_data.dex = dex;

end
