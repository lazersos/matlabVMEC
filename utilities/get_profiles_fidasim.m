function plot_data = plot_fidasim_profiles(filename,bes_range, fida_range,dispersion,lambda_dat)




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

lsave = 0;
plot_type = {};
figs = {};
fac = 1;
name = filename(1:end-2);
if nargin > 2
    i = 1;
    while i < nargin-1
        switch varargin{i}
            case {'FIDA','BES','FIDABES','fida','bes','fidabes'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
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

        end
        i=i+1;
    end
end



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
    end
    plot(ax,R(dex), tmp ,'DisplayName',['FIDASIM ', name], 'LineWidth',2.0);
    xlabel(ax,'R [cm]')
    ylabel(ax,ystr)
    if lsave
        legend(ax);
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

end
