function [plot_data,sim_data] = get_bes_fida_aug_data(filename,t_point,varargin)
%This function expects an hdf5 file with CXRS data from AUG shotfiles
%exported as follows (for shot 38581):
% idl
% load_cfr,38581,data,err=err,raw=0,/newstructure,exp=exp
% write_hdf5,data,filename='38581_CFR.h5'
% It is used to compare FIDASIM to experimental data


shotid = str2double(filename(1:5));
plot_type = {};
figs = {};
lsave = 0;
avg_time=0;
t_passive =0;

R = h5read(filename,'/r_pos');
spec_in= h5read(filename,'/intens');
spec_err_in= h5read(filename,'/intenserr');
lambdatmp = h5read(filename,'/cor_wavel');
lambda = lambdatmp;
names = h5read(filename,'/los_name');
dispersion_in=h5read(filename,'/dispersion');
time = h5read(filename,'/time_arr');

if nargin > 2
    i = 1;
    while i < nargin-1
        switch varargin{i}
            case {'FIDA','BES','FIDABES','fida','bes','fidabes'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case 'channels'
                i=i+1;
                tmp = cell2mat(names);
                dex_in = strcmp(tmp(:,1:3),string(varargin{i}));
                if strcmp("F50",string(varargin{i}))
                    fida_range = [652.5,653.5];
                else
                    fida_range = [659.5, 660.5];
                end
            case{'spectrum'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                i=i+1;
                channel = varargin{i};
            case 't_passive'
                i = i+1;
                t_passive = varargin{i};
            case 'save'
                lsave = 1;
            case 'figs'
                i=i+1;
                figs = varargin{i};
            case 'name'
                i = i+1;
                name = varargin{i};
            case 'avg_time'
                i = i+1;
                avg_time = varargin{i};

        end
        i=i+1;
    end
end

bg_range = [664.5, 666];

bes_range = [  654.621      655.086
    654.797      655.267
    655.009      655.487
    656.204      656.691
    656.438      656.934
    656.669      657.174
    656.907      657.419
    657.135      657.656
    654.844      655.320
    656.124      656.611
    656.453      656.949
    656.787      657.291
    657.115      657.628
    657.442      657.964
    658.290      658.815
    658.200      658.717
    658.107      658.616
    658.013      658.514
    657.908      658.400
    657.812      658.296
    657.716      658.193
    657.622      658.091
    657.526      657.987
    657.430      657.884
    658.256      658.731
    658.839      659.327];

if shotid < 32000
    bes_range = [658.906      658.906
        658.802      658.802
        658.703      658.703
        658.487      658.487
        658.270      658.270
        658.063      658.063
        655.332      655.332
        656.676      656.676
        657.014      657.014
        657.352      657.352
        657.689      657.689
        658.021      658.021
        658.354      658.354
        658.872      658.872
        658.774      658.774
        658.670      658.670
        658.566      658.566
        658.446      658.446
        658.341      658.341
        658.237      658.237
        658.132      658.132
        658.025      658.025]; %#30841
    %fidasim_names = {'CER-01', 'CER-04', 'CER-07', 'CER-13', 'CER-19', 'CER-25', 'F50-1', 'F50-2', 'F50-3', 'F50-4', 'F50-5', 'F50-6', 'F50-7', 'CER-02', 'CER-05', 'CER-08', 'CER-11', 'CER-14', 'CER-17', 'CER-20', 'CER-23', 'CER-26'};
    bes_range=bes_range(1:22,:);
    bes_range(:,1) = bes_range(:,1)-0.5;
end


time_dex = (t_point + avg_time/2 >= time) & (t_point -avg_time/2 <= time);
time_dex = permute(repmat(time_dex',size(spec_in,2),1,size(spec_in,1)),[3,1,2]);
disp(['Max. Frames used for averaging of ', filename,': ', num2str(max(sum(time_dex(1,:,:),3)))]);
if t_passive~=0
    time_dex_passive = (t_passive + avg_time/2 >= time) & (t_passive -avg_time/2 <= time);
    time_dex_passive = permute(repmat(time_dex_passive',size(spec_in,2),1,size(spec_in,1)),[3,1,2]);
end
dex = permute(repmat(dex_in,1,numel(time),size(spec_in,1)),[3,1,2]);

bes_dex = (lambda > repmat(bes_range(:,1)',size(lambda,1),1)) & (lambda < repmat(bes_range(:,2)',size(lambda,1),1));
fida_dex = (lambda > fida_range(1)) & (lambda < fida_range(2));
bg_dex = (lambda > bg_range(1)) & (lambda < bg_range(2));

bes_dex = repmat(bes_dex,1,1,numel(time));
fida_dex = repmat(fida_dex,1,1,numel(time));
bg_dex = repmat(bg_dex,1,1,numel(time));
dispersion = repmat(dispersion_in,1,1,numel(time));

spec = spec_in - repmat(sum(spec_in.*bg_dex,1)./sum(bg_dex,1),size(spec_in,1),1); %Background subtraction

if t_passive~=0
    spec_passive = repmat(sum(spec.*dex.*time_dex_passive,3)./sum(dex.*time_dex_passive,3),1,1,size(spec_in,3));
    spec = spec - spec_passive;
end



%plot(squeeze(spec(:,16,1000)))

bes = sum(spec.*dispersion.*bes_dex,1);
fida = sum(spec.*dispersion.*fida_dex,1);

bg_err = sqrt(sum(spec_err_in.^2.*bg_dex,1))./sum(bg_dex,1);
spec_err = sqrt(spec_err_in.^2 + bg_err.^2);

bes_err = sqrt(sum((spec_err.*dispersion.*bes_dex).^2,1));
fida_err =sqrt(sum((spec_err.*dispersion.*fida_dex).^2,1));
fida_bes_err = squeeze(sqrt((1.0./bes).^2.*fida_err.^2+(fida./bes.^2).^2.*bes_err.^2));
%fida_bes_err = sqrt(sum((fida_bes_err.*or(bes_dex,fida_dex)).^2,1));
%bg_err = squeeze(bg_err);
bes = squeeze(bes);
fida = squeeze(fida);
bes_err = squeeze(bes_err);
fida_err = squeeze(fida_err);
fida_bes_err = squeeze(fida_bes_err);


chandex =  squeeze(any(dex,[1,3]))'; %only channel dex
% if t_passive~=0
% time_dex_passive = (t_passive + avg_time/2 >= time) & (t_passive -avg_time/2 <= time);
% time_dex_passive = repmat(time_dex_passive',size(dex,1),1);
% bes = sum((bes.*dex.*time_dex)-(bes.*dex.*time_dex_passive),2)./sum(dex.*time_dex,2);
% fida =sum((fida.*dex.*time_dex)-(fida.*dex.*time_dex_passive),2)./sum(dex.*time_dex,2);
% spec = (sum(spec.*dex.*time_dex,3)-sum(spec.*dex.*time_dex_passive,3))./sum(dex.*time_dex,2);
% else

dex = squeeze(dex(1,:,:));
time_dex = squeeze(time_dex(1,:,:));
% if t_passive~=0
% time_dex_passive = squeeze(time_dex_passive(1,:,:));
% time_dex_passive = and(dex,time_dex_passive);
% end
bes_out = sum(bes.*dex.*time_dex,2,'omitnan')./sum(time_dex,2);
fida_out =sum(fida.*dex.*time_dex,2,'omitnan')./sum(time_dex,2);
%end
bes_err_out = sqrt(sum(bes_err.^2.*dex.*time_dex,2,'omitnan').^2 + std(bes_err.*dex.*time_dex,0,2).^2)./sum(time_dex,2);
fida_err_out = sqrt(sum(fida_err.^2.*dex.*time_dex,2,'omitnan').^2 + std(fida_err.*dex.*time_dex,0,2).^2)./sum(time_dex,2);
fida_bes_err_out = sqrt(sum(fida_bes_err.^2.*dex.*time_dex,2,'omitnan').^2 + std(fida_bes_err.*dex.*time_dex,0,2).^2)./sum(time_dex,2);


R_pts = repmat(R,1,size(spec_in,3))*100;
time_dex = and(dex,time_dex);

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
            plot(ax,R_pts(time_dex),bes(time_dex),'o','DisplayName',['Data ', num2str(t_point - avg_time/2),'-',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
            %plot(ax,R_pts(time_dex_passive),bes(time_dex_passive),'o','DisplayName',['Data Passive ', num2str(t_passive- avg_time/2),'-',num2str(t_passive+ avg_time/2), 's'], 'LineWidth',2.0);
            tmp =bes_out;
            tmp_err = bes.*dex.*time_dex;%bes_err_out;
                        tmp_err(tmp_err==0) = NaN;
            tmp_err = std(tmp_err,0,2,'omitnan');
            ystr = 'BES';

        case 'fida'
            plot(ax,R_pts(time_dex),fida(time_dex),'o','DisplayName',['Data ', num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
            %plot(ax,R_pts(time_dex_passive),fida(time_dex_passive),'o','DisplayName',['Data Passive ', num2str(t_passive- avg_time/2),'-',num2str(t_passive+ avg_time/2), 's'], 'LineWidth',2.0);
            tmp = fida_out;
            tmp_err = fida.*dex.*time_dex;%fida_err_out;
                        tmp_err(tmp_err==0) = NaN;
            tmp_err = std(tmp_err,0,2,'omitnan');
            ystr = 'FIDA';
        case 'fidabes'
            plot(ax,R_pts(time_dex),fida(time_dex)./bes(time_dex),'o','DisplayName',['Data ', num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
            % plot(ax,R_pts(time_dex_passive),fida(time_dex_passive)./bes(time_dex_passive),'.','DisplayName',['Points Passive ', num2str(t_passive), 's'], 'LineWidth',2.0);
            tmp = fida_out./bes_out;
            tmp_err = (fida./bes).*dex.*time_dex;
            tmp_err(tmp_err==0) = NaN;
            tmp_err = std(tmp_err,0,2,'omitnan');%fida_bes_err_out;
            ystr = 'FIDA/BES';
            ylim([0 0.16])
            
        case 'spectrum'
            tmp =squeeze(sum(spec.*permute(repmat(time_dex,1,1,size(spec,1)),[3,1,2]),3)./sum(permute(repmat(time_dex,1,1,size(spec,1)),[3,1,2]),3));
             plot(ax,lambda(:,channel),tmp(:,channel), 'DisplayName',['Data ',  num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's, Chan: ', names{channel} ], 'LineWidth',2.0);
             %plot(ax,lambda(:,channel),tmp(:,channel), 'DisplayName',['Data ',  num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's, Chan: ', names{channel} ], 'LineWidth',2.0);
             %tmp_in =squeeze(sum(spec_in.*permute(repmat(time_dex,1,1,size(spec_in,1)),[3,1,2]),3)./sum(permute(repmat(time_dex,1,1,size(spec_in,1)),[3,1,2]),3));
             %plot(ax,lambda(:,channel),tmp_in(:,channel), 'DisplayName',['Data, bo BG sub ', num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
             xline(bes_range(channel,:),'DisplayName', 'BES Range')
    end

    if ~strcmp(plot_type{i},'spectrum')
errorbar(ax,R(chandex)*100, tmp(chandex),tmp_err(chandex),'--','DisplayName',['Avg. ', num2str(t_point), 's'], 'LineWidth',2.0);
    xlabel(ax,'R [cm]')
    ylabel(ax,ystr);
    end
    if lsave
        legend(ax);
        sname = [num2str(shotid), '_', plot_type{i}];
        savefig(figs{i},sname)
        exportgraphics(figs{i},[sname,'.png'],'Resolution',300);
    end
end


plot_data.R = R;
plot_data.bes = bes_out;
plot_data.fida = fida_out;
plot_data.spec = spec;
plot_data.bes_err = bes_err_out;
plot_data.fida_err = fida_err_out;
plot_data.fida_bes_err = fida_bes_err_out;
plot_data.names = names;
plot_data.lambda=lambdatmp;
sim_data.figs = figs;
sim_data.lambda=lambdatmp;
sim_data.dispersion = dispersion_in;
sim_data.bes_range = bes_range;
sim_data.fida_range = fida_range;
sim_data.dex = chandex;
end