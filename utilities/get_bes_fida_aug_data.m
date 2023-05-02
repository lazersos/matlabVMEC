function [plot_data,sim_data] = get_bes_fida_aug_data(filename,varargin)
%GET_BES_FIDA_AUG_DATA plots and analyses the calibrated FIDA data exported
%from ASDEX upgrade. Various canned plots for exploring the data are
%available, i.e. FIDA/BES radial profiles and time-traces. The standard
% FIDA integration range is [660, 661]. The spectrum can
%also be plotted individually.
% Averaging time length and passive subtraction time can be supplied
% optionally and error bars are calculated.
%This function expects an hdf5 file with CXRS data from AUG shotfiles
%exported as follows (for shot 38581):
% idl
% load_cfr,38581,data,err=err,raw=0,/newstructure,exp=exp
% write_hdf5,data,filename='38581_CFR.h5'
%
% Example usage:
%      get_bes_fida_aug_data(filename,'t_point',3.5); %Time for profiles
%      get_bes_fida_aug_data(filename,'t_passive',3.5); %Time for passive subtraction
%      get_bes_fida_aug_data(filename,'avg_time',0.01); %Averaging window length
%      get_bes_fida_aug_data(filename,'channels','CER'); %Channel selector (3 characters)
%      get_bes_fida_aug_data(filename,'X'); %Profile
%      get_bes_fida_aug_data(filename,'timetrace_X'); %Timetrace
%      !!! 'X' can be 'fida', 'bes', or 'fidabes' and 'spectrum'
% Miscellaneous Arguments
%      plot_fidasim(runid,'ax',ax); %Used for overplotting
%      plot_fidasim(runid,'save'); %Export figures (.fig and .png)
%      plot_fidasim(runid,'name', 'test'); %ID Name for legend
%
% Integrated examples:
% Plot FIDA, BES, and FIDA/BES profiles at 3.5s with passive signal at 3.6s
% and 0.01s averaging time. Use the 'CER' channels for the profile
%      get_bes_fida_aug_data(filename,'t_point',3.5,'avg_time',0.010,'t_passive',3.6,'channels','CER','fidabes','fida','bes');


shotid = str2double(filename(1:5));
plot_type = {};
ax = {};
lsave = 0;
avg_time=0;
lmovie=0;
t_passive =0;
t_point = 0;
fida_range = [660, 661];

dex_in = '';

R = h5read(filename,'/r_pos');

spec_err_in= h5read(filename,'/intenserr');
lambda = h5read(filename,'/cor_wavel');%-0.055;
spec_in= h5read(filename,'/intens');
names_unsorted = h5read(filename,'/los_name');
dispersion_in=h5read(filename,'/dispersion');
time = h5read(filename,'/time_arr');
instfu_gamma = h5read(filename,'/instfu_gamma')';
instfu_box_nm = h5read(filename,'/instfu_box_nm')';
instfu_gauss_nm = h5read(filename,'/instfu_gauss_nm')';
cwav_mid = interp1(1:size(lambda,1),lambda,[size(lambda,1)/2.]);
instfu = box_gauss_funct(lambda,0.,1.,cwav_mid,instfu_gamma,instfu_box_nm);


% [R,I] = sort(R);
% spec_in = spec_in(:,I,:);
% spec_err_in = spec_err_in(:,I,:);
% lambda = lambda(:,I);
% lambdatmp = lambdatmp(:,I);
% dispersion_in = dispersion_in(:,I);
names = names_unsorted;%(I);
tmp = cell2mat(names);
dex_in = true(size(names));%strcmp(tmp(:,1:3),"CER");
gray = [.5 .5 .5];
if nargin > 2
    i = 1;
    while i < nargin
        switch varargin{i}
            case {'FIDA','BES','FIDABES','fida','bes','fidabes'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
            case 'channels'
                i=i+1;
                dex_in = strcmp(tmp(:,1:3),string(varargin{i}));
                if strcmp("F50",string(varargin{i}))
                    fida_range = [652.5,653.5];
                else
                    fida_range = [660, 661];
                end
            case{'spectrum', 'timetrace_fida', 'timetrace_bes', 'timetrace_fidabes'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible
                i=i+1;
                if ischar(varargin{i}) %TODO: Allow for multiple channels to be specified
                    channel = find(strcmp(varargin{i},cellstr(deblank(names'))));
                elseif iscell(varargin{i})
                    for j = 1:numel(varargin{i})
                        channel(j) = find(strcmp(varargin{i}{j},cellstr(deblank(names'))));
                    end
                else
                    %channel =find(~(I-varargin{i}));% I(varargin{i}); %Inverse lookup,
                    channel = varargin{i};
                end
            case 't_point'
                i=i+1;
                t_point = varargin{i};
            case 't_passive'
                i = i+1;
                t_passive = varargin{i};
            case 'save'
                lsave = 1;
            case 'ax'
                i=i+1;
                ax = varargin{i};
            case 'name'
                i = i+1;
                name = varargin{i};
            case 'avg_time'
                i = i+1;
                avg_time = varargin{i};
            case 'movie'
                lmovie=1;
                i = i+1;
                t = varargin{i};
            otherwise
                disp(['ERROR: Option ', varargin{i}, ' not found!']);

        end
        i=i+1;
    end
end

bg_range = [664.5, 666];

bes_range = [ 654.56110	655.05627
654.73962	655.24072
655.97443	656.48553
656.21759	656.73773
656.45605	656.98492
656.69189	657.22961
656.93341	657.48016
657.16632	657.72180
654.78632	655.29498
656.13635	656.65588
656.47156	657.00012
656.81158	657.34930
657.14575	657.69293
657.47827	658.03491
658.34113	658.90131
658.24969	658.80139
658.15497	658.69794
658.05914	658.59344
657.95184	658.47662
657.85443	658.37079
657.75659	658.26471
657.66028	658.16034
657.56232	658.05450
657.46442	657.94891
658.30511	658.81262
658.89917	659.41998];

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

%bes_range=bes_range(I,:);
if avg_time==0
    avg_time=time(2)-time(1);
end


if t_point ~=0
    time_dex = (t_point + avg_time/2 >= time) & (t_point -avg_time/2 <= time);
    time_dex_vec = logical(sum(time_dex,2));
    time_dex = permute(repmat(time_dex_vec',size(spec_in,2),1,size(spec_in,1)),[3,1,2]);
    disp(['Frames used for averaging of ', filename,': ', num2str(max(sum(time_dex(1,:,:),3))), ', frames: ', sprintf('%d,',(find(time_dex_vec)))]);
end
if t_passive~=0
    if numel(t_passive)==1
        dt = time(2)-time(1);
        time_dex_passive = and(((t_passive +dt/2) >= time),((t_passive-dt/2) < time));%and(((t_passive + avg_time/2) >= time),((t_passive -avg_time/2) < time));
    elseif numel(t_passive)==2
        time_dex_passive = and((t_passive(2) >= time),(t_passive(1) < time));%and(((t_passive + avg_time/2) >= time),((t_passive -avg_time/2) < time));
    end
    time_dex_passive_vec = time_dex_passive;
    time_dex_passive = permute(repmat(time_dex_passive,1,1,size(spec_in,2),size(spec_in,1)),[4,3,1,2]);
    disp(['Frames used for averaging passive sig. of ', filename,': ', num2str(max(sum(time_dex_passive(1,:,:),3))), ', frames: ', sprintf('%d,',(find(time_dex_passive_vec)))]);
    if size(t_passive,2) > 1
        closest_to_passive = interp1(t_passive,t_passive,time,'nearest','extrap');
        closest_to_passive = (closest_to_passive-t_passive);
        closest_to_passive(abs(closest_to_passive)<1e-3) =0;
        closest_to_passive = sum(~closest_to_passive.*(1:numel(t_passive)),2);
        closest_to_passive(closest_to_passive==0) = numel(t_passive)+1;
    else
        closest_to_passive = ones(size(time));
    end
end
if ~(strcmp(dex_in,''))
    dex = permute(repmat(dex_in,1,numel(time),size(spec_in,1)),[3,1,2]);
    %chandex =  squeeze(any(dex,[1,3]))'; %only channel dex
else
    dex = ones(size(spec_in));
end

bes_dex = (lambda > repmat(bes_range(:,1)',size(lambda,1),1)) & (lambda < repmat(bes_range(:,2)',size(lambda,1),1));
fida_dex = (lambda > fida_range(1)) & (lambda < fida_range(2));
bg_dex = (lambda > bg_range(1)) & (lambda < bg_range(2));

bes_dex = repmat(bes_dex,1,1,numel(time));
fida_dex = repmat(fida_dex,1,1,numel(time));
bg_dex = repmat(bg_dex,1,1,numel(time));
dispersion = repmat(dispersion_in,1,1,numel(time));

spec = spec_in - repmat(sum(spec_in.*bg_dex,1)./sum(bg_dex,1),size(spec_in,1),1); %Background subtraction
 %spec=spec_in;
spec(spec<0) = 0;
if t_passive~=0
    for k = 1:size(time_dex_passive,4)
        spec_passive(:,:,k) = sum(spec.*squeeze(time_dex_passive(:,:,:,k)),3)./sum(squeeze(time_dex_passive(:,:,:,k)),3);
    end
    spec_passive(:,:,k+1) = zeros(size(lambda));
    spec_passive(spec_passive<0) = 0;
    spec = spec - spec_passive(:,:,closest_to_passive);
end

spec(spec<0) = 0;

%plot(squeeze(spec(:,16,1000)))

bes = sum(spec.*dispersion.*bes_dex,1);
fida = sum(spec.*dispersion.*fida_dex,1);

bg_err = sqrt(sum(spec_err_in.^2.*bg_dex,1))./sum(bg_dex,1);
spec_err = sqrt(spec_err_in.^2 + bg_err.^2);

bes_err = sqrt(sum((spec_err.*dispersion.*bes_dex),1).^2);
fida_err =sqrt(sum((spec_err.*dispersion.*fida_dex),1).^2);
fida_bes_err = squeeze(sqrt((1.0./bes).^2.*fida_err.^2+(fida./bes.^2).^2.*bes_err.^2));
%fida_bes_err = sqrt(sum((fida_bes_err.*or(bes_dex,fida_dex)).^2,1));
%bg_err = squeeze(bg_err);
bes = squeeze(bes);
fida = squeeze(fida);
bes_err = squeeze(bes_err);
fida_err = squeeze(fida_err);
fida_bes_err = squeeze(fida_bes_err);


if ~(strcmp(dex_in,''))
    dex = squeeze(dex(1,:,:));
    %disp(names(any(dex,2)));
end
if t_point ~=0
    time_dex = squeeze(time_dex(1,:,:));
    time_dex = and(dex,time_dex);
    time_dex_inds = find(time_dex_vec);
end

R_pts = repmat(R,1,size(spec_in,3))*100;


for i = 1:size(plot_type,2)
    if i>numel(ax)
        figure;
        ax{i} = gca;
        hold on;
    else
        %allAxesInFigure = findall(figs{i},'type','axes');
        %ax = allAxesInFigure(~ismember(get(allAxesInFigure,'Tag'),{'legend','Colobar'}));
    end
    %legend(ax{i},'Location','best');
    switch lower(plot_type{i})
        case 'bes'
            %plot(ax{i},R_pts(time_dex),bes(time_dex),'o','DisplayName',['Data ', num2str(t_point - avg_time/2),'-',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
            %plot(ax{i},R_pts(time_dex_passive),bes(time_dex_passive),'o','DisplayName',['Data Passive ', num2str(t_passive- avg_time/2),'-',num2str(t_passive+ avg_time/2), 's'], 'LineWidth',2.0);
            tmp_frames=bes(time_dex);
            tmp = sum(bes.*time_dex,2)./sum(time_dex,2);
            tmp_err=sqrt(sum(bes_err.*time_dex,2).^2)./sum(time_dex,2);
            ystr = 'BES';
        case 'fida'
            %plot(ax{i},R_pts(time_dex),fida(time_dex),'o','DisplayName',['Data ', num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
            tmp_frames=fida(time_dex);
            tmp = sum(fida.*time_dex,2)./sum(time_dex,2);
            tmp_err=sqrt(sum(fida_err.*time_dex,2).^2)./sum(time_dex,2);
            ystr = 'FIDA';
%             t = annotation('textbox',[0.15 0.61 0.3 0.3],'String',['FIDA Int. Range: [', num2str(fida_range),'] nm'],'FitBoxToText','on','LineStyle','none');
%             t.LineWidth = 0.01;
        case 'fidabes'
            %compose('Data %.3f - %.3f s', (t_point - avg_time/2)',(t_point + avg_time/2)')
%             for j = 1:sum(time_dex_vec)
%                 plot(ax{i},R_pts(dex_in,time_dex_inds(j)),fida(dex_in,time_dex_inds(j))./bes(dex_in,time_dex_inds(j)),'o','Color', gray,'DisplayName',sprintf('Data %.3f s', time(time_dex_inds(j))),'LineWidth',2.0);
%             end
            tmp_frames=fida(time_dex)./bes(time_dex);
            tmp = sum(fida./bes.*time_dex,2,'omitnan')./sum(time_dex,2,'omitnan');
            %plot(ax,R_pts(time_dex),fida(time_dex)./bes(time_dex),'o','DisplayName',['Data ', num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
            %tmp = sum(fida./bes.*time_dex,2,'omitnan')./sum(time_dex,2);
            tmp_err=sqrt(sum(fida_bes_err.*time_dex,2,'omitnan').^2)./sum(time_dex,2,'omitnan');
            %             tmp_err = (fida./bes).*dex.*time_dex;
            %             tmp_err(tmp_err==0) = NaN;
            %             tmp_err = std(tmp_err,0,2,'omitnan');%fida_bes_err_out;
            ystr = 'FIDA/BES';
            ylim([0 0.08])
%             t = annotation('textbox',[0.15 0.15 0.3 0.05],'String',['FIDA Int. Range: [', num2str(fida_range),'] nm'],'FitBoxToText','on','LineStyle','none');
%             t.LineWidth = 0.01;
        case 'timetrace_fida'
            plot(ax{i},time(2:end),fida(channel,2:end))
            xlabel(ax{i},'Time [s]')
            ylabel(ax{i},['Integrated FIDA: ', num2str(fida_range(1)), '-',num2str(fida_range(end)), 'nm']);
            legend(ax{i},deblank(names(channel)'))
        case 'timetrace_bes'
            plot(ax{i},time(2:end),fida(channel,2:end))
            xlabel(ax{i},'Time [s]')
            ylabel(ax{i},['Integrated BES: ', num2str(bes_range(channel,1)), '-',num2str(bes_range(channel,end)), 'nm']);
            legend(ax{i},deblank(names(channel)'))
        case 'timetrace_fidabes'
            %plot(ax,time(2:end),fida(channel,2:end)./bes(channel,2:end), 'DisplayName', ['FIDA/BES Chan: ', names{channel}])
            %plot(ax,repmat(time(2:end),1,numel(channel)),movmean(fida(channel,2:end)./bes(channel,2:end),8,2)') %['FIDA/BES smoothed ',
            plot(ax,repmat(time(2:end),1,numel(channel)),(fida(channel,2:end)./bes(channel,2:end))') %['FIDA/BES smoothed ',
            tmp = (fida(channel,:)./bes(channel,:))';
            plot(ax,repmat(time(:),1,numel(channel)),tmp) %['FIDA/BES smoothed ',
            xlabel(ax,'Time [s]')
            ylabel(ax,['Integrated FIDA/BES: ', num2str(fida_range(1)), '-',num2str(fida_range(end)), 'nm']);
            ylim([0,0.16])
            legend(ax,deblank(names(channel)'))
            sname = [filename(1:end-3), '_', plot_type{i},'.txt'];
            fid=fopen(sname,'w');
            fprintf(fid,'%s created on %s\n',plot_type{i}, datestr(now));
            fprintf(fid,'Time %s\n', strjoin(deblank(names(channel)),' '));
            format = repmat('%.3f ',1,size(tmp,2)+1);
            fprintf(fid,[format, '\n'],[time, tmp]');
            %for l = 1:numel(time)
                %fprintf(fid,'%.3f %.3f\n',time(l),tmp(l,:));
            %end
            fclose(fid);
        case 'spectrum'
            time_dex_spec = permute(repmat(time_dex,1,1,size(spec,1)),[3,1,2]);
            tmp =squeeze(sum(spec.*time_dex_spec,3)./sum(time_dex_spec,3));
            plot(ax{i},lambda(:,channel),tmp(:,channel), 'DisplayName',['Data ',  num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);%, Chan: ', names{channel} 
            tmp_err=sqrt(sum(spec_err.*time_dex_spec,3).^2)./sum(time_dex_spec,3);
            %plot(ax{i},lambda(:,channel),squeeze(spec_passive(:,channel,1)),'LineWidth',2.0,'HandleVisibility','off');
            tmp_in =squeeze(sum(spec_in.*time_dex_spec,3)./sum(time_dex_spec,3));
            %plot(ax{i},lambda(:,channel),tmp_in(:,channel),'LineWidth',2.0,'HandleVisibility','off'); 
            %errorbar(ax{i},lambda(:,channel),tmp(:,channel),tmp_err(:,channel),'--','DisplayName',['Avg. ', num2str(t_point), 's'], 'LineWidth',2.0);
            %plot(ax{i},lambda(:,channel),tmp(:,channel), 'DisplayName',['Data ',  num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's, Chan: ', names{channel} ], 'LineWidth',2.0);
            %tmp_in =squeeze(sum(spec_in.*permute(repmat(time_dex,1,1,size(spec_in,1)),[3,1,2]),3)./sum(permute(repmat(time_dex,1,1,size(spec_in,1)),[3,1,2]),3));
            %plot(ax{i},lambda(:,channel),tmp_in(:,channel), 'DisplayName',['Data, bo BG sub ', num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
            xline([bes_range(channel,:)],'HandleVisibility','off');%'BES Range')
            set(ax{i},'YScale','log');
        case 'movie'
            v = VideoWriter([num2str(shotid), '_fidabes_movie','.mp4'], 'MPEG-4');
            open(v);
            for j=1:numel(t)
                plot(ax{i},R_pts(time_dex),fida(time_dex)./bes(time_dex),'o')
                ylim([0 0.24])
                frame = getframe(gcf); %for writing
                writeVideo(v,frame);
            end
            close(v);
    end

    if ~strcmp(plot_type{i},'spectrum') && ~strcmp(plot_type{i}(1:3),'tim')
        plot(ax{i},R_pts(time_dex),tmp_frames,'o','Color', gray,'DisplayName',['Data ', num2str(t_point - avg_time/2),' - ',num2str(t_point + avg_time/2), 's'], 'LineWidth',2.0);
        errorbar(ax{i},R(dex_in)*100, tmp(dex_in),tmp_err(dex_in),'--','DisplayName',['Avg. ', num2str(t_point), 's'], 'LineWidth',2.0);
        xlabel(ax{i},'R [cm]')
        ylabel(ax{i},ystr);
    end
    if lsave
        legend(ax{i});
        sname = [num2str(shotid), '_', num2str(t_point*1000), '_' plot_type{i}];
        savefig(gcf,sname)
        %set(gcf,'renderer','Painters');
        exportgraphics(gcf,[sname,'.eps'],'Resolution',300);
    end
end

plot_data.R = R;
plot_data.spec = spec;
plot_data.names = names;
plot_data.lambda=lambda;
if t_point ~=0
    bes_out = sum(bes.*dex.*time_dex,2,'omitnan')./sum(time_dex,2);
    fida_out =sum(fida.*dex.*time_dex,2,'omitnan')./sum(time_dex,2);
    bes_err_out = sqrt(sum(bes_err.^2.*dex.*time_dex,2,'omitnan').^2 + std(bes_err.*dex.*time_dex,0,2).^2)./sum(time_dex,2);
    fida_err_out = sqrt(sum(fida_err.^2.*dex.*time_dex,2,'omitnan').^2 + std(fida_err.*dex.*time_dex,0,2).^2)./sum(time_dex,2);
    fida_bes_err_out = sqrt(sum(fida_bes_err.^2.*dex.*time_dex,2,'omitnan').^2 + std(fida_bes_err.*dex.*time_dex,0,2).^2)./sum(time_dex,2);
    plot_data.bes = bes_out;
    plot_data.fida = fida_out;
    plot_data.bes_err = bes_err_out;
    plot_data.fida_err = fida_err_out;
    plot_data.fida_bes_err = fida_bes_err_out;
end

sim_data.ax = ax;
sim_data.R = R;
sim_data.names = names;
sim_data.names_unsorted = names_unsorted;
sim_data.lambda=lambda;
sim_data.dispersion = dispersion_in;
sim_data.bes_range = bes_range;
sim_data.fida_range = fida_range;
if ~(strcmp(dex_in,''))
    sim_data.dex = dex_in;
end
sim_data.instfu=instfu;
sim_data.instfu_gamma=instfu_gamma;
sim_data.instfu_box_nm=instfu_box_nm;


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