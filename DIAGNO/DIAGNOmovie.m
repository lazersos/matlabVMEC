function DIAGNOmovie(varargin)
%DIAGNOMOVIE Creates a movie from DIAGNO outputs
%   Detailed explanation goes here

%Setup Defaults
plottype='surf';
ntheta=180;
nzeta=180;
suffix='*';
fluxloopfile='';
fig_color='black';
title_color='white';
title_fontsize=20;
axis_color='white';
axis_bgcolor='black';
axis_fontsize=20;
%Handle Input arguments
if nargin>0
    for i=1:nargin
        switch varargin{i}
            case {'surf'}
                plottype=varargin{i};
            case 'suffix'
                i=i+1;
                data_file=varargin{i};
            case 'ntheta'
                i=i+1;
                ntheta=varargin{i};
            case 'outputname'
                i=i+1;
                outputname=varargin{i};
            case 'nzeta'
                i=i+1;
                nzeta=varargin{i};
            case 'fluxloops'
                i=i+1;
                fluxloopfile=varargin{i};
        end
    end
end
% Get list of filenames
filelist=dir(['diagno_in.' suffix]);
nsfiles=size(filelist,1);
filenames=filelist(1).name;
for i=2:nsfiles
    filenames=[filenames; filelist(i).name];
end
% Setup Figure Environment
fig=figure('DoubleBuffer','on','Position',[1 1 1920 1080],...
    'Color',fig_color,'BackingStore','on','MenuBar','none',...
    'InvertHardcopy','off');
set(gca,'nextplot','replacechildren','XColor',axis_color,...
    'YColor',axis_color,'ZColor',axis_color,'Color',axis_bgcolor,...
    'FontSize',axis_fontsize);
colormap jet;
% Cycle through surfaces
for i=1:nsfiles
    cla;
    % Plot Diagnostics
    if ~isempty(fluxloopfile)
        plot_fluxloops(fluxloopfile,'color','blue');
    end
    hpatch=DIAGNOplot(filenames(i,:),'nzeta',nzeta,'ntheta',ntheta);
    title(filenames(i,:),'Color',title_color,'FontSize',title_fontsize);
    pause(.1);
end

end

