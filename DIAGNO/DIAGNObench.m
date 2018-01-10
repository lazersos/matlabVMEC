function DIAGNObench(varargin)
%DIAGNObench Plots benchmarks for DIAGNO
%   Plots benchmarks for DIAGNO
%   There exist various options for DIAGNO.  Default is the dual axis plot.
%   Options:
%       'dual'  :   Left Axis (exp) Right Axis (sim)
%       'exp'   :   Experimental Measurements only
%       'same'  :   Plotted on the Same Axis
%       'dphi'  :   Remove offset from simulated results (same axis)
%       'error' :   Extract Errorbars from data
%       'outputname'    :   First section of filename string.
%       'suffix':   Last section of filename.

xmin=-1;
xmax=6;
scale=1;
plottype='dual';
diagno_prefix='diagno_flux.val*';    % DIAGNO Results
data_file='/Volumes/slazerso/Sims/LHD/shots/LHDrecondata_93787.nc';     % Experimental Measurements
outputname='LHD93787_';
suffixname='diagno';
if nargin>0
    for i=1:nargin
        switch varargin{i}
            case {'dual','exp','same','dphi','error'}
                plottype=varargin{i};
            case 'filename'
                i=i+1;
                data_file=varargin{i};
            case 'scale'
                i=i+1;
                scale=varargin{i};
            case 'outputname'
                i=i+1;
                outputname=varargin{i};
            case 'suffix'
                i=i+1;
                suffixname=varargin{i};
        end
    end
end
if strcmp(plottype,'error')
    suffixname=[suffixname '_error'];
end
line_color='black';
line_width=4;
marker_line_size=2;
marker_size=10;
font_size=30;
font_weight='bold';
marker_color='black';
% Import Data
eFlux=read_netcdf(data_file);       % Read experimental data
if ~strcmp(plottype,'exp')
    diagno_scan=dir(diagno_prefix);     % Get List of DIAGNO FILES
    nsfiles=size(diagno_scan,1);        % Get Number of DIAGNO FILES
    diagno_files=diagno_scan(1).name;   % Get First Filename
    st=str2num(diagno_scan(1).name(16:18)); % Create Simulation Time Vector
    for i=2:nsfiles                     % Loop Over Files
        diagno_files=[diagno_files; diagno_scan(i).name];
        st=[st; str2num(diagno_scan(i).name(16:18))];
    end
    st=st./100;                         % Renormalize Simulation Time VectorsFlux=read_diagno_flux(diagno_files(1,:));
    for i=2:nsfiles
        sFlux=read_diagno_flux(diagno_files(i,:),sFlux);
    end
    sFlux.data=sFlux.data.*scale;
end
et=(1:eFlux.tdim).*eFlux.dt+eFlux.t0;% Create experimental Time Vector
% Define the order of the plots Array 1 then Array 2
order=[1 2 3 4 5 6 19 20 21 22 23 24 7 8 9 10 11 12 13 14 15 16 17 18];
% Make Plots
for i=1:4
    fhandle=figure;                     % Get a Figure Window
    phandle=subplot(2,3,1);             % Make Subplot
    if ~strcmp(plottype,'exp')
        firstname=sFlux.names{order((i-1)*6+1),:};
    end
    for j=(i-1)*6+1:(i-1)*6+6       % Plot 6 graphs at a time
        pdex=j-(i-1)*6;             % Index for subplot
        subplot(2,3,pdex);          % Choose proper plot
        plottitle='';
        if strcmp(plottype,'exp')
            plot(et,eFlux.(strcat('FluxLoop',num2str(order(j)-1))),...
                'LineWidth',line_width,...
                'Color',line_color);
        elseif strcmp(plottype,'same') 
            plot(et,eFlux.(strcat('FluxLoop',num2str(order(j)-1))),...
                'LineWidth',line_width,...
                'Color',line_color);
            hold on
            plot(st,sFlux.data(order(j),:),'o',...
                'LineWidth',line_width,...
                'Color',line_color);
            hold off
        elseif strcmp(plottype,'dphi')
            plot(et,eFlux.(strcat('FluxLoop',num2str(order(j)-1))),...
                'LineWidth',line_width,...
                'Color',line_color);
            hold on
            plot(st,(sFlux.data(order(j),:)-sFlux.data(order(j),1)),'o',...
                'LineWidth',line_width,...
                'Color',line_color);
            hold off
        elseif strcmp(plottype,'dual')
            [AX, H1, H2]=plotyy(et,eFlux.(strcat('FluxLoop',num2str(order(j)-1))),...
                st,sFlux.data(order(j),:));
            set(H1,'LineWidth',line_width,...
                'Color',line_color);
            set(AX(1),'YColor','black');
            set(AX(2),'YColor','black','FontSize',font_size-10,...
                'FontWeight',font_weight);
            set(H2,'Linestyle','none','Marker','o','Color','black');
            xlim(AX(1),[xmin xmax]);
            xlim(AX(2),[xmin xmax]);
        end
        xlim([xmin xmax]);
        set(gca,'FontSize',font_size-10,'FontWeight',font_weight);
        xlabel('Time [s]','FontSize',font_size,'FontWeight',font_weight,'Color','black');
        ylabel('Flux [Wb]','FontSize',font_size,'FontWeight',font_weight,'Color','black');
        if ~strcmp(plottype,'exp')
            plottitle=['Saddle Coil ',num2str(sFlux.names{order(j),:}),plottitle];
        end
        if strcmp(plottype,'error')
            sdata=sFlux.data(order(j),:);
            edata=eFlux.(strcat('FluxLoop',num2str(order(j)-1)));
            dex=max(find(et<=st(1)));
            for k=2:max(size(sdata))
                dex=[dex max(find(et<=st(k)))];
            end
            % Now we just have the simulated datapoints
            edata=edata(dex);
            et=et(dex);
            % We convert experimental data to simulated via initial offset
            dedata=edata(:)+sdata(1);
            % Get the Error
            error=(dedata'-sdata)./sdata;
            disp([plottitle ' sigma:' num2str(std(error))]);
            plot(et,error*100.,'o','LineWidth',line_width,'Color',line_color);
            ylabel('Error [%]','FontSize',font_size,'FontWeight',font_weight,...
                'Color','black');
            % Now we generate some statistics
            pause(0.1);
        end
        title(plottitle,'FontSize',font_size,'FontWeight',font_weight);
    end
    %saveas(fhandle,[outputname firstname '_' sFlux.names{order(j),:} '_' suffixname '.fig']);
    pause(0.1);
end

end

