function [ output_args ] = LHDplot_current(data,varargin)
%LHDPLOT_CURRENT(data) Creates current plots from data.
%   No Detailed explanation yet

shotnum=0000;
line_width=4;
size=[1920 1080];
font_size=30;
font_weight='bold';
% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'SHOT','shot'};
                i=i+1;
                shotnum=varargin{i};
        end
    end
end
%Setup Figure and time axis
fig=figure('Position',[0 0 size],'Units','pixels');
t=(0:data.tdim-1).*data.dt+data.t0;
%Plot
plot(t,data.CurrentHCI2,'LineWidth',line_width,'Color','yellow')
hold on
plot(t,data.CurrentHCM0,'LineWidth',line_width,'Color','blue')
plot(t,data.CurrentHCO1,'LineWidth',line_width,'Color','green')
plot(t,data.CurrentPOV3,'LineWidth',line_width,'Color','cyan')
plot(t,data.CurrentPIS4,'LineWidth',line_width,'Color','red')
plot(t,data.CurrentPIV5,'LineWidth',line_width,'Color','magenta')
axis tight
vals=[max(data.CurrentHCI2) max(data.CurrentHCM0) max(data.CurrentHCO1)...
      max(data.CurrentPOV3) max(data.CurrentPIS4) max(data.CurrentPIV5)];
ylim([0 max(ceil(vals))]);
hold off
ahandle=gca;
title(['LHD Coil Currents (shot:' num2str(shotnum) ')'],...
    'FontSize',font_size,'FontWeight',font_weight);
xlabel('Time [s]',...
    'FontSize',font_size,'FontWeight',font_weight);
ylabel('Current [kA]',...
    'FontSize',font_size,'FontWeight',font_weight);
lhandle=legend('I_{HCI}','I_{HCM}','I_{HCO}','I_{POV}','I_{PIS}','I_{PIV}');
set(ahandle,'FontSize',font_size,'FontWeight',font_weight);
set(lhandle,'Location','NorthEastOutside');
end

