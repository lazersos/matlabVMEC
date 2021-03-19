function ha = wall_colorbar
%WALL_COLORBAR Plots a wall heat flux colorbar
%   The WALL_COLORBAR function is a wrapper for quickly setting a nice
%   colorbar when plotting wall elements.
%
% Example usage
%      ha=wall_colorbar;
%
% Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
% Version:       1.00

ha=colorbar;
cmax = max(caxis);
scale = 1.0;
units = '[W/m^2]';
if cmax>1E6
    scale = 1E-6;
    units = '[MW/m^2]';
elseif cmax > 1E3
    scale = 1E-3;
    units = '[kW/m^2]';
end
set(ha,'FontSize',24,'Color','white');
ymax = max(ha.Ticks);
ytick = (0:0.25:1).*ymax;

yticklabel={};
for i=1:length(ytick)
    yticklabel{i} = num2str(round(ytick(i).*scale),'%4i');
end
set(ha,'YTick',ytick,'YTickLabel',yticklabel);
ylabel(ha,['Heat Flux ' units]);
end

