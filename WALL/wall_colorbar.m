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
caxis([0 1E6]);
set(ha,'FontSize',24,'Color','white');
ytick=0:250E3:1E6;
yticklable={};
for i=1:length(ytick)
    yticklabel{i} = num2str(ytick(i)./1E6,'%3.2f');
end
set(ha,'YTick',ytick,'YTickLabel',yticklabel);
ylabel(ha,'Heat Flux [MW/m^2]');
end

