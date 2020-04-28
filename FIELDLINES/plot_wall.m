function phandle = plot_wall(data,varargin)
%PLOT_WALL wrapper for plot_divertor
%   Wrapper for PLOT_DIVERTOR
%
%  See also plot_divertor.
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           4/28/20

if nargin>1
    phandle = plot_divertor(data,varargin{:});
else
    phandle = plot_divertor(data);
end
return
end

