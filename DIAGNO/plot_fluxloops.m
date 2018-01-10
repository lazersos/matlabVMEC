function plot_fluxloops(data,varargin)
%PLOT_FLUXLOOPS(loop_data[,'loopname'...]) Plots the DIAGNO fluxloops file.
%   This function plots the data read by the READ_FLUXLOOPS program.  This
%   function accepts any argument which is accepted by the PLOT3 routine.
%
%   Options:
%       READ_FLUXLOOPS(data,'loopname')  This will display the loop name
%       at the geometrical center of each loop.
%   
%   Example:
%       loop_data=read_fluxloops('diagno_fluxloops.machine','scale',0.001);
%       plot_fluxloops(loop_data,'loopname','Color','black',LineWidth',4);
%
%   See also read_fluxloops, plot3.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/04/11

%Set Some defaults
loopnames=0;
l2d = 0;
plot3argin=varargin;
dex = zeros(1,length(varargin));
%Handle Input arguments
offset=1;
if nargin>offset
    i=1;
    while (i <= nargin-offset)
        switch varargin{i}
            case 'loopname'
                loopnames=1;
            case '2D'
                l2d=1;
            otherwise
                dex(i) = 1;
        end
        i = i + 1;
    end
end
plot3argin = varargin(dex==1);
axis vis3d;
hold on
%   Plots Saddle Coils
for i=1:data.nloops
    xl=data.loops{i}(1,:);
    yl=data.loops{i}(2,:);
    zl=data.loops{i}(3,:);
    if (l2d)
        plot(xl(1), zl(1),'s',plot3argin{:});
        if loopnames
            text(xl(1)+0.01,zl(1)-0.01,data.loopname{i}','Color','black');
        end
    else
        plot3(xl,yl,zl,plot3argin{:});
        if loopnames
            x=mean(data.loops{i}(1,:));
            y=mean(data.loops{i}(2,:));
            z=mean(data.loops{i}(3,:));
            text(x,y,z,data.loopname{i}','Color','black');
        end
    end
end
hold off
end

