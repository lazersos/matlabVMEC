function pixplot(varargin)
%PIXPLOT([x,y],value) Pseudocolor (checkerboard) plot (on centered grid)
%   PIXPLOT creates a pseudocolor (checkerboard) plot on a centered grid.
%   This addresses the issue with pcolor where the last column and row are
%   clipped and ticklabels are along grid lines rather than centered on the
%   pixels.
%
%   X and Y are optional arrays (which may be vectors) specifying the x and
%   y value for each datapoint in VALUE.  If they are not supplied, axis
%   labels are assumed to be the indexes of VALUE.
%
%   Example
%     val=magic(100); % Create a 2D array to plot
%     pixplot(val); % Plot the array
%
%   Version 1.05
%   Maintained by: Samuel Lazerson (lazerson@pppl.gov)
%   Date  09/27/2010

numticks=7; % Number of tickmarks
if nargin==1
    val=varargin{1};
    % Create Axes
    nx=size(val,1);
    ny=size(val,2);
    x=repmat((1:nx)',[1 ny]);
    y=repmat((1:ny),[nx 1]);
elseif nargin==3
    x=varargin{1};
    y=varargin{2};
    val=varargin{3};
else
    disp('ERROR: Incorrect number of arguments!');
    return
end
% We need to add a dimension to the arrays
nx=size(val,1);
ny=size(val,2);
maxx=max(max(x));
maxy=max(max(y));
maxval=max(max(val));
newx=ones(nx+1,ny+1).*maxx;
newy=ones(nx+1,ny+1).*maxy;
newval=ones(nx+1,ny+1).*maxval;
% We need to handle x and y carefully.
if size(x,1) == 1
    newx(1:nx,1)=x(1,1:nx)';
    newx=repmat(newx(1:nx+1,1),[1 ny+1]);
elseif size(x,2) == 1
    newx(1:nx,1)=x(1:nx,1);
    newx=repmat(newx(1:nx+1,1),[1 ny+1]);
elseif ((size(x,1) == nx) && (size(x,2) == ny))
    newx(1:nx,1:ny)=x;
else
    disp('ERROR: Size mismatch between x and val!');
    return;
end
if size(y,1) == 1
    newy(1,1:ny)=y(1,1:ny);
    newy=repmat(newy(1,1:ny+1),[nx+1 1]);
elseif size(y,2) == 1
    newy(1,1:ny)=y(1:ny,1)';
    newy=repmat(newy(1,1:ny+1),[nx+1 1]);
elseif ((size(y,1) == nx) && (size(y,2) == ny))
    newy(1:nx,1:ny)=y;
else
    disp('ERROR: Size mismatch between y and val!');
    return;
end
newval(1:nx,1:ny)=val;
% The X and Y axis must have values that make sense
for i=1:nx
    newx(i,ny+1)=2.*newx(i,ny)-newx(i,ny-1);
    newy(i,ny+1)=2.*newy(i,ny)-newy(i,ny-1);
end
for j=1:ny
    newx(nx+1,j)=2.*newx(nx,j)-newx(nx-1,j);
    newy(nx+1,j)=2.*newy(nx,j)-newy(nx-1,j);
end
newx(nx+1,ny+1)=2.*newx(nx,ny+1)-newx(nx-1,ny+1);
newy(nx+1,ny+1)=2.*newy(nx+1,ny)-newy(nx+1,ny-1);
hpcolor=pcolor(newx,newy,newval);
set(hpcolor,'EdgeColor','none'); % Get rid of grids
% Now we create the axes
% We assume the axes are not equidistant but are cartesian
xticks=zeros(1,nx);
yticks=zeros(1,ny);
for i=1:nx-1
        xticks(i)=newx(i,1)+(newx(i+1,1)-newx(i,1))*.5;
end
xticks(nx)=newx(nx,1)+(newx(nx,1)-newx(nx-1,1))*.5;
for j=1:ny-1
    yticks(j)=newy(1,j)+(newy(1,j+1)-newy(1,j))*.5;
end
yticks(ny)=newy(1,ny)+(newy(1,ny)-newy(1,ny-1))*.5;
if nx > numticks
    xstep=round(nx/(numticks));
    xtic=xticks(1:xstep:nx);
    xlab=newx(1:xstep:nx,1);
else
    xtic=xticks;
    xlab=newx(1:nx,1);
end
% Set Y Axis Ticks
if ny > numticks
    ystep=round(ny/(numticks));
    ytic=yticks(1:ystep:ny);
    ylab=newy(1,1:ystep:ny);
else
    ytic=yticks;
    ylab=newy(1,1:ny);
end
set(gca,'XTick',xtic,'XTickLabel',round(xlab,2));
set(gca,'YTick',ytic,'YTickLabel',round(ylab,2));
end

