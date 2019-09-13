function fig = make_focus_movie( data, varargin)
%MAKE_FOCUS_MOVIE(data) Creates a summary animiation of an optimization
%   MAKE_FOCUS_MOVIE(data)  Creates a summary animiation of an optimziation
%   as produced by the Focus code.  It takes the data structure returned by
%   the READ_FOCUS routine.
%
%   Usage:
%   focus_data=read_focus('test.fo.h5');
%   make_focus_movie(focus_data);
%   make_focus_movie(focus_data,'surface'); % displays plasma surface
%   make_focus_movie(focus_data,'movie'); % Save as avi file
%   make_focus_movie(focus_data,'final'); % only plot initial and final
%
%   See also read_focus.
%
%   Written by:     S. Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           02/08/17


%Set some defaults
lsurface=0;
lmovie=0;
movie_frames=1:size(data.coilspace,1);


% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'surface'}
                lsurface=1;
            case {'movie'}
                lmovie=1;
            case {'final'}
                movie_frames=[1 size(data.coilspace,1)];
        end
    end
end

% Setup Movie File
if (lmovie)
    mov=VideoWriter([data.file_ext '_movie.avi']);
    mov.Quality= 100;
    mov.FrameRate=1;
    open(mov);
end

% Setup the harmonic data
mn    = data.NFcoil+1;
nharm = (mn*6)-3+1; % includes current
theta=0:2*pi/128:2*pi;
%theta=0:2*pi/double(data.NDcoil-1):2*pi;
for i=1:mn
    cth(:,i) = cos(theta'.*double(i-1));
end
for i=2:mn
    sth(:,i-1) = sin(theta'.*double(i-1));
end

% Indexing for coil
dxc = [2 5:6:nharm];
dxs = 8:6:nharm;
dyc = [3 6:6:nharm];
dys = 9:6:nharm;
dzc = [4 7:6:nharm];
dzs = 10:6:nharm;

% Create the figure
fig=figure('Position',[1 1 1024 768],'Color','white');
% Plot surface
if lsurface
    subplot(2,2,[1 3]);
    temp=size(data.rsurf);
    e1=1:temp(1); e2 = 1:temp(2);
    if (temp(1) > 64)
        e1=1:round(temp(1)/(64-1)):temp(1);
    end
    if (temp(2) > 64)
        e2=1:round(temp(2)/(64-1)):temp(2);
    end
    r = zeros(1,length(e1),length(e2));
    r(1,:,:) = data.rsurf(e1,e2);
    z = zeros(1,length(e1),length(e2));
    z(1,:,:) = data.zsurf(e1,e2);
    if isfield(data,'curBn')
        bn = zeros(1,length(e1),length(e2));
        bn(1,:,:) = data.curBn(e1,e2);
        isotoro(r,z,data.psurf(1,e2),1,bn);
    else
        isotoro(r,z,data.psurf(1,e2),1);
        camlight left;
    end
    hold on;
end
% Plot the functions
subplot(2,2,2);
plot(data.evolution(:,1),data.evolution(:,2),'k');
hold on;
plot(data.evolution(:,1),data.evolution(:,3),'Color',[0.5 0.5 0.5]);
axis tight;
ylim([min(ylim) 3]);
set(gca,'FontSize',24);
xlabel('Iteration');
ylabel('F, dF/dx');
legend('Function','Derivative');
hred=plot([1 1].*1,ylim,'r');
if (abs(data.evolution(end,2)-data.evolution(1,2)) > 10.), set(gca,'YScale','log');end
hold off;
hbar=[];
%Loop over coils
for i=movie_frames
    subplot(2,2,[1 3]);
    % Extract the coil
    farr=data.coilspace(i,:);
    for j=1:data.Ncoils
        d1=(nharm)*(j-1)+1;
        d2=d1+nharm-1;
        carr=farr(d1:d2);
        current(i,j) = carr(1);
        % Format xc0 yc0 zc0 xc1 xs1 y
        x = cth*carr(dxc)'+sth*carr(dxs)';
        y = cth*carr(dyc)'+sth*carr(dys)';
        z = cth*carr(dzc)'+sth*carr(dzs)';
        hold on;
        lcolor='k';
        if i==1, lcolor='r'; end;
        htemp(j)=plot3(x,y,z,lcolor);
        hold off;
    end
    axis equal;
    view(3);
    subplot(2,2,2);
    delete(hred);
    hold on; hred=plot([1 1].*i,ylim,'r'); hold off;
    subplot(2,2,4);
    cla;
    if i >1
        bar(current([1 i],:)','grouped')
        xlim([0 data.Ncoils+1]);
        set(gca,'FontSize',24);
        xlabel('Current Group');
        ylabel('Current');
    end
    if (i>1 && i<size(data.coilspace,1))
        for j=1:data.Ncoils
            set(htemp(j),'Color',[0.5 0.5 0.5]);
        end
    end
    pause(0.1);
    if (lmovie)
        frame=getframe(fig);
        writeVideo(mov,frame);
    end
    
end

if (lmovie)
    close(mov);
end

end

