function pcmdata=coil_poincm(coildata,varargin)
%COIL_POINCM Creates Poinciare Map using biot-savart


% Set Defaults
sr=0:.025:0.5;  % Starting Points for map
rmaj=5.75;   % Offset for starting points
dt=0.01;     % Timestep
extcur=1.45e6.*[1 1 1 1 1 0 0 0];
phistart=0.0;
nfp=5;
maxnumhits=100;
rhomax=2.0;
cx=(sr+rmaj)*cos(phistart);
cy=(sr+rmaj)*sin(phistart);
cz=0.0*cx;
plotprogress=0;
% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin > numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'plotprogress'}
                plotprogress=1;
            case {'extcur'}
                i=i+1;
                extcur=varargin{i};
            case {'cx'}
                i=i+1;
                cx=varargin{i};
            case {'cy'}
                i=i+1;
                cy=varargin{i};
            case {'cz'}
                i=i+1;
                cz=varargin{i};
            case {'rmaj'}
                i=i+1;
                rmaj=varargin{i};
            case {'dt'}
                i=i+1;
                dt=varargin{i};
            case {'nfp'}
                i=i+1;
                nfp=varargin{i};
            case {'phistart'}
                i=i+1;
                phistart=varargin{i};
            case {'maxnumhits'}
                i=i+1;
                maxnumhits=varargin{i};
        end
    end
end
% Create arrays
philim=phistart+2*pi/nfp;
philinexy(1,:)=(0:(max(sr)+rmaj)/25:(max(sr)+rmaj))*cos(philim);
philinexy(2,:)=(0:(max(sr)+rmaj)/25:(max(sr)+rmaj))*sin(philim);
pcmdata=zeros(max(size(sr)),maxnumhits,2);   % Store the poiniare points
pcmdata(:,1,1)=cx(:).*cos(phistart);
pcmdata(:,1,2)=cz(:);
% Plotting info
xmin=0; xmax=7; ymin=-7; ymax=0; zmin=-2; zmax=2;
% Setup Initial Positions
bx=0.0*cx;
by=0.0*cy;
bz=0.0*cz;
modb=0.0*cz;
linside=0.0*cz;
numhits=0.0*cz;
maxiter=1000;
sx=0.0*cx;
sy=0.0*cy;
sz=0.0*cz;
r=0.0*cx;
z=0.0*cy;
rho=0.0*cz;
theta=0.0*cz;
if plotprogress
    xlim([0 7]);
    ylim([0 7]);
    zlim([-1.5 1.5]);
end
% Read Coil file
%coildata=read_coils('coils.w7x_np10');
% Setup for biot_savart calculation
coildata=coil_biot_prep(coildata);
% Loop over starting points
for i=1:numel(cx)
    linside(i)=1;
    numhits(i)=1;
    if plotprogress
        clf
        subplot(2,2,1)
        hold off
        subplot(2,2,2)
        hold off
        subplot(2,2,3)
        hold off
        subplot(2,2,4)
        hold off
    end
    disp(['  -- Launching Fieldline: ' num2str(i)]);
    pause(0.01);
    tic;
    n=0;
    while linside(i)
        n=n+1;
        % Get B-Field
        [bx(i) by(i) bz(i)]=coil_biot(coildata,cx(i),cy(i),cz(i),extcur);
        % Calculate s
        modb(i)=sqrt(bx(i)*bx(i)+by(i)*by(i)+bz(i)*bz(i));
        sx(i)=bx(i)./modb(i);
        sy(i)=by(i)./modb(i);
        sz(i)=bz(i)./modb(i);
        % Move particle
        cx(i)=cx(i)+sx(i)*dt;
        cy(i)=cy(i)+sy(i)*dt;
        cz(i)=cz(i)+sz(i)*dt;
        % Handle rho>1.5;
        r(i)=sqrt(cx(i).*cx(i)+cy(i).*cy(i))-rmaj;
        z(i)=cz(i);
        rho(i)=sqrt(r(i).*r(i)+z(i).*z(i));
        %if rho(i) > rhomax
        %    linside(i)=0;
        %    disp(['     Point: ' num2str(i) '/' num2str(numel(cx))...
        %        ' Out of Domain! @ numhits:' num2str(numhits(i))]);
        %end
        if numhits(i) > maxnumhits
            linside(i)=0;
        end
        % Handle particle hitting a field period
        theta(i)=atan(cy(i)/cx(i));
        if plotprogress && (linside(i) > 0)
            %Plot Progress
            % X-Y plot
            subplot(2,2,2)
            plot(cx(i),cy(i),'.');
            xlim([xmin xmax]);
            ylim([ymin ymax]);
            hold on
            plot(philinexy(1,:),philinexy(2,:),'red');
            hold off
            title('X-Y');
            % X-Z Plot
            subplot(2,2,3)
            plot(cx(i),cz(i),'.');
            xlim([xmin xmax]);
            ylim([zmin zmax]);
            hold on
            title('X-Z');
            % R-Z plot
            subplot(2,2,1)
            plot(r(i)+rmaj,cz(i),'.');
            hold on
            plot(pcmdata(i,1:numhits(i),1),pcmdata(i,1:numhits(i),2),'.r');
            hold off
            xlim([zmin+rmaj zmax+rmaj]);
            ylim([zmin zmax]);
            axis square
            hold on
            title('R-Z');
            % 3d PLOTS
            subplot(2,2,4);
            %plot x-y
            plot3(cx(i),cy(i),zmin,'.r')
            hold on
            plot3([cx(i) cx(i)],[cy(i) cy(i)],[zmin cz(i)],'b')
            %plot y-z
            plot3(xmax,cy(i),cz(i),'.r')
            plot3([xmax cx(i)],[cy(i) cy(i)],[cz(i) cz(i)],'b')
            %plot x-z
            plot3(cx(i),ymax,cz(i),'.r')
            plot3([cx(i) cx(i)],[ymax cy(i)],[cz(i) cz(i)],'b')
            %plot3D
            plot3(cx(i),cy(i),cz(i),'.')
            grid on
            hold off
            xlim([xmin xmax]);
            ylim([ymin ymax]);
            zlim([zmin zmax]);
            pause(.0001);
        end
        if (n > maxiter)
            linside(i)=0;
        end
        if theta(i)>philim
            %Store the location
            pcmdata(i,numhits(i),1)=r(i)+rmaj;
            pcmdata(i,numhits(i),2)=z(i);
            %plot(r(i)+rmaj,z(i),'.b');
            %title(['Poincm (i=' num2str(i) ' n=' num2str(numhits(i)) ')']);
            cx(i)=(r(i)+rmaj)*cos(phistart);
            cy(i)=(r(i)+rmaj)*sin(phistart);
            disp(['     Point: ' num2str(i) '/' num2str(numel(cx))...
                '   Hit: ' num2str(numhits(i))]);
            numhits(i)=numhits(i)+1;
        end
    end
    disp(['     Took: ' num2str(toc) ' [s]']);
end
return
end

