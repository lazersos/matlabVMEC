function plot_coils(coildata,varargin)
%PLOT_COILS(coildata[,'simple'/'plane'/'other']) Plots coils
%   PLOT_COILS plots the data read by READ_COILS from a VMEC XGRID file.
%
%   Usage:
%   coil_data=read_coils('coils.test');
%   plot_coils(coil_data);  % Plots each fillament
%   plot_coils(coil_data,'plane');  %Plots intersection of coils (phi=0)
%   plot_coils(coil_data,'simple'); %Plots volumetric rendering of groups
%   plot_coils(coil_data,'field_period'); %Plots only for field period
%   plot_coils(coil_data,'tube'); %Plots coils as volumetric tubes.
%   plot_coils(coil_data,'tube','nedges',20,'tubewidth',0.2); %Plots coils 
%       as volumetric tubes of width tubewidth and having nedges edges.
%       Defaults are 4 edges (square) and tubewidth 0.15.
%   
%
%   Note: The simple and plane option only work properly if all fillaments
%   pass through the phi=0 plane (helical system like LHD).
%
%   See also read_coils.
%
%   Written by:     S. Lazerson (lazerson@pppl.gov)
%   Version:        1.5
%   Date:           12/21/10


%%%%Set Defaults%%%%%%%%%%
calcdex=1;
machine_string='';
plottype='full';
colors{1}=[1 1 0];
colors{2}=[1 0 0];
colors{3}=[0 0 1];
colors{4}=[0 1 1];
colors{5}=[0 1 0];
colors{6}=[.5 .5 .1];
colors{7}=[.5 .1 .1];
colors{8}=[.1 .1 .5];
colors{9}=[.1 .5 .5];
colors{10}=[.1 .5 .1];
lfield_period=0;
nedges=4;
tubewidth=0.15;
for i=11:100
    colors{i}='black';
end
bgcolor='white';
pcolor=[];
phi=0.0;
% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'simple','plane','tube','plane_xz','plane_name'}
                plottype=varargin{i};
            case {'machine'}
                i=i+1;
                machine_string=varargin{i};
            case {'color'}
                i=i+1;
                bgcolor=varargin{i};
            case {'pcolor'}
                i=i+1;
                pcolor=varargin{i};
            case {'phi'}
                i=i+1;
                phi=varargin{i};
            case {'nedges'}
                i=i+1;
                nedges=varargin{i};
            case {'tubewidth'}
                i=i+1;
                tubewidth=varargin{i};
            case {'field_period'}
                lfield_period=1;
        end
    end
end
% Handle Simplification
coildata_save=coildata;
if (lfield_period)
    nfp=coildata.periods;
    ncoils=max(coildata.vert(5,:));
    vert=[];
    name={};
    j=1;
    for i=1:ncoils
        dex = coildata.vert(5,:)==i;
        x=coildata.vert(1,dex);
        y=coildata.vert(2,dex);
        z=coildata.vert(3,dex);
        c=coildata.vert(4,dex);
        d=coildata.vert(5,dex);
        p=atan2(y,x);
        if (mean(p) < 2*pi/nfp && mean(p) >= 0.0)
            vert=[vert [x; y; z; c; d.*0.0+j]];
            name=[name; coildata.current_name{i}];
            j=j+1;
        end
    end
    coildata.vert=vert;
    coildata.current_name=name;
end
data=coildata.vert;
% We should sort the Array by Coil Type
data=sortrows(data',5)';
% Handle Coloring
if calcdex == 1
    disp(' - System Defined Coloring');
    numgroups=max(data(5,:));
    base={[1 1 0] [0 0 1] [0 1 0] [0 1 1] [1 0 0] [1 0 1]};
    if numgroups > 6
        numrep=round(numgroups/6);
        cdex=base;
        for i=1:numrep
            cdex=[cdex base];
        end
        cdex=[cdex base{1:(numgroups-numrep*6)}];
    else
        cdex=base(1:numgroups);
    end
end
%  Now Plot the Coils
%fig=figure('Position',[0 0 1920 1080],'Color',bgcolor);
axis([-5 5 -5 5 -5 5]);
ahandle=gca;
set(ahandle,'Color',bgcolor);
hold on;
firstel=1;
if strcmp(plottype,'simple')    %Really only applies to LHD
    %First we wish to find the coils which define the outter verticies
    %this can be done entirly in the zeta=0 plane.  This must be done
    %for each group
    for i=1:max(data(5,:)) %Each group
        for k=1:1 % 1:1 (most coils), 1:2 (LHD)
            groupdata=data(:,data(5,:)==i);                     % Subsample the data so we only work with one group at a time.
            linedex=find(squeeze(groupdata(4,:))==0)+1;   % Get the ending points for each line
            linedex=[1 linedex];                                % Index for start of each line in group
            linedex=linedex(1:max(size(linedex))-1);            % Remove the last point
            grouplength=find(squeeze(groupdata(4,:))==0);
            grouplength=grouplength(1);
            temp=[];
            for l=1:max(size(linedex))
                if k==1 && i<4                                      %Helical Coils
                    if groupdata(1,linedex(l))<4.0
                        temp=[temp groupdata(:,linedex(l):linedex(l)+grouplength-1)];
                    end
                elseif k==2 && i<4
                    if groupdata(1,linedex(l))>=4.0
                        temp=[temp groupdata(:,linedex(l):linedex(l)+grouplength-1)];
                    end
                elseif k==1 && i>=4                                  %Poloidal Coils
                    if groupdata(3,linedex(l))<0.0
                        temp=[temp groupdata(:,linedex(l):linedex(l)+grouplength-1)];
                    end
                else
                    if groupdata(3,linedex(l))>=0.0
                        temp=[temp groupdata(:,linedex(l):linedex(l)+grouplength-1)];
                    end
                end
            end
            groupdata=temp;
            linedex=find(squeeze(groupdata(4,:))==0)+1;   % Get the ending points for each line
            linedex=[1 linedex];                                % Index for start of each line in group
            linedex=linedex(1:max(size(linedex))-1);            % Remove the last point
            % Set initial guess equal to first line
            xxzx=linedex(1);
            xxzm=linedex(1);
            xmzx=linedex(1);
            xmzm=linedex(1);
            for j=2:max(size(linedex))        % For each line
                %Check for xmin zmin
                if (groupdata(1,linedex(j)) <= groupdata(1,xmzm)) && (groupdata(3,linedex(j)) <= groupdata(3,xmzm))
                    xmzm=linedex(j);
                end
                %Check for xmax zmin
                if (groupdata(1,linedex(j)) >= groupdata(1,xxzm)) && (groupdata(3,linedex(j)) <= groupdata(3,xxzm))
                    xxzm=linedex(j);
                end
                %Check for xmin zmax
                if (groupdata(1,linedex(j)) <= groupdata(1,xmzx)) && (groupdata(3,linedex(j)) >= groupdata(3,xmzx))
                    xmzx=linedex(j);
                end
                %Check for xmax zmax
                if (groupdata(1,linedex(j)) >= groupdata(1,xxzx)) && (groupdata(3,linedex(j)) >= groupdata(3,xxzx))
                    xxzx=linedex(j);
                end
            end
            % Now we have the indicies for each box
            % We must create a patch surface out of them.
            % Create the first eight vertices
            verts=[groupdata(1,xmzm) groupdata(2,xmzm) groupdata(3,xmzm)];
            verts=[verts; groupdata(1,xxzm) groupdata(2,xxzm) groupdata(3,xxzm)];
            verts=[verts; groupdata(1,xxzx) groupdata(2,xxzx) groupdata(3,xxzx)];
            verts=[verts; groupdata(1,xmzx) groupdata(2,xmzx) groupdata(3,xmzx)];
            verts=[verts; groupdata(1,xmzm+1) groupdata(2,xmzm+1) groupdata(3,xmzm+1)];
            verts=[verts; groupdata(1,xxzm+1) groupdata(2,xxzm+1) groupdata(3,xxzm+1)];
            verts=[verts; groupdata(1,xxzx+1) groupdata(2,xxzx+1) groupdata(3,xxzx+1)];
            verts=[verts; groupdata(1,xmzx+1) groupdata(2,xmzx+1) groupdata(3,xmzx+1)];
            % Now construct the first four faces
            faces=[1 2 6 5];
            faces=[faces; 2 6 7 3];
            faces=[faces; 4 3 7 8];
            faces=[faces; 1 4 8 5];
            for j=2:grouplength-1
                verts=[verts; groupdata(1,xmzm+j) groupdata(2,xmzm+j) groupdata(3,xmzm+j)];
                verts=[verts; groupdata(1,xxzm+j) groupdata(2,xxzm+j) groupdata(3,xxzm+j)];
                verts=[verts; groupdata(1,xxzx+j) groupdata(2,xxzx+j) groupdata(3,xxzx+j)];
                verts=[verts; groupdata(1,xmzx+j) groupdata(2,xmzx+j) groupdata(3,xmzx+j)];
                % Now create four faces
                l=(j-2)*4+1;
                faces=[faces; l l+1 l+5 l+4];
                faces=[faces; l+1 l+5 l+6 l+2];
                faces=[faces; l+3 l+2 l+6 l+7];
                faces=[faces; l l+3 l+7 l+4];
            end
            % Now add last last face
            faces=[faces; l+4 l+5 2 1];
            faces=[faces; l+5 2 3 l+6];
            faces=[faces; l+7 l+6 3 4];
            faces=[faces; l+4 l+7 4 1];
            patch('Vertices',verts,'Faces',faces,'FaceColor',cdex{i},...
                'EdgeColor','none');
        end
    end
    %legend(leg);
    camlight left;
    lighting phong;
    view(3);
elseif strcmp(plottype,'plane')
    dex = find(data(4,:)==0);
    dex = [1 dex];
    ndex = length(dex);
    temp_phi = atan2(data(2,:),data(1,:));
    temp_phi(temp_phi < 0) = temp_phi(temp_phi < 0) + 2*pi;
    temp_r   = sqrt(data(1,:).*data(1,:)+data(2,:).*data(2,:));
    temp_z   = data(3,:);
    phi_targ=phi;
    for i=1:ndex-1
        d1=dex(i);
        d2=dex(i+1);
        for j=d1:d2-1
            if ((temp_phi(j) < phi_targ) && (temp_phi(j+1) > phi_targ)) ...
                    || (temp_phi(j) == phi_targ) || (temp_phi(j+1) == phi_targ) ...
                    || ((temp_phi(j) > phi_targ) && (temp_phi(j+1) < phi_targ))
                r=pchip(temp_phi(j:j+1),temp_r(j:j+1),phi_targ);
                z=pchip(temp_phi(j:j+1),temp_z(j:j+1),phi_targ);
                hold on;
                plot(r,z,'sk');
                hold off;
            end
        end
    end
    view(2)
elseif strcmp(plottype,'plane_name')
    temp_phi=atan(data(2,:)./data(1,:));
    for i=1:max(data(5,:)) %Each group
        groupdata=data(:,data(5,:)==i);                     % Subsample the data so we only work with one group at a time.
        linedex=find(squeeze(groupdata(4,:))==0);   % Get the ending points for each line
        linedex=[1 linedex];                                % Index for start of each line in group
        %linedex=linedex(1:max(size(linedex))-1);            % Remove the last point
        % Now find phi elements
        %         x=zeros(1,numel(linedex)-1);
        %         y=zeros(1,numel(linedex)-1);
        %         z=zeros(1,numel(linedex)-1);
        %         r=zeros(1,numel(linedex)-1);
        %         for j=1:numel(linedex)-1;
        %             phidex=find(temp_phi(linedex(j):linedex(j+1)-2)>phi,1,'first')+linedex(j)-1;
        %             if temp_phi(phidex)==phi;
        %                 x(j)=data(1,phidex);
        %                 y(j)=data(2,phidex);
        %                 z(j)=data(3,phidex);
        %             else
        %                 dphi=(temp_phi(phidex)-phi)/(temp_phi(phidex+1)-temp_phi(phidex));
        %                 x(j)=data(1,phidex)+(data(1,phidex+1)-data(1,phidex))*dphi;
        %                 y(j)=data(2,phidex)+(data(2,phidex+1)-data(2,phidex))*dphi;
        %                 z(j)=data(3,phidex)+(data(3,phidex+1)-data(3,phidex))*dphi;
        %             end
        %             r(j)=sqrt(x(j)*x(j)+y(j)*y(j));
        %         end
        hold on
        %plot(r,z,'.','Color',cdex{i});
        plot(groupdata(1,linedex),groupdata(3,linedex),'s','Color',cdex{i});
        text(groupdata(1,linedex),groupdata(3,linedex),num2str(i));
        hold off
    end
    %axis([0 6 -3 3]);
    view(2)
elseif strcmp(plottype,'plane_xz')
    temp_phi=atan(data(2,:)./data(1,:));
    for i=1:max(data(5,:)) %Each group
        groupdata=data(:,data(5,:)==i);                     % Subsample the data so we only work with one group at a time.
        linedex=find(squeeze(groupdata(4,:))==0);   % Get the ending points for each line
        linedex=[1 linedex];                                % Index for start of each line in group
        hold on
        plot3(groupdata(1,linedex),0.0.*groupdata(1,linedex),groupdata(3,linedex),'s','Color',cdex{i});
        hold off
    end
    %axis([0 6 -3 3]);
    view(2)
elseif strcmp(plottype,'tube')    %Render as tubes
    for i=1:max(data(5,:))
        index=(data(5,:)==i);
        startpoints=find(index==1,1)-1;
        endpoints=find(data(5,:)==i & data(4,:)==0);
        startpoints=[startpoints endpoints(1:max(size(endpoints))-1)]+1;
        for j=1:max(size(startpoints))
            [xv yv zv] = tubecoords(data(1:3,startpoints(j):endpoints(j))',...
                tubewidth,nedges);
            if isempty(pcolor)
                h=surface(xv,yv,zv,'LineStyle','none','Facecolor',colors{i});  
            else
                h=surface(xv,yv,zv,'LineStyle','none','Facecolor',pcolor);
            end
        end
    end
    axis equal
else
    colors = prism(max(data(5,:)));
    for i=1:max(data(5,:))
        index=(data(5,:)==i);
        startpoints=find(index==1,1)-1;
        endpoints=find(data(5,:)==i & data(4,:)==0);
        startpoints=[startpoints endpoints(1:max(size(endpoints))-1)]+1;
        %endpoints=find(data(4,index)==0);
        %plot First segment
        for j=1:max(size(startpoints))
            plot3(data(1,startpoints(j):endpoints(j)),...
                data(2,startpoints(j):endpoints(j)),...
                data(3,startpoints(j):endpoints(j)),...
                'Color',colors(i,:),...
                'LineWidth',2);
        end
        %This section shows the direction of the coils with arrows
        %x=data(1,startpoints);
        %y=data(2,startpoints);
        %z=data(3,startpoints);
        %u=data(1,startpoints+1)-data(1,startpoints);
        %v=data(2,startpoints+1)-data(2,startpoints);
        %w=data(3,startpoints+1)-data(3,startpoints);
        %quiver3(x,y,z,u,v,w,'Color','black');
    end
    axis equal
end
hold off
title([machine_string 'Coils']);
xlabel('X');
ylabel('Y');
zlabel('Z');
coildata=coildata_save;
end

% The following two functions are shamelessly lifted from streamtube
function [x,y,z]=tubecoords(verts,width,n)
    
    d1 = diff(verts);
    zindex = find(~any(d1,2));
    verts(zindex,:) = [];
    
    if size(verts,1)<2
        x = []; y = []; z = [];
        return;
    end
    
    d1 = diff(verts);
    
    numverts = size(verts,1);
    unitnormals = zeros(numverts,3);
    
    % Radius of the tube.
    if length(width)==1
        width = repmat(width, [numverts,1]);
    else
        width(zindex) = [];
    end
    R = width/2;
    
    d1(end+1,:) = d1(end,:);
    
    x10 = verts(:,1)';
    x20 = verts(:,2)';
    x30 = verts(:,3)';
    x11 = d1(:,1)';
    x21 = d1(:,2)';
    x31 = d1(:,3)';
    
    a = verts(2,:) - verts(1,:);
    b = [0 0 1];
    c = crossSimple(a,b);
    if ~any(c)
        b = [1 0 0];
        c = crossSimple(a,b);
    end
    b = crossSimple(c,a);
    normb = norm(b); if normb~=0, b = b/norm(b); end
    %b = b*R(1);
    
    unitnormals(1,:) = b;
    
    for j = 1:numverts-1;
        
        a = verts(j+1,:)-verts(j,:);
        c = crossSimple(a,b);
        b = crossSimple(c,a);
        normb = norm(b); if normb~=0, b = b/norm(b); end
        %  b = b*R(j);
        unitnormals(j+1,:) = b;
        
    end
    
    unitnormal1 = unitnormals(:,1)';
    unitnormal2 = unitnormals(:,2)';
    unitnormal3 = unitnormals(:,3)';
    
    speed = sqrt(x11.^2 + x21.^2 + x31.^2);
    
    % And the binormal vector ( B = T x N )
    binormal1 = (x21.*unitnormal3 - x31.*unitnormal2) ./ speed;
    binormal2 = (x31.*unitnormal1 - x11.*unitnormal3) ./ speed;
    binormal3 = (x11.*unitnormal2 - x21.*unitnormal1) ./ speed;
    
    % s is the coordinate along the circular cross-sections of the tube:
    s = (0:n)';
    s = (2*pi/n)*s;
    
    % Finally, the parametric surface.
    % Each of x1, x2, x3 is an (m+1)x(n+1) matrix.
    % Rows represent coordinates along the tube.  Columns represent coordinates
    % sgcfin each (circular) cross-section of the tube.
    
    xa1 = ones(n+1,1)*x10;
    xb1 = (cos(s)*unitnormal1 + sin(s)*binormal1);
    xa2 = ones(n+1,1)*x20;
    xb2 = (cos(s)*unitnormal2 + sin(s)*binormal2);
    xa3 = ones(n+1,1)*x30;
    xb3 = (cos(s)*unitnormal3 + sin(s)*binormal3);
    
    R = repmat(R(:)',[n+1,1]);
    x1 = xa1 + R.*xb1;
    x2 = xa2 + R.*xb2;
    x3 = xa3 + R.*xb3;
    %x1 = xa1 + xb1;
    %x2 = xa2 + xb2;
    %x3 = xa3 + xb3;
    
    x = x1';
    y = x2';
    z = x3';
    
    %nx = unitnormal1;
    %ny = unitnormal2;
    %nz = unitnormal3;
end

% simple cross product
function c=crossSimple(a,b)
    c(1) = b(3)*a(2) - b(2)*a(3);
    c(2) = b(1)*a(3) - b(3)*a(1);
    c(3) = b(2)*a(1) - b(1)*a(2);
end