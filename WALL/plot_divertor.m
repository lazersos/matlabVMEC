function phandle=plot_divertor(data,varargin)
%PLOT_DIVERTOR(data,[options]) Plots the vessel divertor structure
%   phandles=PLOT_DIVERTOR plots the divertor data structure as read by the
%   READ_DIVERTOR function.  This function returns a handle to the object
%   it plotted.  Will return -1 if there is an error.
%
%   Options:
%
%       PLOT_DIVERTOR(div_data,'phi',0.314)  Plots the cutplane
%       corresponding to phi=0.314 in radians.  If this plane does not
%       exist then linear interpolation is preformed between the nearest 
%       cutplanes.
%
%       PLOT_DIVERTOR(div_data,'nfp',10)  Plots the full torus where nfp
%       correspondes to the number of field periods.  Also preforms a
%       reduction of the number of faces using the REDUCEPATCH function.
%
%       PLOT_DIVERTOR(div_data,'wire')  Plots a 3D wiremesh of the surface.
%       (default)
%
%       PLOT_DIVERTOR(div_data,'solid')  Plots a patch face representation
%       of the divertors.  Default color is grey.
%
%   Example:
%       div_data=read_divertor('divertor.dat');
%       hves=plot_divertor(div_data,'nfp',10,'solid');
%       set(hves,'FaceColor','red');    % Make surface red
%
%   See also read_divertor.
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           1/10/13


% Defaults
plottype='simple';
nfp=1;
hpatch=-1;
varargin_temp={};

% Check Arguments
numdefargs=1;
if (nargin<1)
    disp('ERROR: plot_divertor require divertor data!');
    return
elseif (nargin>1)
    i=1;
    while (i<=nargin-numdefargs)
        switch varargin{i}
            case 'nfp'
                i=i+1;
                nfp=varargin{i};
            case 'phi'
                plottype='phi';
                i=i+1;
                phi=varargin{i};
            case 'z'
                plottype='z';
                i=i+1;
                zcut=varargin{i};
            case {'wire','simple','solid'}
                plottype=varargin{i};
            otherwise
                varargin_temp=[varargin_temp varargin{i}];
        end
        i=i+1;
    end
end


% Convert to cartesian space
if strcmp(data.datatype,'divertor')
    nplates = max(data.coords(5,:));  % Number of discrete plates
    x=data.coords(1,:).*cos(data.coords(3,:));
    y=data.coords(1,:).*sin(data.coords(3,:));
    z=data.coords(2,:);
    dex_arr=data.coords(5,:);
elseif strcmp(data.datatype,'limiter_trimesh')
    x=data.coords(1,:);
    y=data.coords(2,:);
    z=data.coords(3,:);
end

% Control plotting
switch plottype
    case 'simple'
        % Simple plot
        if isempty(varargin_temp)
            hpatach=plot3(x,y,z,'+k');
        else
            hpatach=plot3(x,y,z,varargin_temp{:});
        end
    case 'phi'
        if exist('SurfaceIntersection','file')~=2
            disp(' Please download the SurfaceIntersection Routine from:')
            disp('     https://www.mathworks.com/matlabcentral/fileexchange/48613-surface-intersection');
            return;
        end
        r = sqrt(x.*x+y.*y);
        rmax=max(r).*1.2;
        rmin=min(r)*.8;
        zmax=max(z).*1.2;
        zmin=min(z).*1.2;
        surf1.vertices = data.coords';
        surf1.faces = data.faces';
        surf2.vertices=[rmin.*cos(phi) rmin.*sin(phi) zmax; rmax.*cos(phi) rmax.*sin(phi) zmax;...
            rmin.*cos(phi) rmin.*sin(phi) zmin; rmax.*cos(phi) rmax.*sin(phi) zmin];
        surf2.faces = [1 2 3; 3 2 4];
        [intMatrix, intSurface] = SurfaceIntersection(surf1, surf2);
        x2 = intSurface.vertices(:,1);
        y2 = intSurface.vertices(:,2);
        z2 = intSurface.vertices(:,3);
        r2 = sqrt(x2.*x2+y2.*y2);
        if isempty(varargin_temp)
            plot(r2(intSurface.edges'),z2(intSurface.edges'),'k')
        else
            plot(r2(intSurface.edges'),z2(intSurface.edges'),varargin_temp{:})
        end
    case 'z'
        if exist('SurfaceIntersection','file')~=2
            disp(' Please download the SurfaceIntersection Routine from:')
            disp('     https://www.mathworks.com/matlabcentral/fileexchange/48613-surface-intersection');
            return;
        end
        xmax=max(x).*1.2;
        xmin=min(x)*.8;
        ymax=max(y).*1.2;
        ymin=min(y).*1.2;
        surf1.vertices = data.coords';
        surf1.faces = data.faces';
        surf2.vertices=[xmin ymin zcut; xmax ymin zcut;...
            xmin ymax zcut; xmax ymax zcut];
        surf2.faces = [1 2 3; 3 2 4];
        [intMatrix, intSurface] = SurfaceIntersection(surf1, surf2);
        x2 = intSurface.vertices(:,1);
        y2 = intSurface.vertices(:,2);
        z2 = intSurface.vertices(:,3);
        r2 = sqrt(x2.*x2+y2.*y2);
        if isempty(varargin_temp)
            plot(x2(intSurface.edges'),y2(intSurface.edges'),'k')
        else
            plot(x2(intSurface.edges'),y2(intSurface.edges'),varargin_temp{:})
        end
    case 'phi2'
        if strcmp(data.datatype,'limiter_trimesh')
            % Take from (http://www.mathworks.com/matlabcentral/newsreader/view_thread/29075)
            p0 = [cos(phi)  sin(phi)  0];
            p1 = [cos(phi); -sin(phi); 0.0]; % Vector points in phi direction
            a = data.coords(:,data.faces(1,:))';
            b = data.coords(:,data.faces(2,:))';
            c = data.coords(:,data.faces(3,:))';
            for i=1:size(a,1)
                a(i,:) = a(i,:) - p0;
                b(i,:) = b(i,:) - p0;
                c(i,:) = c(i,:) - p0;
            end
            da = a*p1;
            db = b*p1;
            dc = c*p1;
            k=(da>0)+(db>0)+(dc>0);
            k=find((k==1)|(k==2)|(k==3));
            edgelist=[];
            xyzp=[];
            j=0;
            tri = data.faces';
            for i=k(1):k(end)
                edgei=[];
                % did we cross edge ab?
                if (da(i)*db(i))<=0
                    j=j+1;
                    edgei=j;
                    t=abs(da(i))/(abs(da(i))+abs(db(i)));
                    xyz0=data.coords(:,tri(i,1))'.*(1-t)+data.coords(:,tri(i,2))'.*t;
                    xyzp=[xyzp;xyz0];
                end
                
                % did we cross edge ac?
                if (da(i)*dc(i))<=0
                    j=j+1;
                    edgei=[edgei,j];
                    t=abs(da(i))/(abs(da(i))+abs(dc(i)));
                    xyz0=data.coords(:,tri(i,1))'.*(1-t)+data.coords(:,tri(i,3))'.*t;
                    xyzp=[xyzp;xyz0];
                end
                
                % did we cross edge bc?
                if (db(i)*dc(i))<=0
                    j=j+1;
                    edgei=[edgei,j];
                    t=abs(db(i))/(abs(db(i))+abs(dc(i)));
                    xyz0=data.coords(:,tri(i,2))'.*(1-t)+data.coords(:,tri(i,3))'.*t;
                    xyzp=[xyzp;xyz0];
                end
                
                edgelist=[edgelist;edgei];
            end
            
            %hold on
            %plot3(data.coords(1,:),data.coords(2,:),data.coords(3,:),'.')
            %plot3(p0(1),p0(2),p0(3),'ms')
            %plot3(xyzp(:,1),xyzp(:,2),xyzp(:,3),'r.')
            %hold off
            r = sqrt(xyzp(:,1).^2+xyzp(:,2).^2);
            hold on;
            plot(r,xyzp(:,3),'r');
            hold off;
        else
            r2 = [];
            z2 = [];
            % Plot a toroidal slice
            maxphi = max(data.coords(3,:));
            minphi = min(data.coords(3,:));
            deltaphi = maxphi-minphi;
            if phi > maxphi
                phi = phi - minphi;
                phi = mod(phi,deltaphi);
                phi = phi + minphi;
            end
            for i=1:nplates
                dex1=find(dex_arr==i,1,'first');
                dex2=find(dex_arr==i,1,'last');
                phi_temp = data.coords(3,dex1:dex2);
                np = data.coords(4,dex1);
                if (phi >= min(phi_temp)) && (phi <= max(phi_temp))
                    x1=reshape(x(dex1:dex2),[np (dex2-dex1+1)/np]);
                    z1=reshape(z(dex1:dex2),[np (dex2-dex1+1)/np]);
                    phi_temp=reshape(phi_temp,[np (dex2-dex1+1)/np]);
                    r_temp=[];
                    z_temp=[];
                    for j = 1 : np
                       r_temp(j)=interp1(phi_temp(j,:),x1(j,:),phi,'pchip');
                       z_temp(j)=interp1(phi_temp(j,:),z1(j,:),phi,'pchip');
                    end
                    hold on
                    if isempty(varargin_temp)
                        plot(r_temp,z_temp,'k');
                    else
                        plot(r_temp,z_temp,varargin_temp{:});
                    end
                    r2=[r2 pchip(1:length(r_temp),r_temp,1:length(r_temp)/39.:length(r_temp))];
                    z2=[z2 pchip(1:length(r_temp),z_temp,1:length(r_temp)/39.:length(r_temp))];
                    %hold off
                end
            end
        end
        %disp(['R_START   = ' num2str(r2,' %20.10E %20.10E %20.10E %20.10E %20.10E \n')]);
        %disp(['Z_START   = ' num2str(z2,' %20.10E %20.10E %20.10E %20.10E %20.10E \n')]);
        %disp(['PHI_START = ' num2str(z2.*0.0+phi,' %20.10E %20.10E %20.10E %20.10E %20.10E \n')]);
    case {'wire','solid'}
        % Patch Plot over data
        % Now we calculate the vertex and face data for the patch surface
        if strcmp(data.datatype,'limiter_trimesh')
            faces=data.faces';
            vertex=data.coords';
        else
            faces=[];
            vertex=[x' y' z'];
            i=1;
            for i=1:nplates
                dex1=find(dex_arr==i,1,'first');
                dex2=find(dex_arr==i,1,'last');
                %vertex=[x(dex1:dex2)' y(dex1:dex2)' z(dex1:dex2)'];
                np = data.coords(4,dex1);
                j=dex1;
                %faces = [faces; j j+1 j+np+1 j+np]; %old 4 point
                faces = [faces; j j+np+1 j+1];
                faces = [faces; j j+np j+np+1];
                while (j < dex2-np-1)
                    j=j+1;
                    if mod(j,np)
                        %faces = [faces; j j+1 j+np+1 j+np]; %old 4 point
                        faces = [faces; j j+np+1 j+1];
                        faces = [faces; j j+np j+np+1];
                    end
                end
            end
        end
        hold on;
        hpatch=patch('Vertices',vertex,'Faces',faces);
        if strcmp(plottype,'wire')
            set(hpatch,'EdgeColor','black','FaceColor','none');
        elseif strcmp(plottype,'solid')
            set(hpatch,'EdgeColor','none','FaceColor',[0.314 0.314 0.314]);
        end
        hold off;
end
phandle=hpatch;

end

