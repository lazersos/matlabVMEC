function phandle=plot_vessel(data,varargin)
%PLOT_VESSEL(data,[options]) Plots the vessel data structure
%   phandles=PLOT_VESSEL plots the vessel data structure as read by the
%   READ_VESSEL function.  This function returns a handle to the object it
%   plotted.  Will return -1 if there is an error.
%
%   Options:
%
%       PLOT_VESSEL(ves_data,'phi',0.314)  Plots the cutplane corresponding
%       to phi=0.314 in radians.  If this plane does not exist then linear
%       interpolation is preformed between the nearest cutplanes.
%
%       PLOT_VESSEL(ves_data,'phi3d',0.314)  Same plot as phi but places
%       curve in the X-Z plane so it can be plotted against 3D data.
%
%       PLOT_VESSEL(ves_data,'nfp',10)  Plots the full torus where nfp
%       correspondes to the number of field periods.  Also preforms a
%       reduction of the number of faces using the REDUCEPATCH function.
%
%       PLOT_VESSEL(ves_data,'wire')  Plots a 3D wiremesh of the surface.
%       (default)
%
%       PLOT_VESSEL(ves_data,'solid')  Plots a patch face representation of
%       the surface.  Default color is grey.
%
%   Example:
%       ves_data=read_vessel('vessel.dat');
%       hves=plot_vessel(ves_data,'nfp',10,'solid');
%       set(hves,'FaceColor','red');    % Make surface red
%
%   See also read_vessel.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/11/11


%
disp('Please use plot_divertor');
return;

% Defaults
plottype='wire';
nfp=1;
hpatch=-1;

% Check Arguments
numdefargs=1;
if (nargin<1)
    disp('ERROR: plot_vessel require vessel data!');
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
            case 'phi3d'
                plottype='phi3d';
                i=i+1;
                phi=varargin{i};
            case {'wire','simple','solid','trimesh','trimesh_small'}
                plottype=varargin{i};
        end
        i=i+1;
    end
end

% Extract some data
temp=size(data.coords);
npt=temp(1);
np=max(find(data.coords(:,4)==1));
ncuts=max(data.coords(:,4));  

% Make full torus
if ~(nfp==1)
    dphi=2*pi/nfp;
    data.coords;
    tempcoords=data.coords;
    for i=2:nfp
        tempcoords(:,3)=tempcoords(:,3)+dphi;
        tempcoords(:,4)=tempcoords(:,4)+ncuts;
        data.coords=[data.coords; tempcoords];
    end
    % Now re-extract some data
    temp=size(data.coords);
    npt=temp(1);
    np=max(find(data.coords(:,4)==1));
    ncuts=max(data.coords(:,4));
end

% Convert to cartesian space
data.coords(data.coords(:,3)==2*pi,3) = 0;
x=data.coords(:,1).*cos(data.coords(:,3));
y=data.coords(:,1).*sin(data.coords(:,3));
z=data.coords(:,2);

% Control plotting
switch plottype
    case 'simple'
        % Simple plot
        hpatch=plot3(x,y,z);
    case 'phi'
        surf1.vertices = data.coords;
        surf1.faces = data.faces;
        rmax = max(surf1.vertices(1,:));
        rmin = min(surf1.vertices(1,:));
        zmax = max(surf1.vertices(3,:));
        zmin = min(surf1.vertices(3,:));
        surf2.vertices=[rmin.*cos(phi) rmin.*sin(phi) zmax; rmax.*cos(phi) rmax.*sin(phi) zmax;...
            rmin.*cos(phi) rmin.*sin(phi) zmin; rmax.*cos(phi) rmax.*sin(phi) zmin];
        surf2.faces = [1 2 3; 3 2 4];
    case 'phi_old'
        % Handle phi > maxphi
        maxphi=max(data.coords(:,3));
        if phi>maxphi
            phi=mod(phi,maxphi);
        end
        % Now find groups
        dex=find(data.coords(:,3)==phi);
        if dex
            hpatch=plot(data.coords(dex,1),data.coords(dex,2));
            xlim([min(data.coords(dex,1)) max(data.coords(dex,1))]);
            ylim([min(data.coords(dex,2)) max(data.coords(dex,2))]);
            axis image
        else
            dex2=find(data.coords(:,3)<phi,1,'last');
            dex1=dex2-np+1;
            dex3=find(data.coords(:,3)>phi,1,'first');
            dex4=dex3+np-1;
            dr=data.coords(dex3:dex4,1)-data.coords(dex1:dex2,1);
            dz=data.coords(dex3:dex4,2)-data.coords(dex1:dex2,2);
            dphi=data.coords(dex3,3)-data.coords(dex1,3);
            r=dr.*(mod(phi,dphi)/dphi)+data.coords(dex1:dex2,1);
            z=dz.*(mod(phi,dphi)/dphi)+data.coords(dex1:dex2,2);
            hpatch=plot(r,z);
            xlim([min(r) max(r)]);
            ylim([min(z) max(z)]);
            axis image
        end
    case 'phi3d'
        % Handle phi > maxphi
        maxphi=max(data.coords(:,3));
        if phi>maxphi
            phi=mod(phi,maxphi);
        end
        % Now find groups
        dex=find(data.coords(:,3)==phi);
        if dex
            x=data.coords(dex,1).*cos(phi);
            y=data.coords(dex,1).*sin(phi);
            z=data.coords(dex,2);
            hpatch=plot3(x,y,z);
            %hpatch=plot3(data.coords(dex,1),0.0.*data.coords(dex,1),data.coords(dex,2));
            xlim([min(data.coords(dex,1)) max(data.coords(dex,1))]);
            ylim([min(data.coords(dex,2)) max(data.coords(dex,2))]);
            axis image
        else
            dex2=find(data.coords(:,3)<phi,1,'last');
            dex1=dex2-np+1;
            dex3=find(data.coords(:,3)>phi,1,'first');
            dex4=dex3+np-1;
            dr=data.coords(dex3:dex4,1)-data.coords(dex1:dex2,1);
            dz=data.coords(dex3:dex4,2)-data.coords(dex1:dex2,2);
            dphi=data.coords(dex3,3)-data.coords(dex1,3);
            r=dr.*(mod(phi,dphi)/dphi)+data.coords(dex1:dex2,1);
            z=dz.*(mod(phi,dphi)/dphi)+data.coords(dex1:dex2,2);
            x=r.*cos(phi);
            y=r.*sin(phi);
            hpatch=plot3(x,y,z);
            %hpatch=plot3(r,0.0.*r,z);
            xlim([min(r) max(r)]);
            ylim([min(z) max(z)]);
            axis image
        end
    case {'wire','solid'}
        % Patch Plot over data
        % Now we calculate the vertex and face data for the patch surface
        vertex=[x y z];
        faces=[1 2 1+np+1 1+np];
        for i=2:npt-np-1
            % Neighbors
            if mod(i,np)
                faces=[faces; i i+1 i+np+1 i+np];
            end
        end
        % Now reduce the number of patch faces by nfp
        if nfp>1
            [faces,vertex]=reducepatch(faces,vertex,1/nfp);
        end
        % Now make the patch
        hpatch=patch('Vertices',vertex,'Faces',faces);
        if strcmp(plottype,'wire')
            set(hpatch,'EdgeColor','black','FaceColor','none');
        elseif strcmp(plottype,'solid')
            set(hpatch,'EdgeColor','none','FaceColor',[0.314 0.314 0.314]);
        end
    case {'trimesh'}
        vertex=[x y z];
        faces=[];
        faces = [faces; 1 1+np+1 1+1];
        faces = [faces; 1 1+np 1+np+1];
        for i=2:npt-np-np
            if mod(i,np)
                faces = [faces; i i+np+1 i+1];
                faces = [faces; i i+np   i+np+1];
            else
                faces = [faces; i i+1 i+1-np];
                faces = [faces; i i+np   i+1];
            end
        end
        offset = i;
        for i=npt-np-np+1:npt-np
            if mod(i,np)
                faces = [faces; i i-offset+1 i+1];
                faces = [faces; i i-offset   i-offset+1 ];
            else
                faces = [faces; i 1           i-np+1];
                faces = [faces; i i-offset+np 1 ];
            end
        end
        hpatch=patch('Vertices',vertex,'Faces',faces);
        set(hpatch,'EdgeColor','black','FaceColor','none');
        fid=fopen('vessel_trimesh.dat','w');
        fprintf(fid,['MACHINE:  ' data.machine '\n']);
        fprintf(fid,['DATE:' datestr(now,'mm-dd-yy') '\n']);
        fprintf(fid,'%6.6d %6.6d\n',size(vertex,1),size(faces,1));
        fprintf(fid,'%20.10E %20.10E %20.10E\n',vertex');
        fprintf(fid,'%5.5d %5.5d %5.5d\n',faces');
        fclose(fid);
    case {'trimesh_small'}
        vertex=[x y z];
        faces=[];
        faces = [faces; 1 1+np+1 1+1];
        faces = [faces; 1 1+np 1+np+1];
        for i=2:npt-np-np
            if mod(i,np)
                faces = [faces; i i+np+1 i+1];
                faces = [faces; i i+np   i+np+1];
            else
                faces = [faces; i i+1 i+1-np];
                faces = [faces; i i+np   i+1];
            end
        end
        offset = i;
        for i=npt-np-np+1:npt-np
            if mod(i,np)
                faces = [faces; i i-offset+1 i+1];
                faces = [faces; i i-offset   i-offset+1 ];
            else
                faces = [faces; i 1           i-np+1];
                faces = [faces; i i-offset+np 1 ];
            end
        end
        [faces2, vertex2]=reducepatch(faces,vertex,0.5/nfp);
        hpatch=patch('Vertices',vertex2,'Faces',faces2);
        set(hpatch,'EdgeColor','black','FaceColor','none');
        fid=fopen('vessel_trimesh.dat','w');
        fprintf(fid,['MACHINE:  ' data.machine '\n']);
        fprintf(fid,['DATE:' datestr(now,'mm-dd-yy') '\n']);
        fprintf(fid,'%6.6d %6.6d\n',size(vertex2,1),size(faces2,1));
        fprintf(fid,'%20.10E %20.10E %20.10E\n',vertex2');
        fprintf(fid,'%5.5d %5.5d %5.5d\n',faces2');
        fclose(fid);
end
phandle=hpatch;

end

