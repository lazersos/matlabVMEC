function sense_fluxloops(data,varargin)
%SENSE_FLUXLOOPS(filename) Calculate the spatial sensativity of a fluxloop.

%Set Some defaults
%Handle Input arguments
offset=1;
if nargin>offset
    for i=1:nargin-offset
        switch varargin{i}
            case 'color'
                i=i+1;
                linecolor=varargin{i};
            case 'linewidth'
                i=i+1;
                linewidth=varargin{i};
            case 'loopname'
                loopnames=1;
            case 'scale'
                i=i+1;
                scale=varargin{i};
        end
    end
end
% Grid Parameters
rmin=3.5;
rmax=5;
zmin=-1.5;
zmax=1.5;
phimin=2*pi*2/4;
phimax=2*pi*3/4;
nr=40;
nz=40;
nphi=36;
% Create Axes
raxis=rmin:(rmax-rmin)/(nr-1):rmax;
zaxis=zmin:(zmax-zmin)/(nz-1):zmax;
phiaxis=phimin:(phimax-phimin)/(nphi-1):phimax;
% Read Space Coordinates
x=zeros(nr,nphi,nz);
y=zeros(nr,nphi,nz);
z=zeros(nr,nphi,nz);
r=zeros(nr,nphi,nz);
phi=zeros(nr,nphi,nz);
for i=1:nr
    for j=1:nphi
        for k=1:nz
            x(i,j,k)=raxis(i)*cos(phiaxis(j));
            y(i,j,k)=raxis(i)*sin(phiaxis(j));
            z(i,j,k)=zaxis(k);
            r(i,j,k)=raxis(i);
            phi(i,j,k)=phiaxis(j);
        end
    end
end
%plot3(squeeze(x(:,:,1)),squeeze(y(:,:,1)),squeeze(z(:,:,1)),'.')
% Here we convert the fluxloop data into a form similar to the coils data
coil_data.verts=zeros(5,data.nloops*max(data.nels));
k=0;
for i=1:data.nloops
    for j=1:data.nels(i)
        k=k+1;
        coil_data.vert(1,k)=data.loops{i}(1,j);
        coil_data.vert(2,k)=data.loops{i}(2,j);
        coil_data.vert(3,k)=data.loops{i}(3,j);
        coil_data.vert(4,k)=1;
        coil_data.vert(5,k)=i;
    end
    coil_data.vert(4,data.nels(i))=0;
end
% External currents array
extcur=zeros(1,data.nloops);
extcur(1)=1;
extcur(2)=1;
extcur(3)=1;
% Magnetic Field
bx=zeros(nr,nphi,nz);
by=zeros(nr,nphi,nz);
bz=zeros(nr,nphi,nz);
coil_data=coil_biot_prep(coil_data);
for i=1:nr
    for j=1:nphi
        for k=1:nz
            [bx(i,j,k) by(i,j,k) bz(i,j,k)]=coil_biot(coil_data,...
                x(i,j,k),y(i,j,k),z(i,j,k),extcur);
        end
    end
end
bphi=by.*cos(phi)-bx.*sin(phi);
br=bx.*cos(phi)+by.*sin(phi);
% Plot Coil and grid
%xlim([-rmax rmax]);
%ylim([-rmax rmax]);
%zlim([-zmax zmax]);
hold on
hsurf=zeros(1,nphi);
for i=1:data.nloops
    if extcur(i)>0
%        plot3(data.loops{i}(1,:),...
%            data.loops{i}(2,:),...
%            data.loops{i}(3,:),'k')
        
        %pixplot(squeeze(r(:,26,:)),squeeze(z(:,26,:)),squeeze(bz(:,26,:)))
    end
end
for j=1:nphi
    %hsurf(j)=surf(squeeze(x(:,j,:)),...
    %    squeeze(y(:,j,:)),...
    %    squeeze(z(:,j,:)),...
    %    squeeze(bz(:,j,:)));
    %set(hsurf(j),'EdgeColor','none');
    axis image
    vecslice(squeeze(r(:,j,:)),...
        squeeze(z(:,j,:)),...
        squeeze(br(:,j,:)),...
        squeeze(bphi(:,j,:)),...
        squeeze(bz(:,j,:)));
    title(['Phi=' num2str(phiaxis(j)) ' radians'])
    colorbar
    pause(1.0);
end
view(3)
hold off
clf
colormap jet
xlim([rmin rmax]);
ylim([phimin phimax]);
zlim([zmin zmax]);
p=patch(isosurface(r,phi,z,bz,.75*max(max(max(bz)))),'EdgeColor','red','FaceColor','none');
p=patch(isosurface(r,phi,z,bz,.50*max(max(max(bz)))),'EdgeColor','green','FaceColor','none');
p=patch(isosurface(r,phi,z,bz,.25*max(max(max(bz)))),'EdgeColor','blue','FaceColor','none');
p=patch(isosurface(r,phi,z,bz,-.75*max(max(max(bz)))),'EdgeColor','red','FaceColor','none');
p=patch(isosurface(r,phi,z,bz,-.50*max(max(max(bz)))),'EdgeColor','green','FaceColor','none');
p=patch(isosurface(r,phi,z,bz,-.25*max(max(max(bz)))),'EdgeColor','blue','FaceColor','none');
return
end