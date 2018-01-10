function mgrid_poincm(mgrid_data)
%   MGRID_POINCM Does field line following

% Set some defaults
extcur=1:mgrid_data.nextcur;
x=zeros(mgrid_data.nr,mgrid_data.nz,mgrid_data.nphi);
y=zeros(mgrid_data.nr,mgrid_data.nz,mgrid_data.nphi);
z=zeros(mgrid_data.nr,mgrid_data.nz,mgrid_data.nphi);
cosphi=cos(mgrid_data.phi);
sinphi=sin(mgrid_data.phi);
for i=1:mgrid_data.nr
    for j=1:mgrid_data.nz
        for k=1:mgrid_data.nphi
            x(i,j,k)=mgrid_data.raxis(i).*cosphi(k);
            y(i,j,k)=mgrid_data.raxis(i).*sinphi(k);
            z(i,j,k)=mgrid_data.zaxis(j);
        end
    end
end
% Handle input arguments

% Recompose total bx by bz
bx=zeros(mgrid_data.nr,mgrid_data.nz,mgrid_data.nphi);
by=zeros(mgrid_data.nr,mgrid_data.nz,mgrid_data.nphi);
bz=zeros(mgrid_data.nr,mgrid_data.nz,mgrid_data.nphi);
for i=1:mgrid_data.nextcur
    bx(:,:,:)=bx+mgrid_data.bx(:,:,:,i).*extcur(i);
    by(:,:,:)=by+mgrid_data.by(:,:,:,i).*extcur(i);
    bz(:,:,:)=bz+mgrid_data.bz(:,:,:,i).*extcur(i);
end
% Create normal vectors
mbinv=1./sqrt(bx.*bx+by.*by+bz.*bz);
nx=bx.*mbinv;
ny=by.*mbinv;
nz=bz.*mbinv;
% Create Starting Points
dx=.5*(x(2,1,1)-x(1,1,1));
dy=0;
dz=.5*(z(1,2,1)-z(1,1,1));
startx=x(1:mgrid_data.nr-1,1:mgrid_data.nz-1,1)+dx;
starty=zeros(mgrid_data.nr-1,mgrid_data.nz-1)+dy;
startz=z(1:mgrid_data.nr-1,1:mgrid_data.nz-1,1)+dz;
% Create a bounding domain
bdr0=3.8;
bdr=0.6;
bdz0=0.0;
% Calculate the r location of all points
startr=sqrt((startx-bdr0).*(startx-bdr0)+startz.*startz);
index=(startr<bdr);
% Get 2D interpolant
Fbx=TriScatteredInterp(x(:),y(:),bx(:));
Fby=TriScatteredInterp(x(:),y(:),by(:));
Fbz=TriScatteredInterp(x(:),y(:),bz(:));
% Initialize lines
xline=startx;
yline=starty;
zline=startz;
bxline=Fbx(startx,starty);
byline=Fby(startx,starty);
bzline=Fbz(startx,starty);
% Now begin loop
for n=2:mgrid_data.nphi-1
    % Define next plane
    x1=x(1,1,n); y1=y(1,1,n); z1=z(1,1,n);
    x2=x(1,2,n); y2=y(1,2,n); z2=z(1,2,n);
    x3=x(2,1,n); y3=y(2,1,n); z3=y(1,2,n);
    tempx=xline(:,:,n-1); tempy=yline(:,:,n-1); tempz=zline(:,:,n-1);
    x4=tempx(index); y4=tempy(index); z4=tempz(index);
    tempx=bxline(:,:,n-1); tempy=byline(:,:,n-1); tempz=bzline(:,:,n-1);
    x5=x4+tempx(index); y5=y4+tempy(index); z5=z4+tempz(index);
    for i=1:mgrid_data.nr-1
        for j=1:mgrid_data.nz-1;
            num=det([ 1 1 1 1 ; x1 x2 x3 x4(i,j); y1 y2 y3 y4(i,j); z1 z2 z3 z4(i,j)]);
            denum=det([ 1 1 1 0 ; x1 x2 x3 x5(i,j)-x4(i,j); y1 y2 y3 y5(:,:)-y4(i,j); z1 z2 z3 z5(i,j)-z4(i,j)]);
            t=-num./denum;
        end
    end
    
end



end