function [ data ] = coil_curv( data )
%COIL_CURV(coil_data) Adds coil curvature/torsion info to coil structure.
%   COIL_CURV takes a coil data structure (as read by READ_COIL) and
%   returns a coil_data structure with torsion and curvature data added.
%   The structure adds the following fields:
%       Tx/y/z:        Tangent vector.
%       Nx/y/z:        Normal vector.
%       Bx/y/z:        Bi-normal vector (B=TxN).
%       kappa:         Curvature
%       tau:           Torsion
%       tl:            Path length
%
%   See also read_coils, plot_coils.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           11/08/16


lplot=1;

coil_ends=find(data.vert(4,:)==0);
nseg =length(coil_ends);
n1=1;
Tx=[]; Ty=[]; Tz=[];
Nx=[]; Ny=[]; Nz=[];
Bx=[]; By=[]; Bz=[];
kappa = []; tau = []; tl=[];
dx = diff(data.vert(1,:));
dy = diff(data.vert(2,:));
dz = diff(data.vert(3,:));
tg = cumsum([0 sqrt(dx.*dx+dy.*dy+dz.*dz)]);
for i=1:nseg
    n2 = coil_ends(i);
    t = tg(n1:n2)-tg(n1);
    tl = [tl t];
    % first and last point are the same
    x= data.vert(1,n1:n2);
    y= data.vert(2,n1:n2);
    z= data.vert(3,n1:n2);
    x(2:end-1) = smooth(x(2:end-1),'sgolay',3);
    y(2:end-1) = smooth(y(2:end-1),'sgolay',3);
    z(2:end-1) = smooth(z(2:end-1),'sgolay',3);
    dx     = [diff(x) x(end)-x(1)];
    dy     = [diff(y) y(end)-y(1)];
    dz     = [diff(z) z(end)-z(1)];
    dl     = sqrt(dx.*dx+dy.*dy+dz.*dz);
    % Calc T
    Txl    = dx./dl;
    Tyl    = dy./dl;
    Tzl    = dz./dl;
    norm   = sqrt(Txl.*Txl+Tyl.*Tyl+Tzl.*Tzl);
    Txl = Txl./norm; Tyl=Tyl./norm; Tzl=Tzl./norm;
    % Calc N
    ddx    = [diff(Txl) Txl(end)-Txl(1)];
    ddy    = [diff(Tyl) Tyl(end)-Tyl(1)];
    ddz    = [diff(Tzl) Tzl(end)-Tzl(1)];
    Nxl    = ddx./dl;
    Nyl    = ddy./dl;
    Nzl    = ddz./dl;
    norm   = sqrt(Nxl.*Nxl+Nyl.*Nyl+Nzl.*Nzl);
    kappal = norm;
    Nxl = Nxl./norm; Nyl=Nyl./norm; Nzl=Nzl./norm;
    % Calc B, B', torsion
    Bxl = Tyl.*Nzl - Tzl.*Nyl;
    Byl = Tzl.*Nxl - Txl.*Nzl;
    Bzl = Txl.*Nyl - Tyl.*Nxl;
    Bxp    = [diff(Bxl) Bxl(end)-Bxl(1)];
    Byp    = [diff(Byl) Byl(end)-Byl(1)];
    Bzp    = [diff(Bzl) Bzl(end)-Bzl(1)];
    Bxp    = Bxp./dl;
    Byp    = Byp./dl;
    Bzp    = Bzp./dl;
    %norm   = sqrt(Bxp.*Bxp+Byp.*Byp+Bzp.*Bzp);
    %Bxp = Bxp./norm; Byp = Byp./norm; Bzp = Bzp./norm;
    taul   = -Nxl.*Bxp-Nyl.*Byp-Nzl.*Bzp;
    % Full arrays
    kappa = [kappa kappal];
    tau = [tau taul];
    Tx = [Tx Txl];
    Ty = [Ty Tyl];
    Tz = [Tz Tzl];
    Nx = [Nx Nxl];
    Ny = [Ny Nyl];
    Nz = [Nz Nzl];
    Bx = [Bx Bxl];
    By = [By Byl];
    Bz = [Bz Bzl];
    n1 = n2 + 1;
end
% Return Arrays
data.Tx = Tx;
data.Ty = Ty;
data.Tz = Tz;
data.Nx = Nx;
data.Ny = Ny;
data.Nz = Nz;
data.Bx = Bx;
data.By = By;
data.Bz = Bz;
data.kappa = kappa;
data.tau = tau;
data.tl = tl;

% Make a few plots
if (lplot)
    cmap=colormap('lines');
    subplot(2,2,1);
    n1=1;
    for i=1:nseg
        n2=coil_ends(i);
        t = ((n1:n2)-n1)./(n2-n1);
        hold on;
        plot(t,kappa(n1:n2),'Color',cmap(mod(i-1,64)+1,:));
        hold off;
    end
    title('Curvature');
    subplot(2,2,3);
    n1=1;
    for i=1:nseg
        n2=coil_ends(i);
        t = ((n1:n2)-n1)./(n2-n1);
        hold on;
        plot(t,tau(n1:n2));
        hold off;
    end
    title('Torsion');
    subplot(2,2,[2 4]);
    quiver3(data.vert(1,:),data.vert(2,:),data.vert(3,:),Tx, Ty, Tz,'b');
    hold on;
    quiver3(data.vert(1,:),data.vert(2,:),data.vert(3,:),Nx, Ny, Nz,'g');
    quiver3(data.vert(1,:),data.vert(2,:),data.vert(3,:),Bx, By, Bz,'r');
    hold off;
    axis equal;
end

return;


end

