function [r_out, v,z_out] = beams3d_suv2rzp(s, u, v, S_ARR,U_ARR, raxis,phiaxis,zaxis)
% Input Parameters:
% s: Normalized Toroidal Flux
% u: Poloidal Angle
% v: Toroidal Angle
% r_out: Major Radius Distance R
% z_out: Vertical Distance Z
% phi_out: Toroidal Angle (Output)


[s,u,v]=ndgrid(s,u,v);
[R_ARR,PHI_ARR,Z_ARR]=ndgrid(raxis,phiaxis,zaxis);
hr= raxis(2)-raxis(1);
hphi= phiaxis(2)-phiaxis(1);
hz= zaxis(2)-zaxis(1);

% Constants for the splines
nphi = numel(phiaxis);
pi2 = 2 * pi;
fnorm_max = 1e5;
max_iterations = 1000;
tolerance = 1e-9;

% Begin Newton Method
residual = 1.0;
factor = 1.0;

data.nphi=nphi;
data.S_ARR=S_ARR;
data.raxis=raxis;
data.zaxis=zaxis;
[r_out,z_out] = beams3d_magaxis(data);

%  r_out = ones(size(phiaxis')).*(raxis(1)+(raxis(numel(raxis))-raxis(1))*.75);
%  z_out=zeros(size(phiaxis'));

% PHI does not change
phi_out = mod(v, max(phiaxis));

r_out=interp1(phiaxis,r_out,squeeze(phi_out(1,1,:)));
z_out=interp1(phiaxis,z_out,squeeze(phi_out(1,1,:)));
r_out=permute(repmat(r_out,1,size(phi_out,1),size(phi_out,2)),[2 3 1]);
z_out=permute(repmat(z_out,1,size(phi_out,1),size(phi_out,2)),[2 3 1]);

% Adjust u
u = mod(u, pi2);

x0 = s .* cos(u);
y0 = s .* sin(u);

%fnorm = x0.^2 + y0.^2;
%fnorm = min(1 ./ fnorm, fnorm_max);
n = 1;
X_ARR=S_ARR.*cos(U_ARR);
Y_ARR=S_ARR.*sin(U_ARR);
X=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,X_ARR,'linear','nearest');
Y=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,Y_ARR,'linear','nearest');

[~, dXdR_ARR, dXdZ_ARR] = gradient(X_ARR,hr,hphi,hz);
[~,dYdR_ARR, dYdZ_ARR] = gradient(Y_ARR,hr,hphi,hz);

dXdR=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dXdR_ARR,'linear','nearest');
dXdZ=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dXdZ_ARR,'linear','nearest');
dYdR=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dYdR_ARR,'linear','nearest');
dYdZ=griddedInterpolant(R_ARR,PHI_ARR,Z_ARR,dYdZ_ARR,'linear','nearest');
%figure
%hold on
% Loop (Newton's Method)
while ( max(residual,[],'all') > tolerance && n < max_iterations) %residual > tolerance &&

    x_term = x0 - X(r_out,phi_out,z_out);
    y_term = y0 - Y(r_out,phi_out,z_out);

    detJ = dXdR(r_out,phi_out,z_out) .* dYdZ(r_out,phi_out,z_out) - dYdR(r_out,phi_out,z_out) .* dXdZ(r_out,phi_out,z_out);
    detJ = max(abs(detJ), 0.0001); % Upper bound for step size as detJ enters in denominator
    delR = (x_term .* dYdZ(r_out,phi_out,z_out) - y_term .* dXdZ(r_out,phi_out,z_out)) ./ detJ;
    delZ = (-x_term .* dYdR(r_out,phi_out,z_out) + y_term .* dXdR(r_out,phi_out,z_out)) ./ detJ;

    delR = min(max(delR, -hr(1)), hr(1));
    delZ = min(max(delZ, -hz(1)), hz(1));

    residual = (x_term.^2 + y_term.^2);%.* fnorm;

    dex=residual < 0.01;
    delR(dex) = delR(dex) * 0.5;
    delZ(dex) = delZ(dex) * 0.5;


    r_out = max(min(r_out + delR * factor, raxis(end)), raxis(1));
    z_out = max(min(z_out + delZ * factor, zaxis(end)), zaxis(1));

    n = n + 1;
    
    %scatter3(v(1:5:end),r_out(1:5:end),z_out(1:5:end),10.0,residual(1:5:end),'o');

end
%xlim([0 0.1])

% If the maximum number of iterations is reached, print a warning.
if n >= max_iterations
    disp(['Residual max; ',num2str(max(residual,[],'all'))])
    figure
    scatter3(v(:),r_out(:),z_out(:),10.0,residual(:),'o');
    hold on
    [x_b3d,y_b3d,z_b3d] = meshgrid(phiaxis,raxis,zaxis);
    %tmp=permute(S_ARR,[3 1 2]);
    contourslice(x_b3d,y_b3d,z_b3d,S_ARR,[v(1,1,:)],[],[],linspace(0,1,10),'linear');%,linspace(-1,1,20)
    contourslice(x_b3d,y_b3d,z_b3d,U_ARR,[v(1,1,:)],[],[],linspace(0,2*pi,10),'linear');
    caxis([0 1])
    xlim([0 1.1])
    warning('Maximum number of iterations reached.');
end
end
