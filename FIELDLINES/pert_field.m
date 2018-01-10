function [ bnmn, bn2 ] = pert_field( field,nline, coil0, extcur0, coilp, extcurp ,varargin)
%PERT_FIELD(field,....) Calculates the perturbed spectrum.
%   PERT_FIELD(field,nline,coil0,extcur0,coilp,extcurp) Calculated the
%   perturbed magnetic field specturm given a field line trace and
%   perturbing coil structures.
%   PARAMETERS
%       field:  A fieldline data structure such as those returned by the
%               READ_FIELDLINES code.  This trace should contain a magnetic
%               axis and at least one good flux surface.
%       nline:  Index of the flux surface in the R_lines field of the field
%               structure.
%       coil0:  Coil data structure for background field subtraction (as
%               returned by the READ_COILS function).
%       extcur0:Array of currents for the coil0 background coil.
%       coilp:  Coil data structure for the perturbing field (as returned
%               by the READ_COILS function).
%       extcurp:Array of currents for the coilp perturbing coil.
%
%   Example usage
%      fieldline = read_fieldlines('fieldlines_test.h5');
%      coil_data = read_coils('coils.test');
%      EXTCUR    = zeros(1,22);
%      EXTCURP   = zeros(1,22);
%      EXTCURP(18) = 48E3;
%      data = pert_field(fieldline,1000,coil_data,EXTCUR,coil_data,EXTCURP);
%      [m,n]=size(data);
%      pixplot((-m/2:m/2-1),(-n/2:n/2-1),abs(data))
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00
%
%   See also read_fieldlines, read_coils, pixplot.

% Handle optional arguments
lplots=0;
ltest=0;
if nargin>6
    for i=1:nargin-6
        switch varargin{i}
            case 'plot'
                lplots=1;
            case 'test'
                ltest=1;
        end
    end
end

% extract fieldline data
R = field.R_lines(nline,:);
Z = field.Z_lines(nline,:);
PHI = field.PHI_lines(nline,:);
X = R.*cos(PHI);
Y = R.*sin(PHI);
R0 = field.R_lines(1,:);
Z0 = field.Z_lines(1,:);
PHI0 = field.PHI_lines(1,:);
X0 = R0.*cos(PHI0);
Y0 = R0.*sin(PHI0);

% Calculate the stright fieldline angle
r_temp = R-R0;
z_temp = Z-Z0;
th_temp = atan2(z_temp,r_temp);
dtheta = diff(th_temp);
dtheta(dtheta<-pi) = dtheta(dtheta<-pi)+2*pi;
dtheta(dtheta>pi) = dtheta(dtheta>pi)-2*pi;
th = [0 cumsum(dtheta)];
f0 = fit(PHI',th','poly1');
iota = f0.p1;
th2 = iota.*PHI;


% Calculate field along line
BX0 = zeros(1,length(R));
BY0 = zeros(1,length(R));
BZ0 = zeros(1,length(R));
B0 = zeros(1,length(R));
if ~ltest
    disp('  - Calculating Main Field');
    for i=1:length(R)
        x = R(i).*cos(PHI(i));
        y = R(i).*sin(PHI(i));
        [BX0(i) BY0(i) BZ0(i)] = coil_biot(coil0,x,y,Z(i),extcur0);
    end
    B0 = sqrt(BX0.*BX0+BY0.*BY0+BZ0.*BZ0);
end

% Calculate perturbation field along B
BXp = zeros(1,length(R));
BYp = zeros(1,length(R));
BZp = zeros(1,length(R));
Bp = zeros(1,length(R));
if ~ltest
    disp('  - Calculating Pert. Field');
    for i=1:length(R)
        x = R(i).*cos(PHI(i));
        y = R(i).*sin(PHI(i));
        [BXp(i) BYp(i) BZp(i)] = coil_biot(coilp,x,y,Z(i),extcurp);
    end
    Bp = sqrt(BXp.*BXp+BYp.*BYp+BZp.*BZp);
end

% Now calculate perturbing field
BX = BXp-BX0;
BY = BYp-BY0;
BZ = BZp-BZ0;
B  = sqrt(BX.*BX+BY.*BY+BZ.*BZ);

% Calculate the geometry
dx = diff(X);
dy = diff(Y);
dz = diff(Z);
dl = sqrt(dx.*dx+dy.*dy+dz.*dz);
dx = dx./dl;
dy = dy./dl;
dz = dz./dl;

% Handle field direction
field_sign = 1;
if (BY(1) < 0), field_sign = -1; end

% Construct grid
t_temp = -pi:2*pi/1000:pi;
%t_temp = 0:2*pi/1000:2*pi;
npoinc = field.npoinc*field.nfp;
bn2=[];
for i=1:npoinc
    %Extract a single Poincare Section
    bx_temp = BX(i:npoinc:end);
    by_temp = BY(i:npoinc:end);
    bz_temp = BZ(i:npoinc:end);
    x_temp = X(i:npoinc:end);
    y_temp = Y(i:npoinc:end);
    z_temp = Z(i:npoinc:end);
    dx_temp = dx(i:npoinc:end);
    dy_temp = dy(i:npoinc:end);
    dz_temp = dz(i:npoinc:end);
    % Sort theta in accending order (remaped to 2*pi)
    theta = mod(th2(i:npoinc:end),2*pi);
    theta(theta > pi) = theta(theta > pi) - 2*pi; %-pi to pi
    [theta, idex] = sort(theta);
    % Handle issue with i=npoinc
    if max(idex) > length(dx_temp)
        idex=idex(idex <= length(dx_temp));
        theta = theta(1:end-1);
    end
    % Resort arrays
    bx_temp = bx_temp(idex);
    by_temp = by_temp(idex);
    bz_temp = bz_temp(idex);
    x_temp  = x_temp(idex);
    y_temp  = y_temp(idex);
    z_temp  = z_temp(idex);
    dx_temp  = dx_temp(idex);
    dy_temp  = dy_temp(idex);
    dz_temp  = dz_temp(idex);
    % Calculate poloidal direction
    px      = diff(x_temp);
    py      = diff(y_temp);
    pz      = diff(z_temp);
    pl      = sqrt(px.*px+py.*py+pz.*pz);
    px      = px./pl;
    py      = py./pl;
    pz      = pz./pl;
    % Calculate normal direction
    nx      = dy_temp(1:end-1).*pz - dz_temp(1:end-1).*py;
    ny      = dz_temp(1:end-1).*px - dx_temp(1:end-1).*pz;
    nz      = dx_temp(1:end-1).*py - dy_temp(1:end-1).*px;
    nl      = field_sign.*sqrt(nx.*nx + ny.*ny + nz.*nz);
    bn_temp = (bx_temp(1:end-1).*nx + by_temp(1:end-1).*ny + bz_temp(1:end-1).*nz)./nl;
    bn_temp = smooth(theta(1:end-1),bn_temp,'rloess');
    bn2(1:length(t_temp),i) = pchip(theta(1:end-1),bn_temp,t_temp);
end

if ltest
    ntemp = size(bn2,1);
    for i=1:ntemp
        bn2(i,:) = cosd(15.*360.*(i-1)./(ntemp-1));
    end
end

% Make array periodic
bn2(:,end+1) = bn2(:,1);

% Now fft
disp('  - Fourier Transform');
mt = pow2(nextpow2(20));
nt = pow2(nextpow2(60));
scale = numel(bn2);
bnmn = fft2(bn2)./scale;
bnmn = fftshift(bnmn);
m = size(bnmn,1);
n = size(bnmn,2);
% subsample
m0 = floor(m/2)+1;
n0 = floor(n/2)+1;
bnmn = bnmn(m0-mt/2:m0+mt/2,n0-nt/2:n0+nt/2);
n = nt;
m = mt;

if lplots
    p_temp = double(0:360.0/npoinc:360);
    figure('Color','white','Position',[1 -100 2048 768]);
    subplot(1,2,1);
    set(gca,'FontSize',24);
    pixplot(p_temp,180.*(t_temp)/pi,bn2');
    set(gca,'XTick',[1.5 90.5 180.5 270.5 359.5],...
        'XTickLabel',{'0' '90' '180' '270' '360'});
    %set(gca,'YTick',[1.5 90.5 180.5 270.5 359.5],...
    %    'YTickLabel',{'0' '90' '180' '270' '360'});
    set(gca,'YTick',[-179.5 -89.5 0 89.5 179.5],...
        'YTickLabel',{'-180' '-90' '0' '90' '180'});
    hold on;
    i_temp = iota*pi*p_temp./180.0;
    plot(p_temp,180*i_temp/pi-180,'w','LineWidth',2.0)
    hold off;
    xlabel('\phi [^o]');
    ylabel('\theta [^o]');
    title('B-Normal');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'B.n');
    subplot(1,2,2);
    set(gca,'FontSize',24);
    pixplot(-n/2:n/2,-m/2:m/2,abs(bnmn)')
    set(gca,'XTick',[-30 -15 0 15 30]+0.5,...
        'XTickLabel',{'-30' '-15' '0' '15' '30'});
    set(gca,'YTick',[-15 -10 -5  0 5 10 15]+0.5,...
        'YTickLabel',{'-15' '-10' '-5' '0' '5' '10' '15'});
    xlabel('Toroidal Mode Number (n)');
    ylabel('Poloidal Mode Number (m)');
    title('B-Normal Spectrum');
    ha = colorbar;
    set(ha,'FontSize',24);
    ylabel(ha,'(B.n)_{mn}');
    hold on;
    for i=-30:5:30
        plot([i i]+0.5,ylim,'w');
    end
    for i=[-15 -10 -5 -1 0 1 5 10 15]
        plot(xlim,[i i]+0.5,'w');
    end
    hold off;
end

return;
end

