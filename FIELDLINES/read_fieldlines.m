function data = read_fieldlines(filename,varargin)
%READ_FIELDLINES Reads the HDF5 file created by FIELDLINES
% This funciton reads the fieldlines file and returns the data from the 
% file in a structure.
%
% Example usage
%      data=read_fieldlines('fieldline_test.h5');  % Reads FIELDLINE HDF5 file
%      data=read_fieldlines('fieldline_test.h5','iota');  % Calculate Iota
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.3

liota=0;
i=1;
if nargin>1
    while i<=1
        switch varargin{i}
            case 'iota'
                liota=1;
        end
        i=i+1;
    end
end

data = read_hdf5(filename);
if ~isstruct(data)
    return;
end
data.datatype='FIELDLINES';
data.X_lines=data.R_lines.*cos(data.PHI_lines);
data.Y_lines=data.R_lines.*sin(data.PHI_lines);
data.phiend = data.PHI_lines(:,data.nsteps);
data.dphi = double(data.phiend)/double(data.nsteps-1);
% Do this for VMECplot
data.ns = data.nlines;
data.mpol = 0;
data.nu = 2;
data.ntor = 0;
data.nzeta_grid = data.npoinc;
if isfield(data,'phiaxis')
    data.nfp= round(2*pi/max(data.phiaxis));
elseif isfield(data,'phi_grid')
    data.nfp= round(2*pi/max(data.phi_grid));
end
data.input_extension=filename;
% If we have B_PHI make B_R and B_Z true B_R and B_Z
if isfield(data,'B_PHI')
    r_temp   = repmat(repmat(data.raxis,[1 data.nphi]),[1 1 data.nz]);
    data.B_R = data.B_R.*data.B_PHI./r_temp;
    data.B_Z = data.B_Z.*data.B_PHI./r_temp;
end

% Handle torlines output
if ~isfield(data,'phiaxis')
    data.nv = double(data.nv);
    data.nu = double(data.nu);
    data.nfp = double(data.nfp);
    data.phiaxis=0:1./(data.nv-1):1;
    data.phiaxis=2*pi.*data.phiaxis./data.nfp;
end

if ~liota, return; end

% Fix values that don't start at phi=0;
phi0 = data.PHI_lines(1,1);
ishift = find(data.PHI_lines(:,1) ~=0);
PHI = data.PHI_lines;
R   = data.R_lines;
Z   = data.Z_lines;
for i = ishift'
    zeta = mod(PHI(i,:),data.phiaxis(end)');
    idex = find(zeta==phi0,1,'first');
    ldex = find(data.R_lines(i,:) == 0,1,'first')-1;
    if isempty(ldex), ldex=double(data.nsteps)-1; end
    dex = idex:ldex;
    R(i,1:length(dex)) = R(i,dex);
    R(i,length(dex)+1:end) = 0.0;
    Z(i,1:length(dex)) = Z(i,dex);
    Z(i,length(dex)+1:end) = 0.0;
    PHI(i,1:length(dex)) = PHI(i,dex);
    PHI(i,length(dex)+1:end) = -1.0;
end

% Calculate rotational transform
R0 = repmat(R(1,:),[data.nlines 1]);
Z0 = repmat(Z(1,:),[data.nlines 1]);
x  = R - R0;
y  = Z - Z0;
theta = atan2(y,x);
dtheta = diff(theta')';
dtheta(dtheta<-pi) = dtheta(dtheta<-pi)+2*pi;
dtheta(dtheta>pi) = dtheta(dtheta>pi)-2*pi;
theta = cumsum(dtheta')';
theta=[zeros(data.nlines,1) theta];
rho=sqrt(x.*x+y.*y);
clear Z x y R0 Z0 dtheta zeta;
data.iota = zeros(1,data.nlines);
data.iota_err = zeros(1,data.nlines);
data.rho = zeros(1,data.nlines);
data.rho_err = zeros(1,data.nlines);
try
    for j=2:data.nlines
        dex = R(j,:) ~=0;
        f0 = fit(PHI(j,dex)',theta(j,dex)','poly1');
        ci = confint(f0)-f0.p1;
        ci = ci(:,1);
        data.iota(j) = f0.p1;
        data.iota_err(j) = max(abs(ci));
        data.rho(j) = mean(rho(j,dex),2);
        data.rho_err(j) = std(rho(j,dex),0,2);
    end
catch
    for j=2:data.nlines
        dex = R(j,:) ~=0;
        f0 = polyfit(PHI(j,dex)',theta(j,dex)',1);
        data.iota(j) = f0(1);
        data.iota_err(j) = 0.0;
        data.rho(j) = mean(rho(j,dex),2);
        data.rho_err(j) = std(rho(j,dex),0,2);
    end
end
%ldex = R >0;
%data.rho=mean(rho,2);
%if isfield(data,'iota0')
%    data.iota(1) = data.iota0;
%    data.iota_err(1) = 1.0E-6;
%else
    data.iota(1)=2.*data.iota(2)-data.iota(3);
    data.iota_err(1) = 2.*data.iota_err(2)-data.iota_err(3);
%end
return;
end

