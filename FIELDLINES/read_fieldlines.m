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


% Calculate rotational transform
for i=1:data.nsteps
    x(:,i) = data.R_lines(:,i)-data.R_lines(1,i);
    y(:,i) = data.Z_lines(:,i)-data.Z_lines(1,i);
end
theta = atan2(y,x);
dtheta = diff(theta')';
dtheta(dtheta<-pi) = dtheta(dtheta<-pi)+2*pi;
dtheta(dtheta>pi) = dtheta(dtheta>pi)-2*pi;
theta = cumsum(dtheta')';
theta=[zeros(data.nlines,1) theta];
rho=sqrt(x.*x+y.*y);
try
    for j=2:data.nlines
        f0 = fit(data.PHI_lines(j,:)',theta(j,:)','poly1');
        ci = confint(f0)-f0.p1;
        ci = ci(:,1);
        data.iota(j) = f0.p1;
        data.iota_err(j) = max(abs(ci));
    end
catch
    for j=2:data.nlines
        f0 = polyfit(data.PHI_lines(j,:)',theta(j,:)',1);
        data.iota(j) = f0(1);
        data.iota_err(j) = 1.0;
    end
end
data.rho=mean(rho,2);
%if isfield(data,'iota0')
%    data.iota(1) = data.iota0;
%    data.iota_err(1) = 1.0E-6;
%else
    data.iota(1)=2.*data.iota(2)-data.iota(3);
    data.iota_err(1) = 2.*data.iota_err(2)-data.iota_err(3);
%end
return;
end

