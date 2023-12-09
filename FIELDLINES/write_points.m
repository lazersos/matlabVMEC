function write_points(data,n,filename)
%WRITE_POINTS Outputs fieldlines data in DESCUR format
%   This routine outputs a FIELDLINES trace in DESCUR format for fitting.
%   The code takes a fieldlines data structure (data), a field line label
%   (n), and an output filename as input.
%
%   Usage:
%       line_data=read_fieldlines('fieldlines_test.h5');
%       write_points(line_data,64,'descur_in.txt');
%
%   See also read_fieldlines.
%
%   Created by: S. Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    1.00
%   Date:       12/4/2023
%
%

phimax = max(data.phiaxis);
nfp = round(2*pi/max(data.phiaxis));
ntheta = round(data.nsteps/data.npoinc);
nzeta = data.npoinc;
fid = fopen(filename,'w');
fprintf(fid,'  %d  %d  %d\n',ntheta,nzeta,nfp);
outarray = zeros(3,ntheta.*nzeta);
k1 = 1;
k2 = ntheta;
for i=1:nzeta
    % Get axis info
    R0 = data.R_lines(1,i);
    Z0 = data.Z_lines(1,i);
    % Get indexes
    dex = i:data.npoinc:data.nsteps-1;
    dex = dex(1:ntheta);
    % Extract R and Z
    R = data.R_lines(n,dex);
    Z = data.Z_lines(n,dex);
    plot(R0,Z0,'+k'); hold on;
    plot(R,Z,'.k');
    % Calc and sort by theta
    theta = atan2(Z-Z0,R-R0);
    [theta,sdex] = sort(theta);
    R = R(sdex);
    Z = Z(sdex);
    R = [R R(1)];
    Z = [Z Z(1)];
    % Spline to dl
    dr = diff(R);
    dz = diff(Z);
    dl = sqrt(dr.*dr+dz.*dz);
    l = cumsum(dl);
    l = [0 l./l(end)];
    l_out = linspace(0,1,ntheta+1);
    l_out = l_out(1:ntheta);
    %R = smooth(R);
    %Z = smooth(Z);
    outarray(1,k1:k2) = pchip(l,R,l_out);
    outarray(3,k1:k2) = pchip(l,Z,l_out);
    %outarray(1,k1:k2) = data.R_lines(n,dex);
    outarray(2,k1:k2) = data.PHI_lines(n,dex);
    %outarray(3,k1:k2) = data.Z_lines(n,dex);
    hold on;
    plot(outarray(1,k1:k2),outarray(3,k1:k2));
    pause(1);
    k1 = k2+1;
    k2 = k1+ntheta-1;
end
outarray(2,:) = mod(outarray(2,:),phimax);
fprintf(fid,'  %20.12E  %20.12E  %20.12E\n',outarray);
fclose(fid);
end

