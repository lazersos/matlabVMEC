function error = boozer_qerror(booz_data,m,n)
%BOOZER_QERROR Returns the quasi-symmetry error
%   The BOOZER_QERROR function returns the quasi-symmetry error given a
%   boozer data structure as returned by read_boozer. The function takes
%   three arguments: the boozer data strucutre, poloidal (m) quasi-symmetry
%   target, and the toroidal (n) quasi-symmetry target. It returns an array
%   over radial gridpoints of the quasi-symmetry error defined as
%   sqrt(sum(B_mn^2))/B_00, where the appropriate modes are excluded.
%
%   Example
%       booz_data  = read_boozer('boozmn_test.nc');
%       qas_error  = boozer_qerror(booz_data,0,1);
%       qps_error  = boozer_qerror(booz_data,1,0);
%       qhs_error  = boozer_qerror(booz_data,2,1);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       2.00

% Defaults
error = zeros(1,booz_data.ns);

% Create mask of modes
if m==0
    maskdex = booz_data.xn == 0;
elseif n==0
    maskdex = booz_data.xm == 0;
else
    d1 = booz_data.xm==m;
    d2 = booz_data.xn==n;
    maskdex = and(d1,d2);
end

% Mask QXS modes
bmnc  = booz_data.bmnc;
b00 = bmnc(1,:);
bmnc(maskdex,:)=0;

% Calc Error
for i=1:booz_data.ns
    if b00(i) == 0, continue; end
    error(i) = sqrt(sum(bmnc(:,i).^2))./b00(i);
end

return;

end

