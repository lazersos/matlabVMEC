function coil_data=focus_to_coil(data)
%FOCUS_TO_COIL(data) Convert Focus coil data to coil structure.
%   FOCUS_TO_COIL(data) takes a Focus data structure and returns a coils
%   data structure.
%
%   Usage:
%   focus_data=read_focus('test.fo.h5');
%   coil_data=focus_to_coil(focus_data);
%   plot_coil(coil_data)
%
%   See also read_focus.
%
%   Written by:     S. Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           02/08/17


% Defaults
ncoil = size(data.coilspace,1);

% Setup the harmonic data
mn    = data.NFcoil+1;
nharm = (mn*6)-3+1; % includes current
theta=0:2*pi/double(data.NDcoil-1):2*pi;
for i=1:mn
    cth(:,i) = cos(theta'.*double(i-1));
end
for i=2:mn
    sth(:,i-1) = sin(theta'.*double(i-1));
end

% Indexing for coil
dxc = [2 5:6:nharm];
dxs = 8:6:nharm;
dyc = [3 6:6:nharm];
dys = 9:6:nharm;
dzc = [4 7:6:nharm];
dzs = 10:6:nharm;

% Extract the coil
farr=data.coilspace(ncoil,:);
vert=[];
for j=1:data.Ncoils
    d1=(nharm)*(j-1)+1;
    d2=d1+nharm-1;
    carr=farr(d1:d2);
    current(j) = carr(1);
    % Format xc0 yc0 zc0 xc1 xs1 y
    x = cth*carr(dxc)'+sth*carr(dxs)';
    y = cth*carr(dyc)'+sth*carr(dys)';
    z = cth*carr(dzc)'+sth*carr(dzs)';
    c = 0.0.*z+current(j);
    c(end) = 0;
    n = 0.0.*z+double(j);
    vert=[vert; x y z c n];
    coil_data.current_name{j}=['COIL_' num2str(j,'%3.3i')];
end
coil_data.vert=vert';
coil_data.periods=1;
coil_data.datatype='coil_data';
return;

end