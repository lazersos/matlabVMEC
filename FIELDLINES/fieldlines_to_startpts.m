function fieldlines_to_startpts( data, varargin)
%FIELDLINES_TO_STARTPTS Calculates startpoints from fieldlines
%   The FIELDLINES_TO_STARTPTS routine takes data produced by the
%   FIELDLINES code (read by READ_FIELDLINES) and creates a set of start
%   points based on a single fieldlines trajectory.
%   Options:
%       'ns':           Fieldlines index to use (default 128).
%       'npts':         Number of startpoints (default 64000)
%       'filename':     Filename for file.    
%       'plot':         Make a plot showing the surface.
%
%   Usage:
%       line_data=read_fieldlines('fieldlines_test.h5');
%       fieldlines_to_startpts(line_data,'ns',64,'npts',32000);
%
%   See also read_fieldlines.
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Version:    1.00
%   Date:       10/24/14

%Defaults go here
filename='PTS0.TXT';
npts = 64000;
ns=200;
lplot=1;

% Handle varargin
if nargin > 1
    for i=1:nargin-1
        switch varargin{i}
            case 'plot'
                lplot=1;
            case 'filename'
                i=i+1;
                filename=varargin{i};
            case 'ns'
                i=i+1;
                ns=varargin{i};
            case 'npts'
                i=i+1;
                npts=varargin{i};
        end
    end
end

% Extract data
phi = data.PHI_lines(ns,:);
r   = data.R_lines(ns,:);
z   = data.Z_lines(ns,:);
phistart = min(phi);
phiend = max(phi);
if phistart == -1, phistart=0; end
phinew = phistart:(phiend-phistart)/(npts-1):phiend;

% Random distribution
%plot(mod(phinew,2*pi),'o');
phinew = phinew + rand(size(phinew))-0.5;
phinew(phinew>phiend) = phinew(phinew>phiend) - phiend;
%hold on;
%plot(mod(phinew,2*pi),'+');

% Spline to new coordinate
Rnew = pchip(phi,r,phinew);
Znew = pchip(phi,z,phinew);

% Plot
if lplot
    figure('Position',[1 1 1024 768],'Color','white');
    plot_fieldlines(data,'color','k');
    hold on;
    plot(r(1:data.npoinc:end),z(1:data.npoinc:end),'.r');
    hold off;
end




% Adjust so no particle starts greater than 2*pi
phinew = mod(phinew,2*pi);

% Output file
fid = fopen(filename,'w');
fprintf(fid,'  PHI_END = %i*%-20.10E \n',length(Rnew),2*pi*1000000);
for i=1:length(Rnew)
    fprintf(fid,'  R_start(%i) = %-20.10E ',i,Rnew(i));
    fprintf(fid,'  Z_start(%i) = %-20.10E ',i,Znew(i));
    fprintf(fid,'  PHI_start(%i) = %-20.10E \n',i,phinew(i));
end
fprintf(fid,'/\n');
fclose(fid);


end

