function data = read_vessel(filename)
%READ_VESSEL(filename) Reads the vessel data into a structure
%   READ_VESSEL reads a vessel data file.  The vessel data file is a file
%   with the format:
%    
%    MACHINE:  LHD
%    DATE:  08-12-11
% 
%    0.5094E+01  0.1016E+01  0.0000E+00
%    0.5065E+01  0.1006E+01  0.0000E+00
%    0.5037E+01  0.9960E+00  0.0000E+00
%
%   Here the string after MACHINE is the name of the machine.  The DATE
%   string contains information about the date of the file.  The array that
%   follows contains a series of datapoints (R, Z, phi) [m, m, rad].  These
%   datapoints define a piecewise defined limiter surface at each phi.
%
%   data:
%       datatype:   Identifies the structure type. 'vessel'
%       machine:    Machine string from file.
%       date:       Date string from file.
%       coords:     Coordinate data (R,Z,Phi,Face)
%
%   Example:
%       ves_data=read_vessel('vessel.dat');
%
%   See also plot_vessel.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/11/11


% Check arguments
if nargin<1
    disp('ERROR: read_vessel requires filename');
    return
end
% Define datatype
data.datatype='vessel';
% Open File
fid=fopen(filename,'r');
% Read Header
header_line1=fgetl(fid);
header_line2=fgetl(fid);
data.machine=strtrim(header_line1(strfind(header_line1,':')+1:numel(header_line1)));
data.date=strtrim(header_line2(strfind(header_line2,':')+1:numel(header_line2)));
% Scan through file
first=1;
group=1;
while ~feof(fid)
    line=fgetl(fid);
    if ~isempty(line)
        temp=sscanf(line,'%e %e %e',[1 3]);
        r=temp(1);
        z=temp(2);
        phi=temp(3);
        if first
            data.coords=[r z phi group];
            first=0;
            groupphi=phi;
            rfirst=r;
            zfirst=z;
            phifirst=phi;
        else
            if phi~=groupphi
                % First add first point as last
                data.coords=[data.coords; rfirst zfirst phifirst group];
                % Now set current point as first
                rfirst=r;
                zfirst=z;
                phifirst=phi;
                % Update group number and groupphi
                group=group+1;
                groupphi=phi;
            end
            data.coords=[data.coords; r z phi group];
        end
    end
end
% Add last point
data.coords=[data.coords; rfirst zfirst phifirst group];
fclose(fid);
return

end

