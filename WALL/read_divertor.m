function data = read_divertor( filename )
%READ_DIVERTOR Read's a divertor file
%   READ_DIVERTOR reads a divertor data file.  The divertor data file is a
%   file with the format:
%
%MACHINE:  W7-X
%DATE:  01-10-13
%
%    5.3802400000E+00     3.6732000000E-01      1.1344640138E-01  00004  00001
%    5.4016900000E+00     3.7863000000E-01      1.1344640138E-01  00004  00001
%    5.3386200000E+00     4.9053000000E-01      1.1344640138E-01  00004  00001
%
%   Here the string after MACHINE is the name of the machine.  The DATE
%   string contains information about the date of the file.  The array that
%   follows contains a series of datapoints (R, Z, phi, np, n) [m, m, rad,
%   N/A, N/A].  These points define cross sections of various divertor
%   plates.  The first index (np) indicates the number of points per
%   cross section, the secon number is a plate identifier.
%   data:
%       datatype:   Identifies the structure type. 'vessel'
%       machine:    Machine string from file.
%       date:       Date string from file.
%       coords:     Coordinate data (R,Z,Phi,np,n)
%
%   Example:
%       div_data=read_divertor('divertor.dat');
%
%   See also plot_vessel.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/10/13
% Check arguments
if nargin<1
    disp('ERROR: read_vessel requires filename');
    return
end
% Open File
fid=fopen(filename,'r');
% Read Header
header_line1=fgetl(fid);
header_line2=fgetl(fid);
data.machine=strtrim(header_line1(strfind(header_line1,':')+1:numel(header_line1)));
data.date=strtrim(header_line2(strfind(header_line2,':')+1:numel(header_line2)));
% Read dataset
data.coords=fscanf(fid,'%E %E %E %d %d',[5 inf]);
fclose(fid);
data.datatype='divertor';


end

