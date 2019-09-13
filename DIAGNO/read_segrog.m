function data=read_segrog(filename,varargin)
%READ_SEGROG(filename[,'scale',1.0]) Reads the DIAGNO flux loop file.
%   This function reads the DIAGNO Fluxloops file and returns the contents
%   as a data structure.  This structure has the following elements:
%
%   data:
%       nloops:     Number of flux loops
%       loops:      cell array of element coordinates (4,nels)
%       loopnames:  cell array of loop names
%       nels:       Array of number of elements per loop.
%       datatype:   'fluxloop'
%
%   Options:
%       READ_SEGROG(filename,'scale',k)  The coordiantes will be
%       multiplied by k to allow for conversion of units.
%   
%   Example:
%       loop_data=read_segrog('diagno_segrog.machine','scale',0.001);
%
%   See also read_fluxloops.
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           9/11/19

%Set Some defaults
%Handle Input arguments
offset=1;
scale=1;
if nargin>offset
    for i=1:nargin-offset
        switch varargin{i}
            case 'scale'
                i=i+1;
                scale=varargin{i};
        end
    end
end
data.datatype='fluxloop';
%   Plots Saddle Coils
fid=fopen(filename,'r');
% Read Header
line=fgetl(fid);
data.nloops=sscanf(line,'%d');
% Allocate cell array for each loop
data.loops=cell(1,data.nloops);
data.loopname=cell(1,data.nloops);
data.nels=zeros(1,data.nloops);
for i=1:data.nloops
    line=fgetl(fid);
    val=sscanf(line,'%d %d %d');
    nels=val(1);
    data.nels(i)=nels;
    t1=val(2);
    t2=val(3);
    data.loopname{i}=char(sscanf(line,'%*d %*d %*d %s'));
    line=fgetl(fid);
    val=sscanf(line,'%f %f %f %f');
    xl=val(1);
    yl=val(2);
    zl=val(3);
    cl=val(4);
    for j=2:nels
        line=fgetl(fid);
        val=sscanf(line,'%f %f %f %f');
        xl=[xl val(1)];
        yl=[yl val(2)];
        zl=[zl val(3)];
        cl=[cl val(4)];
    end
    data.loops{i}=[xl.*scale ; yl.*scale; zl.*scale; cl];
end
fclose(fid);
return
end