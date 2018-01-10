function data = read_diagno_induction(varargin)
%READ_DIAGNO_INDUCTION('filename') Reads the 'diagno_induction' file.
%   Read the diagno_flux file and returns the values as a structure.  It
%   can also append another flux file to the strcuture.  The data structure
%   has the following format:
%   data:
%       _data:      Data for each diagnostic.
%       _names:     Cell array of names for each diagnostic.
%       _header:    Header for each diagnostic
%       datatype:   'diagno_induction'
%
%   Options:
%       flux_data=READ_DIAGNO_INDUCTION('filename')  Reads the 
%       'diagno_induction' file.
%
%   
%   Example:
%       flux_data=read_diagno_flux('diagno_induction.test');
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.1
%   Date:           2/11/11

data.datatype='diagno_induction';
filename=varargin{1};
names={};
% Open the File
fid=fopen(filename,'r');
while ~feof(fid)
    line=strtrim(fgetl(fid));
    switch line
        case {'FLUXLOOPS'}
            type=line;
            data.([type '_header'])=fgetl(fid); % Grab the Header line
            scan_line='%*8c %d %d %e';
            line=fgetl(fid);
            oldname={sscanf(line,'%8c',1)};
            data.([type '_name'])={oldname};
            data.([type '_data'])=sscanf(line,scan_line,inf)';
        case {'BPROBES'}
            type=line;
            data.([type '_header'])=fgetl(fid); % Grab the Header line
            scan_line='%d %d %e %e %e';
            line=fgetl(fid);
            data.([type '_data'])=sscanf(line,scan_line,inf)';
        otherwise
            data.([type '_data'])=[data.([type '_data']); ...
                sscanf(line,scan_line,inf)'];
            if strcmp(type,'FLUXLOOPS')
                newname={sscanf(line,'%8c',1)};
                if ~strcmp(newname,oldname)
                    oldname=newname;
                    data.([type '_name'])=[data.([type '_name']); oldname];
                end
            end
    end     
end
fclose(fid);
return
end