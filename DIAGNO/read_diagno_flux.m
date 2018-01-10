function data = read_diagno_flux(varargin)
%READ_DIAGNO_FLUX('filename'[,flux_data]) Reads the 'diagno_flux' file.
%   Read the diagno_flux file and returns the values as a structure.  It
%   can also append another flux file to the strcuture.  The data structure
%   has the following format:
%   data:
%       data:       The flux data for each loop (loop,numfiles).
%       names:      Cell Array of fluxloop names.
%       filename:   Cell Array of filenames.
%       datatype:   'diagno_flux'
%
%   Options:
%       flux_data=READ_DIAGNO_FLUX('filename')  Reads the flux values in
%       the 
%
%       data=READ_DIAGNO_FLUX('diagno_flux.test',data)  This will read a
%       diagno_flux.test file and append the conents of the file to the
%       existing structure.
%   
%   Example:
%       flux_data=read_diagno_flux('diagno_flux.val000');
%       flux_data=read_diagno_flux('diagno_flux.val020',flux_data);
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/04/11

data.datatype='diagno_flux';
filename=varargin{1};
% Open file
fid=fopen(filename);
nels=fscanf(fid,'%d',1);
flux=fscanf(fid,'%e',[nels 1]);
    line=fgetl(fid);
for i=1:nels
    line=fgetl(fid);
    names{i}=strtrim(line);
end
fclose(fid);
if nargin==2
    data.data=[varargin{2}.data flux];
    data.names=names;
    data.filename=[varargin{2}.filename; filename];
else
    data.data=flux;
    data.names=names;
    data.filename{1}=filename;
end
return
end




% % Get the first line containitng data
% dataline=fgetl(fid);
% temp=sscanf(dataline,'%e');
% vals=[];
% while ~isempty(temp)
%     vals=[vals; temp];
%     dataline=fgetl(fid);
%     line_temp=strrep(upper(dataline),'E','1');
%     line_temp=strrep(line_temp,'+','1');
%     line_temp=strrep(line_temp,'-','1');
%     if any(isletter(line_temp));
%         temp=[];
%     else
%         temp=sscanf(dataline,'%e');
%     end
% end
% % Now read each flux loop name
% names=cell(1,1);
% nameline=dataline;
% names{1}=nameline;
% while ~feof(fid)
%     nameline=fgetl(fid);
%     names=[names; nameline];
% end
% fclose(fid);
% % Now parse the structure.
% n=numel(names);
% %vals=sscanf(dataline,'%f',n);
% if nargin==2
%     data.data=[varargin{2}.data vals];
%     data.names=names;
%     data.filename=[varargin{2}.filename; filename];
% else
%     data.data=vals;
%     data.names=names;
%     data.filename{1}=filename;
% end
% return
% end