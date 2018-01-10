function data = read_diagno_segrog(varargin)
%READ_DIAGNO_SEGROG('filename'[,flux_data]) Reads the 'diagno_segrog' file.
%   Read the diagno_seg file and returns the values as a structure.  It
%   can also append another segrog file to the strcuture.  The data 
%   structure has the following format:
%   data:
%       data:       The segrog data for each loop (loop,numfiles).
%       filename:   Cell Array of filenames.
%       datatype:   'diagno_segrog'
%
%   Options:
%       flux_data=READ_DIAGNO_SEGROG('diagno_seg.test')  Reads the
%       segrog values in the file diagno_seg.test.
%
%       flux_data=READ_DIAGNO_SEGROG('diagno_seg.test',data)  This will 
%       read a diagno_seg.test file and append the conents of the file
%       to the existing structure.
%   
%   Example:
%       segrog_data=read_diagno_segrog('diagno_seg.val000');
%       segrog_data=read_diagno_segrog('diagno_seg.val020',flux_data);
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           7/18/13

data.datatype='diagno_segrog';
filename=varargin{1};
% Open file
fid=fopen(filename);
nels=fscanf(fid,'%d',1);
segrog=fscanf(fid,'%e',[nels 1]);
fclose(fid);
if nargin==2
    data.data=[varargin{2}.data segrog];
    data.filename=[varargin{2}.filename; filename];
else
    data.data=segrog;
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