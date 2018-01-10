function write_mgrid_coil(data,filename)
%WRITE_MGRID_COIL(data,filename) Writes out an makegrid coils file.
%   WRITE_MGRID_COIL is designed to output a coils file from a coils data
%   structure (as read by read_coils) after modification.
%   Example:
%       write_mgrid_coil(data,filename);
%
%   See also read_coils.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/14/11


%%%%Check some Stuff%%%%%%%%%%
if (nargin < 2)
    disp('ERROR: Input Arguments! write_mgrid_coil(data,filename)');
    return
end

% Open the file for output
fid=fopen(filename,'w');

% Write the Header
fprintf(fid,'periods %2d\n',data.periods);
fprintf(fid,['begin filament   ! Created by MATLAB ' date '\n']);
fprintf(fid,'mirror NUL\n');

% Now loop through coils
% for i=1:size(data.vert,2)
%     fprintf(fid,'%+20.10e  %+20.10e  %+20.10e  %+20.10e',data.vert(1,i),data.vert(2,i),data.vert(3,i),data.vert(4,i));
%     if data.vert(4,i)==0
%         fprintf(fid,'  %d  %s',data.vert(5,i),char(data.current_name{data.vert(5,i)}'));
%     end
%     fprintf(fid,'\n');
% end
% fprintf(fid,'end\n');

i=1;
ctemp=1;
enddex=size(data.vert,2);
dex=find(data.vert(4,:)==0);
while (i<enddex)
    dex2=dex(ctemp)-1;
    fprintf(fid,'%+20.10e  %+20.10e  %+20.10e  %+20.10e\n',data.vert(1:4,i:dex2));
    i=dex2+1;
    fprintf(fid,'%+20.10e  %+20.10e  %+20.10e  %+20.10e  %d  %s\n',data.vert(1,i),data.vert(2,i),data.vert(3,i),data.vert(4,i),data.vert(5,i),char(data.current_name{data.vert(5,i)}'));
    ctemp=ctemp+1;
    i=i+1;
end
fprintf(fid,'end\n');

% Close the file
fclose(fid);

end