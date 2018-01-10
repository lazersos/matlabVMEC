function coil2diagno(filename)
%coil2daigno(filename) Converts makegrid coils file to DIAGNO Form
%   This is a script for creating the coils file used by DIAGNO from one
%   created for makegrid (xgrid).  Specifically the header of the DIAGNO
%   file is modified to indicate the number of elements and number of
%   current groups in the file.  Also the current group names are not
%   included in the DIAGNO coils file.  The output filename is
%   coils_diagno.suffix where suffix is the suffix of the makegrid coils
%   file.
%   Example:
%       coil2diagno('coils.machine');
%       % This produces a file 'coils_diagno.machine'
%
%   See also read_coils.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/04/11


%%%%Check some Stuff%%%%%%%%%%
if (nargin < 1)
    disp('read_coils requires a filename.');
    return
end

%First Read the Coil file
disp(' - Reading the CoilFile');
dex=strfind(filename,'.')+1;
outfile=filename(dex:max(size(filename)));
outfile=['coils_diagno.' outfile];
coil_data=read_coils(filename);
%Now sort the data
coil_data.vert=sortrows(coil_data.vert',5)';
%Now write the data
disp(' - Writting DIAGNO Version');
% Open Output file
fid=fopen(strcat(outfile),'w');
% Write the Number of Current Groups in the File
ncurs=max(coil_data.vert(5,:));
fprintf(fid,' %d\n',ncurs);
% Now write the coils
for i=1:ncurs
    disp([' - Working on Current Group:' num2str(i)]);
    sp=find(coil_data.vert(5,:)==i,1,'first');
    ep=find(coil_data.vert(5,:)==i,1,'last');
    % Now write the number of elements
    fprintf(fid,' %d\n',ep-sp+1);
    for j=sp:ep
        fprintf(fid,'%16.10e %16.10e %16.10e %16.10e\n',...
            coil_data.vert(1,j),coil_data.vert(2,j),...
            coil_data.vert(3,j),coil_data.vert(4,j));
    end
    plot3(coil_data.vert(1,sp:ep),...
        coil_data.vert(2,sp:ep),...
        coil_data.vert(3,sp:ep),'.')
    pause(1.0);
end
% Close Output File
fclose(fid);
end