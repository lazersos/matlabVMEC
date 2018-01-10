function data=read_coils(filename,varargin)
%READ_COILS(filename) Reads the coils data from a coils file.
%   READ_COILS reads an VMEC MAKEGRID coils file into a strcuture.  The
%   structure has the following fields:
%       periods:        Number of field periods in the file.
%       current_name:   Names of current groups
%       vert:
%           vert(1,:)  X
%           vert(2,:)  Y
%           vert(3,:)  Z
%           vert(4,:)  current
%           vert(5,:)  group
%       datatype:       'coil_data' For identification by other routines.
%
%   See also read_coils, plot_coils, coil_pies.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.2
%   Date:           12/21/10

%  04/25/12 -SAL  Modified so both lower and upper case END are read
%                 properly


%%%%Check some Stuff%%%%%%%%%%
if (nargin == 0)
    disp('read_coils requires a filename.');
    data=1;
    return
end
fid=fopen(filename,'r');
% First pass get file size
nlines=1;
disp(' - Getting File Size.');
while ~feof(fid);
    line=fgetl(fid);
    nlines=nlines+1;
end
fclose(fid);
vert=zeros(5,nlines-4);
% Second pass gets the data
fid=fopen(filename,'r');
% Read Header Data
disp(' - Parseing Header');
line=fgetl(fid);
data.periods=sscanf(line,'%*s %d');
line=fgetl(fid);
line=fgetl(fid);
% First pass count number of filaments
disp(' - Reading Coils');
% Read first line
vals=fscanf(fid,'%e %e %e %e\n',4)';
vert(1,1)=vals(1);
vert(2,1)=vals(2);
vert(3,1)=vals(3);
vert(4,1)=vals(4);
cgrouplist=1;
cgroup_last=0;
i=1;
ilast=1;
efound=0;
while ~feof(fid) && ~efound
    line=fgetl(fid);
    if ~strcmp(lower(strtrim(line)),'end')
        vline=sscanf(line,'%e',4)';
        i=i+1;
        vert(1,i)=vline(1);
        vert(2,i)=vline(2);
        vert(3,i)=vline(3);
        vert(4,i)=vline(4);
        if vline(4) == 0
            cgroup=sscanf(line,'%*e %*e %*e %*e %d',1);
            if isempty(cgroup)
                cgroup=1;
            end
            vert(5,ilast:i)=cgroup;
            ilast=i+1;
            if (cgroup_last ~= cgroup)
                disp(['     Coil Group:' num2str(cgroup_last)]);
            end
            cgroup_last=cgroup;
            cname=sscanf(line,'%*e %*e %*e %*e %*d %s',1);
            if isempty(cname)
                cname='EXTCUR001';
            end
            cgrouplist=[cgrouplist; cgroup];
            data.current_name{cgroup}=char(cname);
        end
    else
        efound=1;
    end
end
fclose(fid);
data.vert=vert(:,1:max(size(vert))-1);
data.datatype='coil_data';
return;
end