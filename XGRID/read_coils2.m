function data=read_coils2(filename,varargin)
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
%   Version:        2.0
%   Date:           12/21/10

%%%%Check some Stuff%%%%%%%%%%
if (nargin == 0)
    disp('read_coils requires a filename.');
    data=1;
    return
end
fid=fopen(filename,'r');
while ~feof(fid)
    line=strtrim(fgetl(fid));
    first_word=strtrim(sscanf(line,'%s %*s'));
    switch lower(first_word)
        case 'periods'          % Period declaration
            data.periods=sscanf(line,'%*s %d');
        case 'begin'% Coil deffinintion
            type=strtrim(sscanf(line,'%*s %s'));
            % Increment number of coils of a given type
            if isfield(data,['n_' type])
                data.(['n_' type])=data.(['n_' type])+1;
            else
                data.(['n_' type])=1;
            end
            % Now create the strucutre
            coilname=[type num2str(data.(['n_' type]))];
            % Now handle reading the coil
            switch type
                case 'filament'
                    coildata=read_fillament(fid);
                otherwise
                    disp(['ERROR: Unknown coil type detected ' type]);
                    data=-2;
                    return
            end
            % Now add the coildata structure to the data structure
            data.(strtrim(coilname))=coildata;
    end
end
data.datatype='multi_coil_data';
fclose(fid);
return
end

function coil_data=read_fillament(fid)
x=[];
y=[];
z=[];
cur=[];
cgroup=[];
names={};
% Always define rotate and mirror even if not read
coil_data.mirror=0;
coil_data.rotate=0;
line=lower(strtrim(fgetl(fid)));
while ~strcmp(line,'end')
    if strcmp(line(1:6),'mirror')
        mirror=strtrim(sscanf(line,'%*s %s'));
        if strcmp(mirror,'mirror')
            coil_data.mirror=1;
        else
            coil_data.mirror=0;
        end
    elseif strcmp(line(1:6),'rotate')
        coil_data.rotate=1;
    elseif strcmp(line(1),'#') % do nothing
    else %Read a line
        temp=sscanf(line,'%f %f %f %f');
        group=-1;
        x=[x temp(1)];
        y=[y temp(2)];
        z=[z temp(3)];
        cur=[cur temp(4)];
        cgroup=[cgroup group];
        if temp(4)==0.0  % Check for group and name
            group=sscanf(line,'%*f %*f %*f %*f %d');
            if ~isempty(group)
                name=sscanf(line,'%*f %*f %*f %*f %*d %s');
                names=[names char(name)];
                cgroup(cgroup == -1)=group;
            end
        end
    end
    line=strtrim(fgetl(fid));
end
coil_data.x=x;
coil_data.y=y;
coil_data.z=z;
coil_data.cur=cur;
coil_data.group=cgroup;
coil_data.names=names;
return
end