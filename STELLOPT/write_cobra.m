function write_cobra(vmec_data,varargin)
%WRITE_COBRA(vmec_data) Creates a COBRAVMEC input file.
%
%   WRITE_COBRA(vmec_data)  Creates a COBRAVMEC in_cobra input file.
%
%   WRITE_COBRA(vmec_data,'theta',theta)  Creates a COBRAVMEC in_cobra 
%   input file where theta is an array of ballooning stability angles in
%   theta (degrees, default=0).
%
%   WRITE_COBRA(vmec_data,'zeta',zeta)  Creates a COBRAVMEC in_cobra 
%   input file where zeta is an array of ballooning stability angles in
%   zeta (degrees, default=0).
%
%   WRITE_COBRA(vmec_data,'ns',ns)  Creates a COBRAVMEC in_cobra 
%   input file where ns is an array of surfaces on which to calculate
%   ballooning stability. (default=2:ns-1)
%
%   See also read_vmec.
%
%   Example:
%       theta=0:90:360;                        %Angles in Degrees
%       zeta=0:90:360;
%       vmec_data=read_vmec('wout_test.nc');
%       write_cobra(vmec_data,'theta',theta,'zeta',zeta,'ns',2:vmec_data.ns-1);
%
%   Written by: Samuel Lazerson (lazerson@pppl.gov)
%   Version:    1.0
%   Date:       09/28/2011

% Get defaults
filename=['in_cobra.' strtrim(vmec_data.input_extension)];
ntheta=1;
nzeta=1;
zeta=0;
theta=0;
surfs=2:vmec_data.ns-1;
nsurfs=length(surfs);

numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    i=1;
    while i<=(nargin-numdefargs)
        if strcmp(varargin{i},'theta')
            i=i+1;
            theta=varargin{i};
            ntheta=length(theta);
        elseif strcmp(varargin{i},'zeta')
            i=i+1;
            zeta=varargin{i};
            nzeta=length(zeta);
        elseif strcmp(varargin{i},'ns')
            i=i+1;
            surfs=varargin{i};
            nsurfs=length(ns);
        end
        i=i+1;
    end
end


% Write the file
fid=fopen(filename,'w');
fprintf(fid,' %d %d\n',nzeta,ntheta);
fprintf(fid,' %s\n',strtrim(vmec_data.input_extension));
fprintf(fid,' %f ',zeta);
fprintf(fid,'\n');
fprintf(fid,' %f ',theta);
fprintf(fid,'\n');
fprintf(fid,' %d\n',nsurfs);
fprintf(fid,' %d ',surfs');
fprintf(fid,'\n');
fclose(fid);

end

