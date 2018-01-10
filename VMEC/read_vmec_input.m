function data=read_vmec_input(filename)
% READ_VMEC_INPUT(handles) Reads the VMEC input file.
%   This function reads the VMEC input file and returns the contents
%   as a data structure.  This structure has an element for each variable
%   in the VMEC input namelist.
%   
%   Example:
%       loop_data=read_vmec_input('input.demo');
%
%   See also VMECedit.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.5
%   Date:           8/03/11

% Create the input namelist structure
data=read_namelist(filename,'INDATA');
% Merge the namelist with the default one
data=struct_merge(data,vmec_namelist_init('INDATA'));
% Return the data structure;
data.datatype='VMEC_input';
% Add some fields
if ~isfield(data,'ntheta')
    data.ntheta=2*data.mpol+6;
end
if data.ntheta == 0
    data.ntheta=2*data.mpol+6;
end
if ~isfield(data,'nzeta')
    data.nzeta=2*data.ntor+4;
end
if data.nzeta == 0
    data.nzeta=2*data.ntor+4;
end
return;
end

