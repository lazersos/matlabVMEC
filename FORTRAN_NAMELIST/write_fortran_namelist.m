function write_fortran_namelist(fid,data,namelist)
%WRITE_FORTRAN_NAMELIST(data,namelist) Writes a structure to a
%   This function writes a FORTRAN input namelist to an open text file.
%
%   Example:
%       data=read_namelist('input.test','INDATA');
%       fid=fopen('input.test_out','w');
%       write_fortran_namelist(fid,data,'INDATA');
%       fclose(fid);
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           3/22/11

% Check the data field
if ~isstruct(data) || ~ischar(namelist) || (fid == -1)
    disp('Proper usage:');
    disp('   write_fortran_namelist(fid,data,namelist)');
    disp('      fid:       Open File ID');
    disp('      data:      Namelist Strucutre');
    disp('      namelist:  Namelist name');
    return
end
% Write namelist name
fprintf(fid,'%s\n',['&' upper(strtrim(namelist))]);
% Now write each element of the namelist
names=fieldnames(data);
nels=numel(names);
for i=1:nels
    val=data.(names{i});
    if ischar(val)
        fprintf(fid,'  %s = ''%s''\n',upper(names{i}),val);
    elseif (((size(val,1)==1) && (size(val,2)>1))...
            ||((size(val,2)==1) && (size(val,1)>1)))
        if isinteger(val(1))
            write_namelist_vec(fid,names{i},val,'int');
        elseif islogical(val(1))
            write_namelist_vec(fid,names{i},val,'boolean');
        else
            write_namelist_vec(fid,names{i},val);
        end
    elseif (numel(size(val))>=2) && (size(val,1)>1) && (size(val,2)>1)
        write_namelist_arr(fid,names{i},val);
    elseif islogical(val)
        write_namelist_boo(fid,names{i},val);
    elseif isinteger(val)
        write_namelist_int(fid,names{i},val);
    elseif isfloat(val)
        write_namelist_flt(fid,names{i},val);
    else
        disp('Unknown Type');
        disp(['     Variable: ' names{i}]);
        disp(['     Value:    ' num2str(val)]);
    end
end
% Write namelist end
fprintf(fid,'%s\n','/');
end