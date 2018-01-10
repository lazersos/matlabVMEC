function f = write_namelist_str(fid,name,val)
%WRITE_NAMELIST_STR(fid,'name',val) Write var to namelist
%   This function writes a varialbe of type string to a FORTRAN
%   namelist.
%   fid     File ID to output to (opened by fopen)
%   name    String of variable name
%   val     varialbe string value
%   
%   Output  
%           0 = successful;
%          -1 = Error;
%
%   Use
%           fileid=fopen('test');
%           iter='bob.file';
%           write_namelist_str(fileid,'ITER',iter);
%           fclose(fid);

format='  %s = ''%s''\n';
try
    fprintf(fid,format,upper(name),val);
    f=0;
catch fprintf_err
    f=-1;
end
return
end

