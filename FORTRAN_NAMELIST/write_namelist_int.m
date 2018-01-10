function f = write_namelist_int(fid,name,val)
%WRITE_NAMELIST_INT(fid,'name',val) Write var to namelist
%   This function writes a given variable with value to a FORTRAN
%   namelist.  In an integer way.
%   fid     File ID to output to (opened by fopen)
%   name    String of variable name
%   val     Value assigned to variable.
%   
%   Output  
%           0 = successful;
%          -1 = Error;
%
%   Use
%           fileid=fopen('test');
%           iter=100;
%           write_namelist_int(fileid,'ITER',iter);
%           fclose(fid);

format='  %s = %s\n';
try
    fprintf(fid,format,upper(name),num2str(val));
    f=0;
catch fprintf_err
    f=-1;
end
return
end

