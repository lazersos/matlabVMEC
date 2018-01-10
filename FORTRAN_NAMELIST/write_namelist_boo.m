function f = write_namelist_boo(fid,name,val)
%WRITE_NAMELIST_INT(fid,'name',val) Write var to namelist
%   This function writes a given variable with value to a FORTRAN
%   namelist.  In a boolean way.
%   fid     File ID to output to (opened by fopen)
%   name    String of variable name
%   val     Value assigned to variable.
%   
%   Output  
%           0 = successful;
%           -1 = Error;
%
%   Use
%           fileid=fopen('test');
%           iter=0;  %1: TRUE 2: FALSE
%           write_namelist_int(fileid,'ITER',iter);
%           fclose(fid);
f=1;
format='  %s = %s\n';
switch val
    case 1
        fprintf(fid,format,upper(name),'T');
    case 0
        fprintf(fid,format,upper(name),'F');
    case default
        f=-1;
end
end

