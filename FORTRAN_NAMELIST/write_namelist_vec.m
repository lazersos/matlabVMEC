function f = write_namelist_vec(fid,name,val,varargin)
%WRITE_NAMELIST_VEC(fid,'name',val[,'vert']) Write (float) vec to namelist
%   This function writes a given vector to a FORTRAN
%   namelist.  Default type is float.
%   fid         File ID to output to (opened by fopen)
%   name        String of variable name
%   val         Value assigned to variable.
%   Optional Arguments
%   'int'       Write integer array
%   'boolean'   Write Boolean array
%   'vert'      Write vector element by element (ex. TEST(1) = 5.0)
%   Output
%           0 = successful;
%          -1 = Error;
%
%   Use
%           fileid=fopen('test');
%           iter=1:10;
%           write_namelist_vec(fileid,'ITER',iter);
%           fclose(fid);

% Handle varargin
vert=0;
intflag=0;
booflag=0;
if nargin > 3
    for i=1:nargin-3
        switch varargin{i}
            case 'vert'
                vert=1;
            case 'int'
                intflag=1;
            case 'boolean'
                booflag=1;
        end
    end
end
f=1;
length=max(size(val));
j=1;
% Write Variable Name
try
    if intflag
        if vert
            for i=1:length
                fprintf(fid,'  %s(%3d) = %s\n',upper(name),i,num2str(val(i)));
            end
        else
            fprintf(fid,'  %s = ',upper(name));
            for i=1:length-1
                if j == 10
                    fprintf(fid,' %s\n ',num2str(val(i)));
                    j=1;
                else
                    fprintf(fid,' %s ',num2str(val(i)));
                    j=j+1;
                end
            end
            fprintf(fid,' %s\n',num2str(val(length)));
        end
    elseif booflag
        if vert
            for i=1:length
                if (val(i) == 1)
                    fprintf(fid,'  %s(%3d) = %s\n',upper(name),i,'T');
                else
                    fprintf(fid,'  %s(%3d) = %s\n',upper(name),i,'F');
                end
            end
        else
            fprintf(fid,'  %s = ',upper(name));
            for i=1:length-1
                if j == 10
                    if (val(i) == 1)
                        fprintf(fid,' %s\n ','T');
                    else
                        fprintf(fid,' %s\n ','F');
                    end
                    j=1;
                else
                    if (val(i) == 1)
                        fprintf(fid,' %s ','T');
                    else
                        fprintf(fid,' %s ','F');
                    end
                    j=j+1;
                end
            end
            if (val(length) == 1)
                fprintf(fid,' %s\n','T');
            else
                fprintf(fid,' %s\n','F');
            end
        end
    else
        if vert
            for i=1:length
                fprintf(fid,'  %s(%3d) = %21.14E\n',upper(name),i,val(i));
            end
        else
            fprintf(fid,'  %s = ',upper(name));
            for i=1:length-1
                if j == 3
                    fprintf(fid,'%21.14E\n ',val(i));
                    j=1;
                else
                    fprintf(fid,'%21.14E ',val(i));
                    j=j+1;
                end
            end
            fprintf(fid,'%21.14E\n',val(length));
        end
    end
    f=0;
catch fprintf_err
    f=-1;
end
return
end
