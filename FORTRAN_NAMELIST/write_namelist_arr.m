function f = write_namelist_arr(fid,name,val,varargin)
%WRITE_NAMELIST_ARR(fid,name,val,[index1min,index1max,index2min,index2max]) Write Array to namelist
%   This function writes a given array with value to a FORTRAN
%   namelist. x1 and y1 
%   fid     File ID to output to (opened by fopen)
%   name    String of variable name
%   val     Matrix.
%   Optional Arguments
%   index1min   Minimum index value for first index
%   index1max   Maximum index value for first index
%   index2min   Minimum index value for second index
%   index2max   Maximum index value for second index
%   If not specified index 1 starts goes from 1:max(val,1) and index 2
%   goes from 1:max(val,2)
%   
%   Output  
%           0 = successful;
%          -1 = Error;
%
%   Use
%           fileid=fopen('test');
%           iter=ones(10,20);
%           write_namelist_arr(fileid,'ITER',iter);
%           or
%           write_namelist_arr(fileid,'ITER',iter,0,9,-9,10);
%           fclose(fid);

if nargin == 3
    imin=1;
    imax=size(val,1);
    jmin=1;
    jmax=size(val,2);
elseif nargin == 7
    imin=varargin{1};
    imax=varargin{2};
    jmin=varargin{3};
    jmax=varargin{4};
else
    disp('Incorrect number of arguments!\n');
    disp('WRITE_NAMELIST_ARR(fid,name,val,[index1min,index1max,index2min,index2max])\n');
    f=-1;
    return
end
format='  %s(%s,%s) = %21.14E\n';
for i=1:size(val,1)
    for j=1:size(val,2);
        l=i-1+imin;
        m=j-1+jmin;
        fprintf(fid,format,upper(name),num2str(l),num2str(m),val(i,j));
    end
end
f=0;
return
end