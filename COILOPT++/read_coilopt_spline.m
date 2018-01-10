function data = read_coilopt_spline( filename )
%READ_COILOPT_SPLINE Reads a COILOPT++ spline file.
% This funciton reads the spline files generated and utilized by the
% COILOPT++ code.
%
% Example usage
%      data=read_vmec('fd.spline');     % Reads COILOPT++ spline file
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00


data=[];
% Number of Arguments
if (nargin == 0)
    disp('read_coilopt_spline requires a filename.');
    return
end
fid=fopen(filename,'r');
if fid < 0
    disp(['Error reading spline file: ' filename]);
    return;
end
temp=fscanf(fid,'%d',3);
data.nfp = temp(1);
data.ndoub = temp(2);
data.nsing = temp(3);
data.ncoil=temp(2)+temp(3);
for i=1:data.ncoil
    temp=fscanf(fid,'%d %e %e %e %d',[5]);
    data.nels(i) = temp(1);
    data.duu(i) = temp(1);
    data.dul(i) = temp(3);
    data.current(i) = temp(4);
    data.isym(i) = temp(5);
    knots = fscanf(fid,'%e',[data.nels(i)+4]);
    nodes = fscanf(fid,'%e',[2 data.nels(i)]);
    data.sp{i} = spmak(knots,nodes);
end
fclose(fid);

end

