function data = read_bnorm( file )
% READ_BNORM(filename) This function reads the BNORM output file.
% This function reads the BNORM output file into a data structure.
%
% Example usage
%      bnorm_data=read_vmec('bnorm.test');     % Reads VMEC wout file
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

% NOTES:
%


data=-1;
fid = fopen(file,'r');
if (fid < 0), return; end
temp=fscanf(fid,'%d %d %e',[3 inf]);
data.xm = temp(1,:);
data.xn = temp(2,:);
data.bnorm = temp(3,:)';

fclose(fid);

return;

end

