function data = read_pies_stoch(filename)
%READ_PIES_STOCH Reads the PIES 'stoch_data' file into a structure
%   READ_PIES_STOCH(filename)   This subroutine reads the 'stoch_data' file
%   produced by PIES when calculating finite pressure gradients in
%   stochastic regions and returns the contents as a structure.  
%
% Example usage
%      data=read_pies_stoch('pies.stoch_data');
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

% Set some defaults
fid=fopen(filename,'r');
temp=fscanf(fid,'%d',3);
data.k=temp(1);
data.m=temp(2);
data.n=temp(3);
data.rho0=fscanf(fid,'%g',[data.m+1 2*data.n+1]);
data.th0=fscanf(fid,'%g',[data.m+1 2*data.n+1]);
temp=fscanf(fid,'%d',4);
data.k0=temp(1);
data.n0=temp(2);
data.m0=temp(3);
data.ir_grid_active=temp(4);
if data.ir_grid_active
    data.rho_irreg=fscanf(fid,'%g',[1 data.k0+1]);
end
data.iota=fscanf(fid,'%g',[1 data.k0+1]);
data.mu=fscanf(fid,'%g',(data.m0+1)*(2*data.n0+1)*(data.k0+1));
data.mu=reshape(data.mu,[(data.m0+1) (2*data.n0+1) (data.k0+1)]);
data.mustochf=fscanf(fid,'%d',1);
fclose(fid);
data.datatype='pies_stoch';

end

