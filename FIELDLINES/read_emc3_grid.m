function data = read_emc3_grid(file)
%data = read_emc3_grid(file) Returns data from a grid3d.dat file.
%   Detailed explanation goes here


fid=fopen(file,'r');
for j=1:1
    temp=fscanf(fid,'%d %d %d',3);
    data.nr = temp(1);
    data.ntheta = temp(2);
    data.nzeta = temp(3);
    for i=1:data.nzeta
        temp=fscanf(fid,'%f',1);
        data.phi(i,j)=temp;
        temp=fscanf(fid,'%f',[data.nr.*data.ntheta]);
        data.R(:,i,j)=temp;
        temp=fscanf(fid,'%f',[data.nr.*data.ntheta]);
        data.Z(:,i,j)=temp;
    end
end
fclose(fid);
end

