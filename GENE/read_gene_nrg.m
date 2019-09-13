function data = read_gene_nrg( filename )
%READ_GENE_NRG Read a GENE nrg file.
%   Detailed explanation goes here


fid = fopen(filename,'r');
% Read first line
line=fgetl(fid);
time = sscanf(line,'%f',1);
line=fgetl(fid);
temp = sscanf(line,'%f',10);
n_spec = 0;
while numel(temp) > 1
   line=fgetl(fid);
   temp = sscanf(line,'%f',10);
   n_spec = n_spec + 1;
end
fclose(fid);
% Now read the data
fid = fopen(filename,'r');
data=[];
i=1;
while ~feof(fid)
    time = fscanf(fid,'%f',1);
    if isempty(time), break; end
    for j=1:n_spec
        temp=fscanf(fid,'%f',10);
        data(i,1,j) = time;
        data(i,2:11,j) = temp;
    end
    i = i + 1;
end

fclose(fid);

end

