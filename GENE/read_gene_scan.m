function [data_array varargout] = read_gene_scan( input_args )
%READ_GENE_SCAN Reads a scan file.
%   The READ_GENE_SCAN function reads a GENE scan file from the
%   scanfiles0000 directory.  It can return just the growthrate or multiple
%   values:
%       data = read_gene_scan; % Growth Rate
%       [data x y] = read_gene_scan; % Growth Rate, axes
%       [data x y real] = read_gene_scan; % Growth Rate, axes, imaginary

dir = 'scanfiles0000';
fid = fopen([dir '/scan.log'],'r');
temp = fgetl(fid) % Read the header;
temp = fgetl(fid);
ndelim = numel(strfind(temp,'|'));
scan_str = ' %d ';
for i=1:ndelim-1
    scan_str = [scan_str '| %e '];
end
scan_str = [scan_str '| %f %f'];
data = [];
while temp ~= -1
    vals = sscanf(temp,scan_str);
    data = [data; vals'];
    temp = fgetl(fid);
end
fclose(fid);

num_dims=size(data,2)-3;
num_gr=size(data,2)-1;

array_size = [];
for i=1:num_dims
    array_size = [array_size sum(data(:,i+1)==data(1,i+1))];
end
data_array = reshape(data(:,num_gr),fliplr(array_size));
nout = max(nargout,1) - 1;
if (nout > 1)
    varargout{1} = reshape(data(:,2),fliplr(array_size));
    varargout{2} = reshape(data(:,3),fliplr(array_size));
end
if (nout > 2)
    varargout{3} = reshape(data(:,num_gr+1),fliplr(array_size));
end

%for k=1:nout
%    varargout{k} = reshape(data(:,num_gr),fliplr(array_size));
%end

end

