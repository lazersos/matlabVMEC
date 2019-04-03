function f = read_cobra(filename)
%READ_COBRA(filename) Reads a COBRA output file into a COBRA structure.
%
%   READ_COBRA(filename) Reads a COBRA output file into a COBRA structure.
%
%   Example:
%       cobra_data=read_cobra('grate_cobra.test');
%
%   Written by: Samuel Lazerson (lazerson@pppl.gov)
%   Version:    1.5
%   Date:       09/28/2011

fid=fopen(filename);
data=fscanf(fid,'%f %f %d\n',3);
zeta=data(1);
theta=data(2);
ns=data(3);
%data=fscanf(fid,'%d %f\n',[2 ns]); % old way
data=fscanf(fid,'%d %f %f\n',[3 ns]);
surfs=data(1,:);
grate=data(3,:);
while ~feof(fid)
    data=fscanf(fid,'%f %f %d\n',3);
    if isempty(data); break; end
    zeta=[zeta data(1)];
    theta=[theta data(2)];
    data=fscanf(fid,'%d %f %f\n',[3 ns]);
    grate=[grate; data(3,:)];
end
fclose(fid);
f.grate = grate;

% Now parse into ns by theta by zeta arrays
% grate=permute(grate,[2 1]);
% i_axis=unique(theta);
% j_axis=unique(zeta);
% if (zeta(length(j_axis)+1)==zeta(length(j_axis)+2))
%     j_axis=[j_axis 36];
% end
% f.grate=zeros(ns,length(i_axis),length(j_axis));
% for k=1:size(grate,2)
%     f.grate(:,i_axis==theta(k),j_axis==zeta(k))=grate(:,k);
% end
% 
% f.zeta=zeta;
% f.theta=theta;
% f.ns=ns;
% f.surfs=surfs;

end

