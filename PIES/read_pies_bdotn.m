function data = read_pies_bdotn()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Set Defaults
filename='fort.26';
line_format='%*10c %d %d %e\n';
nmin=0;
nmax=0;
mmax=0;
data=[];

fid=fopen(filename,'r');
% Read the array
line=fscanf(fid,line_format,[3 inf]);
mmax=max(line(1,:));
nmax=max(line(2,:));
nmin=min(line(2,:));
% Allocate the bdotn array
bdotn=zeros(mmax+1,2*nmax+1);
% Close the file
fclose(fid);
% Now process the data into bdotn
for i=1:size(line,2)
    bdotn(line(1,i)+1,line(2,i)+1-nmin)=line(3,i);
end
data=bdotn;

end

