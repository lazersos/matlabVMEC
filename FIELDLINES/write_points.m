function write_points(data,n,filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

phimax = max(data.phiaxis);
nfp = round(2*pi/max(data.phiaxis));
ntheta = round(data.nsteps/data.npoinc);
nzeta = data.npoinc;
fid = fopen(filename,'w');
fprintf(fid,'  %d  %d  %d\n',ntheta,nzeta,nfp);
dex = 1:data.npoinc:data.nsteps;
dex_len = length(dex);
outarray = [];
for i=1:nzeta
    dex = i:data.npoinc:data.nsteps-1;
    temp(1,:) = data.R_lines(n,dex);
    temp(2,:) = data.PHI_lines(n,dex);
    temp(3,:) = data.Z_lines(n,dex);
    outarray = [outarray temp];
end
outarray(2,:) = mod(outarray(2,:),phimax);
fprintf(fid,'  %20.12E  %20.12E  %20.12E\n',outarray);
fclose(fid);
end

