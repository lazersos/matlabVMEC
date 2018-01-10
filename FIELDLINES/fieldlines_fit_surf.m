function fieldlines_fit_surf( data, surf )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


rpts = data.R_lines(surf,:);
zpts = data.Z_lines(surf,:);
ppts = data.PHI_lines(surf,:);

% Order by PHI
phimax = 2*pi/data.nfp;
ppts = mod(ppts,phimax);
n    = floor(data.nsteps/data.npoinc-1).*data.npoinc;
rpts = rpts(1:n);
zpts = zpts(1:n);
ppts = ppts(1:n);
order_array = [data.npoinc round(data.nsteps/data.npoinc)-1]
ppts = reshape(ppts,order_array);
rpts = reshape(rpts,order_array);
zpts = reshape(zpts,order_array);

% Order in a poloidal sense
r1(:,1) = rpts(:,1);
z1(:,1) = zpts(:,1);
rpts = rpts(:,2:end);
zpts = zpts(:,2:end);
n = size(rpts,2);
for i=2:n
    rt = r1(:,i-1);
    zt = z1(:,i-1);
    dr =[]; dz =[];
    for j=1:length(rt)
        dr(j,:) = rpts(j,:) - rt(j);
        dz(j,:) = zpts(j,:) - zt(j);
    end
    dl = sqrt(dr.^2+dz.^2);
    mins = min(dl');
    rpts2 =[]; zpts2=[];
    for j=1:length(mins)
        dex(j) = find(dl(j,:) == mins(j));
        r1(j,i) = rpts(j,dex(j));
        z1(j,i) = zpts(j,dex(j));
        if (dex(j) < size(rpts,2))
           rpts2(j,:) = rpts(j,[1:dex(j)-1 dex(j)+1:end]);
           zpts2(j,:) = zpts(j,[1:dex(j)-1 dex(j)+1:end]);
        else
            rpts2(j,:) = rpts(j,1:dex(j)-1);
            zpts2(j,:) = zpts(j,1:dex(j)-1);
        end
    end
    rpts = rpts2;
    zpts = zpts2;
%pause(0.01);
end

% Output a descur file
fid = fopen('descur.input','w');
fprintf(fid,'%3i %3i %3i\n',size(r1,2)-1,size(r1,1),data.nfp);
for i=1:size(rpts,1)
    temp=[r1(i,[1 3:size(r1,2)]); 0.0*r1(i,[1 3:size(r1,2)])+ppts(i,1); z1(i,[1 3:size(r1,2)])];
    fprintf(fid,'%20.10E %20.10E %20.10E\n',temp);
    fprintf(fid,'\n');
end
fclose(fid);

return;

% Now Fourier transform the quantities
rmn = abs(2*fft(r1',24)./size(r1,1));
zmn = abs(2*fft(z1',24)./size(z1,1));

% Now reformulate the plot
theta = 0:2*pi/359:2*pi;
zeta  = 0;
xm    = 0:12;
xn    = 0.*xm;
r2 = cfunct(theta,zeta,rmn(1:13,1),xm,xn);
z2 = sfunct(theta,zeta,zmn(1:13,1),xm,xn);
plot(r1(1,:),z1(1,:),'o-k');
hold on;
plot(r2,z2,'r');




end

