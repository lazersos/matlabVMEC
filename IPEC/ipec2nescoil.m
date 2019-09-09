function ipec2nescoil(iter,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Load data from files
idcon=read_ipec(['idcon_equil.' num2str(iter,'%5.5d') '.out']);
xdata=read_stellopt('xvec.dat');

% Extract Bnorm
iter2 = find(xdata.iter == iter);
bnorm=-xdata.x(:,iter2);
bnorm=[0:(length(bnorm)-1); ones(1,length(bnorm))*n; bnorm'];

% Extract plasma shape
mpol = 32;
r = idcon.r(:,end);
z = idcon.z(:,end);
Yr=fft(r,mpol)/(mpol-1);
Yr=abs(Yr(1:mpol/2+1));
Yz=fft(z,mpol)/(mpol-1);
Yz=abs(Yz(1:mpol/2+1));
xm=0:mpol/2;
ps=[xm' xm'.*0 Yr Yz];
ps(:,5:8)=0.0;



fid=fopen('bnorm.ipec','w');
fprintf(fid,'%4d %4d %22.10E\n',bnorm);
fclose(fid);

fid=fopen('nescin.ipec','w');
fprintf(fid,'------ Spatial dimensions ----\n');
fprintf(fid,'nu, nv, nu1, nv1, npol, ntor, lasym_bn\n');
fprintf(fid,' 64 64 64 64 64 10 F\n\n');
fprintf(fid,'------ Fourier Dimensions ----\n');
fprintf(fid,'mf, nf, md, nd (max in surf and bnorm files)\n');
fprintf(fid,' %3d %3d %3d %3d\n\n',size(bnorm,2)-1,1,20,20);
fprintf(fid,'------ Plasma Information ----\n');
fprintf(fid,'np     iota_edge       phip_edge       curpol\n');
fprintf(fid,' %3d %22.10E 0.0E+00 0.0E+00\n\n',1,idcon.q(end));
fprintf(fid,'------ Current Controls ----\n');
fprintf(fid,'cut  cup  ibex(=1,use fixed background coils)\n');
fprintf(fid,' 0.0E+00 1.0E+00 0\n\n');
fprintf(fid,'------ SVD controls -----\n');
fprintf(fid,'mstrt, mstep, mkeep, mdspw, curwt, trgwt\n');
fprintf(fid,'%2d %2d %2d %2d 0.0E+00 0.0E+00\n\n',0,0,0,4);
fprintf(fid,'------ Output controls -----\n');
fprintf(fid,'w_psurf w_csurf w_bnuv w_jsurf w_xerr w_svd\n');
fprintf(fid,'%2d %2d %2d %2d %2d %2d\n\n',0,0,0,0,0,0);
fprintf(fid,'------ Plasma Surface ----\n');
fprintf(fid,'Number of fourier modes in table\n');
fprintf(fid,'%4d\n',size(ps,1));
fprintf(fid,'Table of fourier coefficients\n');
fprintf(fid,'m,n,crc,czs,cls,crs,czc,clc\n');
fprintf(fid,'%2d %2d %22.10E %22.10E %22.10E %22.10E %22.10E %22.10E\n',ps');
fprintf(fid,'\n');
fprintf(fid,'------ Current Surface -----\n');
fprintf(fid,'Number of fourier modes in table\n');
fprintf(fid,'%4d\n',size(ps,1));
fprintf(fid,'Table of fourier coefficients\n');
fprintf(fid,'m,n,cr2c,czs2,crs2,czc2\n');
fclose(fid);


end

