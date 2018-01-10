function [varargout]=DIAGNOplot(filename,varargin)
%DIAGNOPLOT(filename[,options]) Plots DIAGNO input files
%   DIAGNOPLOT(filename) plots the potential on the flux surface from a 
%   diagno input file.  Returns a handle to the patch surfaces produced by
%   isotoro.
%   Options:
%   DIAGNOPLOT(filename,'ntheta',n)     n is the number of points to use in
%                                       theta.
%   DIAGNOPLOT(filename,'nzeta',n)      n is the number of points to use in
%                                       zeta.
%

ntheta=36;
nzeta=36;
% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case 'ntheta'
                i=i+1;
                ntheta=varargin{i}-1;
            case 'nzeta'
                i=i+1;
                nzeta=varargin{i}-1;
        end
    end
end
theta=0:2*pi/(ntheta-1):2*pi;
zeta=0:2*pi/(nzeta-1):2*pi;
fid=fopen(filename,'rt');   %Open file
% diagno_in files have a specific format
input_type=fgetl(fid);      %Code Type
if ~strcmp(input_type,'vmec2000')
    disp('ERROR:  Not a VMEC2000 file!');
    return
end
fgetl(fid);                 %Header
line=fgetl(fid);            %nfp mpol ntor
temp=sscanf(line,'%d',3);
nfp=temp(1);
mpol=temp(2);
ntor=temp(3);
ns=(mpol-1)*(2*ntor+1)+(ntor+1);
np=(mpol+2)*(2*ntor+1);
rmnc=zeros(1,ns);
zmns=zeros(1,ns);
potsin=zeros(1,np);
fgetl(fid);                 %Header
for i=1:ns                  %RMNC
        line=fgetl(fid);
        rmnc(i)=sscanf(line,'%e');
end
fgetl(fid);                 %Header
for i=1:ns                  %ZMNS
        line=fgetl(fid);
        zmns(i)=sscanf(line,'%e');
end
fgetl(fid);                 %Header
for i=1:np                  %POTSIN
        line=fgetl(fid);
        potsin(i)=sscanf(line,'%e');
end
potsin=potsin*.5;
fgetl(fid);                 %Header
line=fgetl(fid);
phiedge=sscanf(line,'%f');
fgetl(fid);                 %Header
line=fgetl(fid);
ncur=sscanf(line,'%d');
fgetl(fid);                 %Header
line=fgetl(fid);
extcur=sscanf(line,'%f',[ncur]);
fgetl(fid);                 %Header
line=fgetl(fid);
torcur=sscanf(line,'%e');
fgetl(fid);                 %Header
line=fgetl(fid);
raxis=sscanf(line,'%e',ntor+1);
fgetl(fid);                 %Header
line=fgetl(fid);
zaxis=sscanf(line,'%e',ntor+1);
% end reading file
fclose(fid);                %Close file
% Reformat Arrays to rbc zbs potbs style
msize=mpol+2;
nsize=2*ntor+1;
rbc=zeros(msize,nsize,1);
zbs=zeros(msize,nsize,1);
potbs=zeros(msize,nsize,1);
offset=ntor+1;
i=1;
m=1;
for n=offset:nsize
    rbc(m,n)=rmnc(i);
    zbs(m,n)=zmns(i);
    i=i+1;
end
for m=2:mpol
    for n=1:nsize
        rbc(m,n)=rmnc(i);
        zbs(m,n)=zmns(i);
        i=i+1;
    end
end
i=1;
for n=1:(2*ntor+1)
    for m=1:mpol+2
        potbs(m,n)=potsin(i);
        i=i+1;
    end
end
r=cfunct(theta,zeta,rbc,nfp);
z=sfunct(theta,zeta,zbs,nfp);
pot=sfunct(theta,zeta,potbs,nfp);
hpatch=isotoro(r,z,zeta,[1],pot);
title('VMEC-DIAGNO Potential');
% Handle outputs
varargout{1}=hpatch;
end

