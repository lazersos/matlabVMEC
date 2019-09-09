function data = read_pies_magco(m,n,ns,varargin)
%READ_PIES_MAGCO Read the magco.out file and return rbc and zbs arrays
%   READ_PIES_MAGCO(m,n,ns) Reads the 'magco.out' in the current directory and
%   returns the RBC and ZBS like arrays.  Here m is the number of PIES
%   poloidal modes, n is the number to PIES toroidal modes, and ns is the
%   number of PIES surfaces.  These values should be set the same as those
%   found in the PIES input file.  The 'magco.out' file is produced on the
%   laster iteration by setting IWRTMG=1 in the EXLSTA input namelist.
%
%   READ_PIES_MAGCO(m,n,ns,'filename',filename)  Allows the user to specify
%   the specific filename via a string filename.
%
%   READ_PIES_MAGCO(m,n,ns,'rmaj',rmaj)  Allows the user to specify the
%   rmaj value thus returning the rbc array in VMEC form (as opposed to
%   PIES attempting to center the coordinates).  
%
% Example usage
%       pies_data=read_pies_netcdf('pies_save.nc');     % Reads the PIES output file
%       pies_magco=read_pies_maco(pies_data.m,pies_data.n,pies_data.ns,'rmaj',pies_data.Rmajor)
%       ntheta=179;
%       nzeta=30;
%       theta=0:2*pi/ntheta:2*pi;
%       zeta=0:2*pi/(nzeta*pies_data.nfp):2*pi/pies_data.nfp;
%       r_magco=cfunct(theta,zeta,pies_magco.rbc,pies_data.nfp);
%       z_magco=sfunct(theta,zeta,pies_magco.zbs,pies_data.nfp);
%       toroslice(r_magco,1,z_magco,1:pies_data.ns);
%
%   See also:  read_pies_netcdf, cfunct, sfunct, toroslice.
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0

% Defaults
filename='magco.out';
data=[];
offset=3;
rmaj=0.0;
i=1;

% Handle varargin
while i <= nargin-offset
    switch varargin{i}
        case 'filename'
            i=i+1;
            filename=varargin{i};
        case 'rmaj'
            i=i+1;
            rmaj=varargin{i};
        case 'hitsrf'
            i=i+1;
            hitsrf=varargin{i};
    end
    i=i+1;
end

fid=fopen(filename,'r');
if n > 0
    rbc=zeros(m+1,2*n+1,ns+1);
    zbs=zeros(m+1,2*n+1,ns+1);
    for i=1:ns+1
        for j=1:2*n+1
            rbc(:,j,i)=fscanf(fid,'%e',m+1);
        end
        for j=1:2*n+1
            zbs(:,j,i)=fscanf(fid,'%e',m+1);
        end
    end
    rbc(1,n+1,:)=rbc(1,n+1,:)+rmaj;
else
    rbc=zeros(m+1,3,ns+1);
    zbs=zeros(m+1,3,ns+1);
    for i=1:ns+1
        rbc(:,2,i)=fscanf(fid,'%e',m+1);
        zbs(:,2,i)=fscanf(fid,'%e',m+1);
    end
    rbc(1,2,:)=rbc(1,2,:)+rmaj;
    n=1;
end
fclose(fid);

% Need to convert to VMEC convention
mnmax=(m+1)*(2*n+1);
rmnc=zeros(mnmax,ns+1);
zmns=zeros(mnmax,ns+1);
k=1;
for i=1:m+1
    for j=1:2*n+1
        xm(k)=i-1;
        xn(k)=j-n-1;
        rmnc(k,:)=rbc(i,j,:);
        zmns(k,:)=zbs(i,j,:);
        k=k+1;
    end
end
data.rmnc=rmnc;
data.zmns=zmns;
data.xm=double(xm);
data.xn=double(xn);
data.datatype='pies_magco';
data.ns = size(rmnc,2);
data.hitsrf = hitsrf;
return

end

