function f=cfunct_old(theta,zeta,rbc,nfp,varargin)
% CFUNCT_OLD(theta,zeta,rbc,nfp) Cosine Fourier Transform
% This function returns the cosine Fourier tranform of a value for a 
% given theta and zeta (coordinates over the transform).  It assumes a 
% form of the transform where:
% f=\sum_{m=0,M}\sum_{n=-N,N} rbc(m,n)*cos(m*theta-n*nfp*zeta)
% the VMEC Fourier transform.
%
% f=CFUNCT_OLD(theta,zeta,zbc,nfp)
% Inputs
% theta:    poloidal coordinate (0,2pi)
% zeta:     toroidal coordinate (0,2pi)
% rbc:      Array of Fouier Coefficients (ns,mpol,ntor)
% nfp:      Number of field Periods
%
% If zeta is evenly divisible by nfp then we only calc over the field
% period subdomain and copy to the rest of the array.  If theta is a 3D
% array will use it (this is done for VMEC spectral condenstation).
%
% Exmaple Usage
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      r=cfunct_old(theta,zeta,data.rbc,data.nfp);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.50

% 08/18/11 - SAL  Modified to allow Boozer zeta parameterization

index1=size(rbc,1);
index2=size(rbc,2);
index3=size(rbc,3);
zs=size(zeta,2);
ts=size(theta,2);
% Check to see if zs is divisible by nfp otherwise calculate full array
% Also check to see you want a full toroid
if ~(length(size(zeta))==3)
    maxzeta=max(zeta);
else
    maxzeta=max(max(max(zeta)));
end
% Transform Zeta to 3D
if ~(length(size(zeta))==3)
    zeta=repmat(repmat(zeta',[1 ts])',[1 1 index3]);
else
    zs=size(zeta,3);
    zeta=permute(zeta,[2 3 1]);
end
% Calculate zeta bounds
if ~mod(zs-1,nfp) && (nfp ~= 1) && (maxzeta == 2*pi)
    zbound=(zs-1)/nfp;
else
    zbound=zs;
end
% Expand theta and zeta into theta,zeta,ns sized grid
% Check theat to see if theta 3D
if ~(length(size(theta))==3)
    theta=repmat(repmat(theta',[1 zs]),[1 1 index3]);
else
    theta=permute(theta,[2 3 1]);
end
f=zeros(ts,zs,index3);
% Calculate the transform
for m=1:index1
    for n=1:index2
        for i=1:index3
            ndex=(n-(index2-1)/2-1)*nfp;
            mdex=m-1;
            f(:,1:zbound,i)=f(:,1:zbound,i)+rbc(m,n,i).*cos(mdex.*theta(:,1:zbound,i)+ndex.*zeta(:,1:zbound,i));
        end
    end
end
% Now fill f if zs was divisible by nfp copy it to the rest of the array
if ~mod(zs-1,nfp) && (nfp ~= 1) && (maxzeta == 2*pi)
    temp=repmat(permute(f(:,1:zbound,:),[3 1 2]),[1 1 nfp]);
    f=zeros(index3,ts,zs);
    f(:,:,1:zs-1)=temp;
    f(:,:,zs)=f(:,:,1);
elseif (nfp == 1)
    f=permute(f,[3 1 2]);
    f(:,:,zs)=f(:,:,1);
else
    f=permute(f,[3 1 2]);
end
return
end