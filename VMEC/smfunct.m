function f=smfunct(theta,zeta,zbs,nfp)
% SMFUNCT(theta,zeta,zbc,nfp) Sine Fourier Transform for derivatives.
% This function returns the sine Fourier tranform of a value for a given
% theta and zeta (coordinates over the transform).  It assumes a form of
% the transform where:
% f=\sum_{m=0,M}\sum_{n=-N,N} zbs(m,n)*m*sin(m*theta-n*nfp*zeta)
% the VMEC Fourier transform.
%
% f=SMFUNCT(theta,zeta,zbc,nfp)
% Inputs
% theta:    poloidal coordinate (0,2pi)
% zeta:     toroidal coordinate (0,2pi)
% zbs:      Array of Fouier Coefficients (ns,mpol,ntor)
% nfp:      Number of field Periods
% If zeta is evenly divisible by nfp then we only calc over the field
% period subdomain and copy to the rest of the array.  If theta is a 3D
% array will use it (this is done for VMEC spectral condenstation).
%
% Exmaple Usage
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      z=smfunct(theta,zeta,data.zbs,data.nfp);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.05

index1=size(zbs,1);
index2=size(zbs,2);
index3=size(zbs,3);
zs=size(zeta,2);
ts=size(theta,2);
f=zeros(ts,zs,index3);
% Check to see if zs is divisible by nfp otherwise calculate full array
% Also check to see you want a full toroid
maxzeta=max(zeta);
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
zeta=repmat(repmat(zeta',[1 ts])',[1 1 index3]);
% Calculate the transform
for m=1:index1
    for n=1:index2
        for i=1:index3
            ndex=(n-(index2-1)/2-1)*nfp;
            mdex=m-1;
            f(:,1:zbound,i)=f(:,1:zbound,i)-zbs(m,n,i).*mdex.*sin(mdex.*theta(:,1:zbound,i)-ndex.*zeta(:,1:zbound,i));
        end
    end
end
% Now fill f if zs was divisible by nfp copy it to the rest of the array
if ~mod(zs-1,nfp) && (nfp ~= 1) && (maxzeta == 2*pi)
    temp=repmat(permute(f(:,1:zbound,:),[3 1 2]),[1 1 nfp]);
    f=zeros(index3,ts,zs);
    f(:,:,1:zs-1)=temp;
    f(:,:,zs)=f(:,:,1);
else
    f=permute(f,[3 1 2]);
end
return
end