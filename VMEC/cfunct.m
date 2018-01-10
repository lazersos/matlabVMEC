function f=cfunct(theta,zeta,rmnc,xm,xn,varargin)
% CFUNCT(theta,zeta,rmnc,xm,xn) Cosine Fourier Transform
% This function returns the cosine Fourier tranform of a value for a 
% given theta and zeta (coordinates over the transform).  It assumes a 
% form of the transform where:
% f=\sum_{m=0,M}\sum_{n=-N,N} rbc(m,n)*cos(m*theta-n*nfp*zeta)
% the VMEC Fourier transform.
%
% f=CFUNCT(theta,zeta,rmnc,xm,xn)
% Inputs
% theta:    poloidal coordinate (0,2pi)
% zeta:     toroidal coordinate (0,2pi)
% rmnc:     Array of Fouier Coefficients (mn,ns)
% xm:       Poloidal mode number array (mn)
% xn:       Toroidal mode number array (mn)
%
% Exmaple Usage
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      r=cfunct(theta,zeta,data.rmnc,data.xm,data.xn);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       2.00

% 10/17/11 - SAL  Vectorized Version for increased speed.

ns=size(rmnc,2);
lt=length(theta);
lz=length(zeta);
% Create mode x angle arrays
mt=xm'*theta;
nz=xn'*zeta;
% Create Trig Arrays
cosmt=cos(mt);
sinmt=sin(mt);
cosnz=cos(nz);
sinnz=sin(nz);
% Calcualte the transform
f=zeros(ns,lt,lz);
for k=1:ns
    rmn=repmat(rmnc(:,k),[1 lt]);
    f(k,:,:)=(rmn.*cosmt)'*cosnz-(rmn.*sinmt)'*sinnz;
end
return
end