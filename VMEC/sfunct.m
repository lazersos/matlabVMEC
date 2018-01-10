function f=sfunct(theta,zeta,zmns,xm,xn,varargin)
% SFUNCT(theta,zeta,zmns,xm,xn) Sine Fourier Transform
% This function returns the cosine Fourier tranform of a value for a 
% given theta and zeta (coordinates over the transform).  It assumes a 
% form of the transform where:
% f=\sum_{m=0,M}\sum_{n=-N,N} rbc(m,n)*cos(m*theta-n*nfp*zeta)
% the VMEC Fourier transform.
%
% f=SFUNCT(theta,zeta,zmns,xm,xn,)
% Inputs
% theta:    poloidal coordinate (0,2pi)
% zeta:     toroidal coordinate (0,2pi)
% zmns:     Array of odd Fouier Coefficients (mn,ns)
% xm:       Poloidal mode number array (mn)
% xn:       Toroidal mode number array (mn)
%
% Exmaple Usage
%      theta=0:2*pi/36:2*pi;
%      zeta=0:2*pi/36:2*pi;
%      data=read_vmec('wout.test');
%      r=sfunct(theta,zeta,data.rmnc,data.xm,data.xn);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       2.00

% 10/17/11 - SAL  Vectorized Version for increased speed.

ns=size(zmns,2);
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
    zmn=repmat(zmns(:,k),[1 lt]);
    f(k,:,:)=(zmn.*sinmt)'*cosnz+(zmn.*cosmt)'*sinnz;
end
return
end