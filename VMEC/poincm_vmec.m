function poincm_vmec(varargin)
%POINCM_VMEC(cdata) Creates a Poincare plot from VMEC and Coil data
%   POINCM_VMEC(vdata,cdata) Creates a Poincare plot from VMEC and coil
%   data.

% Set defaults
vmecdata=[];
coildata=[];
i=1;
ns=5;
npunk=100;
maxiter=1000;
dt=0.01;
phiend=2*pi/10;
% Handle varargin
while i <= nargin
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch varargin{i}.datatype
                case 'wout'
                    vmecdata=varargin{i};
                case 'coil_data'
                    coildata=varargin{i};
            end
        end
    elseif ischar(varargin{i})
        switch varargin{i}
            case 'filename'
                i=i+1;
                filename=varargin{i};
        end
    end
    i=i+1;
end
% Do some error checking
if isempty(vmecdata) && isempty(coildata)
    disp('ERROR:  You must supply a Coil or VMEC structure');
    return
end
if isempty(coildata)
    ns=vmecdata.ns;
end

% Begin Program
fieldlines=zeros(ns,npunk,2); % Array of punctures
% Setup starting points for field line following
if ~isempty(vmecdata)
    for i=1:ns
        fieldlines(i,1,1)=cfunct(0,0,vmecdata.rbc(:,:,i),vmecdata.nfp);
    end
    fieldlines(:,1,2)=0.0;
end
pause(1.0);
% Follow field lines
for i=1:ns
    zeta=0.0;
    theta=0.0;
    nhit=1;
    r=cfunct(theta,zeta,vmecdata.rbc(:,:,i),vmecdata.nfp)-vmecdata.Rmajor;
    z=sfunct(theta,zeta,vmecdata.zbs(:,:,i),vmecdata.nfp);
    hold on
    plot(r,z,'xb');
    hold off
    while nhit <= npunk
        nstep=0;
        while nstep < maxiter
            bu=cfunct(theta,zeta,vmecdata.buc(:,:,i),vmecdata.nfp);
            bv=cfunct(theta,zeta,vmecdata.bvc(:,:,i),vmecdata.nfp);
            drdu=sfunct(theta,zeta,vmecdata.rus(:,:,i),vmecdata.nfp);
            drdv=sfunct(theta,zeta,vmecdata.rvs(:,:,i),vmecdata.nfp);
            dzdu=cfunct(theta,zeta,vmecdata.zuc(:,:,i),vmecdata.nfp);
            dzdv=cfunct(theta,zeta,vmecdata.zvc(:,:,i),vmecdata.nfp);
            br=bu*drdu+bv*drdv;
            bphi=r*bv;
            bz=bu*dzdu+bv*dzdv;
            bx=br*cos(zeta)-bphi*sin(zeta);
            by=br*sin(zeta)+bphi*cos(zeta);
            bnorm=sqrt(bx*bx+by*by+bz*bz);
            x=r*cos(zeta)+dt*bx/bnorm;
            y=r*sin(zeta)+dt*by/bnorm;
            z=z+dt*bz/bnorm;
            r=sqrt(x*x+y*y);
            theta=atan(z/r);
            zeta=atan(y/x);
            r=cfunct(theta,zeta,vmecdata.rbc(:,:,i),vmecdata.nfp)-vmecdata.Rmajor;
            z=sfunct(theta,zeta,vmecdata.zbs(:,:,i),vmecdata.nfp);
            hold on
            plot(r,z,'.k');
            pause(.1);
            hold off
            nstep=nstep+1;
            if abs(zeta) > phiend
                fieldlines(i,nhit,1)=r+vmecdata.Rmajor;
                fieldlines(i,nhit,2)=z;
                zeta=0.0;
                nstep
                nstep=maxiter;
            end
        end
        nhit=nhit+1;
    end
end
return
end
