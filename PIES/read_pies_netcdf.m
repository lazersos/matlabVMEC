function f = read_pies_netcdf(varargin)
%READ_PIES_NETCDF(filename) Read NetCDF data from PIES into VMEC-like
%structure
%   data=READ_PIES_NETCDF(filename)
%   This function reads the PIES output file and returns the data in a
%   structure similar to that of READ_VMEC.  Note that you must pass -theta
%   to the cfunct and sfunct routines for PIES data.  This is done to
%   handle the PIES kernel which is (n*zeta-m*theta), as we use the NESCOIL
%   convention in the matlabVMEC suite (m*theta+n*zeta).
%   Written by: Samuel Lazerson (lazerson@pppl.gov)
%   Created:    06/11/10
%   Last Edit:  08/02/11 -SAL

%   08/02/11 - SAL Edited so we covert relbup_old into B field not Bup.

%   Do some Error Checking
if nargin ==1
    filename=varargin{1};
else
    return
end

% Read file with read_netcdf
f=read_netcdf(filename);
% Scalars
f.Rmajor=f.rmaj;
f.nfp=f.nper;
f.mpol=f.m;
f.ns=size(f.x,2);
f.ntor=f.n;
f.datatype='pies_out';
f.nu=2.*(f.mpol+1)+6;
f.ierr_vmec=0;
f.input_extension=filename;
f.mgrid_file='coil_data';
f.mode=permute(f.mode,[2 1]);
f.jcurv=f.iprim0.*(2*pi/f.beta);
% Now create VMEC-like arrays
f.rbc=permute(f.x(1,:,:,:),[4 3 2 1]);
mid=(size(f.rbc,2)+1)/2;
f.rbc(1,mid,:)=f.rbc(1,mid,:)+f.rmaj;               % Do this so R not r
f.zbs=permute(f.x(2,:,:,:),[4 3 2 1]);
f.bsc=permute(f.relbup_old(1,:,:,:),[4 3 2 1]);
f.buc=permute(f.relbup_old(2,:,:,:),[4 3 2 1]);
f.bvc=permute(f.relbup_old(3,:,:,:),[4 3 2 1]);
f.rho_mag=permute(f.ooldrh,[3 2 1]);
f.tht_mag=permute(f.ooldth,[3,2,1]);
% Now we preserve the shape of the arrays by padding n=1 if n=0
% if (f.ntor == 0)
%     names={'rbc' 'zbs' 'bsc' 'buc' 'bvc'};
%     for i=1:length(names)
%         temp=f.(names{i});
%         f.(names{i})=zeros(f.mpol+1,3,f.ns);
%         f.(names{i})(:,2,:)=temp;
%     end
%     f.ntor=1;
%     % Handle mode seperately
%     temp=f.mode;
%     f.mode=zeros(f.mpol+1,3);
%     f.mode(:,2)=temp;
% end
% Now create the xm and xn arrays and vector arrays
f.mn=(f.m+1)*(2*f.n+1);
f.xm=zeros(1,f.mn);
f.xn=zeros(1,f.mn);
f.xm_mode=zeros(1,f.mn);
f.xn_mode=zeros(1,f.mn);
k=1;
for i=0:f.m
    for j=-f.n:f.n
        f.xm(k)=i;
        f.xn(k)=j.*f.nfp;
        f.xm_mode(k)=i.*f.mode(i+1,j+f.n+1);
        f.xn_mode(k)=j.*f.nfp.*f.mode(i+1,j+f.n+1);
        k=k+1;
    end
end
i=1;
f.rmnc=zeros(f.mn,f.ns);
f.zmns=zeros(f.mn,f.ns);
f.bsmnc=zeros(f.mn,f.ns);
f.bumnc=zeros(f.mn,f.ns);
f.bvmnc=zeros(f.mn,f.ns);
for j=1:f.mpol+1
    for k=1:2*f.ntor+1
        f.rmnc(i,:)=f.rbc(j,k,:);
        f.zmns(i,:)=f.zbs(j,k,:);
        f.bsmnc(i,:)=f.bsc(j,k,:);
        f.bumnc(i,:)=f.buc(j,k,:);
        f.bvmnc(i,:)=f.bvc(j,k,:);
        i=i+1;
    end
end
f.rumns=-f.rmnc.*repmat(f.xm',[1 f.ns]);
f.rvmns=-f.rmnc.*repmat(f.xn',[1 f.ns]);
f.zumnc=f.zmns.*repmat(f.xm',[1 f.ns]);
f.zvmnc=f.zmns.*repmat(f.xn',[1 f.ns]);
% J||
if isfield(f,'mu')
    f.jlc=permute(f.mu,[3 2 1]);
    f.jlc_title='J_{||}';
    f.jlmnc=zeros(f.mn,f.ns);
    i=1;
    for j=1:f.mpol+1
        for k=1:2*f.ntor+1
            f.jlmnc(i,:)=f.jlc(j,k,:)./f.beta;
            i=i+1;
        end
    end
end
if isfield(f,'mu_norm')
    f.jlc_norm=permute(f.mu_norm,[3 2 1]);
    f.jlc_norm_title='J_{||} (norm)';
    f.jlmnc_norm=zeros(f.mn,f.ns);
    i=1;
    for j=1:f.mpol+1
        for k=1:2*f.ntor+1
            f.jlmnc_norm(i,:)=f.jlc_norm(j,k,:);
            i=i+1;
        end
    end
end
if isfield(f,'muold')
    f.jlc_old=permute(f.muold,[3 2 1]);
    f.jlc_old_title='J_{||} (old)';
    f.jlmnc_old=zeros(f.mn,f.ns);
    i=1;
    for j=1:f.mpol+1
        for k=1:2*f.ntor+1
            f.jlmnc_old(i,:)=f.jlc_old(j,k,:);
            i=i+1;
        end
    end
end
% Pressure
if isfield(f,'pm_norm')
    f.press_norm=permute(f.pm_norm,[3 2 1]);
    f.press_title='Plasma Pressure (norm)';
    f.pressmnc=zeros(f.mn,f.ns);
    i=1;
    for j=1:f.mpol+1
        for k=1:2*f.ntor+1
            f.pressmnc_norm(i,:)=f.press_norm(j,k,:);
            i=i+1;
        end
    end
end
if isfield(f,'pm')
    f.press=permute(f.pm,[3 2 1])./f.beta;
    f.press_title='Plasma Pressure';
    f.pressmnc=zeros(f.mn,f.ns);
    i=1;
    for j=1:f.mpol+1
        for k=1:2*f.ntor+1
            f.pressmnc(i,:)=f.press(j,k,:);
            i=i+1;
        end
    end
end
if isfield(f,'pmold')
    f.press_old=permute(f.pmold,[3 2 1])./f.beta;
    f.press_title='Plasma Pressure (old)';
    f.pressmnc_old=zeros(f.mn,f.ns);
    i=1;
    for j=1:f.mpol+1
        for k=1:2*f.ntor+1
            f.pressmnc_old(i,:)=f.press_old(j,k,:);
            i=i+1;
        end
    end
end
% Create the Resonance Array
f.M=1:max(f.xm);
f.N=min(f.xn(f.xn > 0)):min(f.xn(f.xn > 0)):max(f.xn);
Minv=1./f.M;
f.iota_res=f.N'*Minv;
iota_spline=pchip(1:f.ns,abs(f.iota));
f.s_res=0.0.*f.iota_res;
for i=1:size(f.s_res,1)
    for j=1:size(f.s_res,2)
        if (f.iota_res(i,j) > min(abs(f.iota))) && (f.iota_res(i,j) < max(abs(f.iota)))
            f_iota=@(x)ppval(iota_spline,x)-f.iota_res(i,j);
            f.s_res(i,j)=fzero(f_iota,50);
        end
    end
end
% Clean Arrays
for i=1:f.ns
    f.rbc(:,:,i)=f.rbc(:,:,i).*double(f.mode);
    f.zbs(:,:,i)=f.zbs(:,:,i).*double(f.mode);
end
% Others
f.mu=single(f.mda);  % Moved here to protect mu
f.nv=single(f.nda);
f.aspect=-1;
f.Aminor=-1;
f.Volume=-1;
f.betatot=-1;
f.betapol=-1;
f.betator=-1;
f.betaxis=-1;
f.VolAvgB=-1;
f.b0=-1;
f.Itor=-1;
f.IonLarmor=-1;
return
end

