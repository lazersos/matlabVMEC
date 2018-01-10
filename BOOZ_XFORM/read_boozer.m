function booz_data = read_boozer(filename,varargin)
%READ_BOOZER Reads a boozer_xform output
%   The READ_BOOZER function reads the output of a the booz_xform code.  It
%   currently only reads the binary file produced by version 1.0.
%   The code returns a structure with the following fields:
%
%   data:
%       nfp:        Number of field periods
%       ns:         Number of surfaces
%       aspect:     Aspect ratio.
%       rmax:       Maximum radial extent
%       rmin:       Minimum radial extent
%       betaxis:    Plasma beta on axis
%       iota:       Rotational transform profile
%       pres:       Pressure profile
%       beta:       Beta profile
%       phip:       Poloidal flux profile
%       phi:        Toroidal flux profile
%       bvco:       Toroidal field profile
%       buco:       Polodial field profile
%       idx:        Index of transformed surfaces
%       m:          Number of boozer poloidal modes
%       n:          Number of boozer toroidal modes (per field period)
%       version:    Version of booz_xform code
%       bc:         ModB Fourier modes
%       rbc:        R Fourier modes
%       zbs:        Z Fourier modes
%       ps:         p transform Fourier modes (phi transform)
%       gs:         g Fourier modes
%       datatype:   Identifies the structure type. 'boozer'
%
%   Example:
%       booz_data=read_boozer('boozmn.test');
%
%
%   See also.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           4/06/11

% Set some defaults
numdefarg=1;
int='int';
float='double';
datatype='boozer';

% Handle varargin
if nargin > numdefarg
    i=1;
    while i <= (nargin-numdefarg)
        switch varargin{i}
            otherwise
        end
        i=i+1;
    end
end

% Open the file
filename=strtrim(filename);
nname=length(filename);

if strcmp('.nc',filename(nname-2:nname))
    try
        booz_data=read_netcdf(filename,'flipdim');
        booz_data.nfp=booz_data.nfp_b;
        booz_data.ns=booz_data.ns_b;
        booz_data.aspect=booz_data.aspect_b;
        booz_data.rmax=booz_data.rmax_b;
        booz_data.rmin=booz_data.rmin_b;
        booz_data.betaxis=booz_data.betaxis_b;
        booz_data.iota=booz_data.iota_b;
        booz_data.pres=booz_data.pres_b;
        booz_data.phip=booz_data.phip_b;
        booz_data.phi=booz_data.phi_b;
        booz_data.bvco=booz_data.bvco_b;
        booz_data.buco=booz_data.buco_b;
        booz_data.m=booz_data.mboz_b;
        booz_data.n=booz_data.nboz_b;
        booz_data.mn=booz_data.mnboz_b;
        booz_data.idx=zeros(1,booz_data.ns);
        booz_data.idx(booz_data.jlist)=1;
        % Version is correct
        % Now allocate the 2D arrays
        booz_data.bmnc=zeros(booz_data.mn,booz_data.ns);
        booz_data.rmnc=zeros(booz_data.mn,booz_data.ns);
        booz_data.zmns=zeros(booz_data.mn,booz_data.ns);
        booz_data.pmns=zeros(booz_data.mn,booz_data.ns);
        booz_data.gmnc=zeros(booz_data.mn,booz_data.ns);
        for i=1:booz_data.pack_rad
            booz_data.bmnc(:,booz_data.jlist(i))=booz_data.bmnc_b(:,i);
            booz_data.rmnc(:,booz_data.jlist(i))=booz_data.rmnc_b(:,i);
            booz_data.zmns(:,booz_data.jlist(i))=booz_data.zmns_b(:,i);
            booz_data.pmns(:,booz_data.jlist(i))=booz_data.pmns_b(:,i);
            booz_data.gmnc(:,booz_data.jlist(i))=booz_data.gmn_b(:,i);       % Note it's called gmn not gmnc
        end
        % Handle the Asymetric File
        booz_data.lasym=0;
        if booz_data.lasym__logical__
            booz_data.lasym=1;
            booz_data.bmns=zeros(booz_data.mn,booz_data.ns);
            booz_data.rmns=zeros(booz_data.mn,booz_data.ns);
            booz_data.zmnc=zeros(booz_data.mn,booz_data.ns);
            booz_data.pmnc=zeros(booz_data.mn,booz_data.ns);
            booz_data.gmns=zeros(booz_data.mn,booz_data.ns);
            for i=1:booz_data.pack_rad
                booz_data.bmns(:,booz_data.jlist(i))=booz_data.bmns_b(:,i);
                booz_data.rmns(:,booz_data.jlist(i))=booz_data.rmns_b(:,i);
                booz_data.zmnc(:,booz_data.jlist(i))=booz_data.zmnc_b(:,i);
                booz_data.pmnc(:,booz_data.jlist(i))=booz_data.pmnc_b(:,i);
                booz_data.gmns(:,booz_data.jlist(i))=booz_data.gmns_b(:,i);
            end
        end
        booz_data.xm=double(booz_data.ixm_b);
        booz_data.xn=double(booz_data.ixn_b);
        booz_data=rmfield(booz_data,'nfp_b');
        booz_data=rmfield(booz_data,'ns_b');
        booz_data=rmfield(booz_data,'aspect_b');
        booz_data=rmfield(booz_data,'rmax_b');
        booz_data=rmfield(booz_data,'rmin_b');
        booz_data=rmfield(booz_data,'betaxis_b');
        booz_data=rmfield(booz_data,'iota_b');
        booz_data=rmfield(booz_data,'pres_b');
        booz_data=rmfield(booz_data,'phip_b');
        booz_data=rmfield(booz_data,'phi_b');
        booz_data=rmfield(booz_data,'bvco_b');
        booz_data=rmfield(booz_data,'buco_b');
        booz_data=rmfield(booz_data,'mboz_b');
        booz_data=rmfield(booz_data,'nboz_b');
        booz_data=rmfield(booz_data,'mnboz_b');
        booz_data=rmfield(booz_data,'bmnc_b');
        booz_data=rmfield(booz_data,'rmnc_b');
        booz_data=rmfield(booz_data,'zmns_b');
        booz_data=rmfield(booz_data,'pmns_b');
        booz_data=rmfield(booz_data,'gmn_b');
        if booz_data.lasym__logical__
            booz_data=rmfield(booz_data,'bmns_b');
            booz_data=rmfield(booz_data,'rmns_b');
            booz_data=rmfield(booz_data,'zmnc_b');
            booz_data=rmfield(booz_data,'pmnc_b');
            booz_data=rmfield(booz_data,'gmns_b');
        end
        booz_data=rmfield(booz_data,'lasym__logical__');
        booz_data=rmfield(booz_data,'ixm_b');
        booz_data=rmfield(booz_data,'ixn_b');
        booz_data.input_extension=filename(8:nname-3);
    catch booz_error
        booz_data=booz_error;
        return
    end
else
    try
        fid=fopen(filename,'r');
        booz_data.input_extension=filename;
    catch booz_error
        booz_data=fid;
        return
    end
    % Read 1D Header
    fread(fid,1,int);
    data=fread(fid,2,int);
    booz_data.nfp=data(1);
    booz_data.ns=data(2);
    data=fread(fid,4,float);
    booz_data.aspect=data(1);
    booz_data.rmax=data(2);
    booz_data.rmin=data(3);
    booz_data.betaxis=data(4);
    % Create 1D Arrays
    booz_data.iota=zeros(1,booz_data.ns);
    booz_data.pres=zeros(1,booz_data.ns);
    booz_data.beta=zeros(1,booz_data.ns);
    booz_data.phip=zeros(1,booz_data.ns);
    booz_data.phi=zeros(1,booz_data.ns);
    booz_data.bvco=zeros(1,booz_data.ns);
    booz_data.buco=zeros(1,booz_data.ns);
    booz_data.idx=zeros(1,booz_data.ns);
    fread(fid,1,int);
    % Read 1D Vars
    for i=2:booz_data.ns
        fread(fid,1,int);
        data=fread(fid,7,float);
        booz_data.iota(i)=data(1);
        booz_data.pres(i)=data(2);
        booz_data.beta(i)=data(3);
        booz_data.phip(i)=data(4);
        booz_data.phi(i) =data(5);
        booz_data.bvco(i)=data(6);
        booz_data.buco(i)=data(7);
        fread(fid,1,int);
    end
    % Read 2D Header
    fread(fid,1,int);
    data=fread(fid,3,int);
    booz_data.m=data(1);
    booz_data.n=data(2);
    booz_data.mn=data(3);
    fread(fid,1,int);
    fread(fid,1,int);
    booz_data.version=fread(fid,38,'*char')';
    fread(fid,1,int);
    % Create 2D Arrays
    booz_data.bmnc=zeros(booz_data.mn,booz_data.ns);
    booz_data.rmnc=zeros(booz_data.mn,booz_data.ns);
    booz_data.zmns=zeros(booz_data.mn,booz_data.ns);
    booz_data.pmns=zeros(booz_data.mn,booz_data.ns);
    booz_data.gmnc=zeros(booz_data.mn,booz_data.ns);
    booz_data.xm=zeros(1,booz_data.mn);
    booz_data.xn=zeros(1,booz_data.mn);
    % Read 2D First Arrays
    fread(fid,1,int);
    nsval=fread(fid,1,int);
    booz_data.idx(nsval)=1;
    fread(fid,1,int);
    for i=1:booz_data.mn
        fread(fid,1,int);
        data=fread(fid,2,int);
        booz_data.xn(i)=data(1);
        booz_data.xm(i)=data(2);
        fread(fid,1,int);
    end
    for i=1:booz_data.mn
        fread(fid,1,int);
        data=fread(fid,5,float);
        booz_data.bmnc(i,nsval)=data(1);
        booz_data.rmnc(i,nsval)=data(2);
        booz_data.zmns(i,nsval)=data(3);
        booz_data.pmns(i,nsval)=data(4);
        booz_data.gmnc(i,nsval)=data(5);
        fread(fid,1,int);
    end
    % Read 2D Arrays (do it this way because not all surfaces may have been
    % calculated.
    fread(fid,1,int);
    nsval=fread(fid,1,int);
    while ~isempty(nsval)
        fread(fid,1,int);
        booz_data.idx(nsval)=1;
        for i=1:booz_data.mn
            fread(fid,1,int);
            data=fread(fid,5,float);
            booz_data.bmnc(i,nsval)=data(1);
            booz_data.rmnc(i,nsval)=data(2);
            booz_data.zmns(i,nsval)=data(3);
            booz_data.pmns(i,nsval)=data(4);
            booz_data.gmnc(i,nsval)=data(5);
            fread(fid,1,int);
        end
        fread(fid,1,int);
        nsval=fread(fid,1,int);
    end
    % Close the file
    fclose(fid);
    booz_data.lasym=0;
end
% Now we transform to our VMEC represenataion
booz_data.rbc=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
booz_data.zbs=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
booz_data.bc=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
booz_data.ps=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
booz_data.gc=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
for k=1:booz_data.ns
    if (booz_data.idx(k) == 1)
        for j=1:booz_data.mn
            n=booz_data.xn(j)/booz_data.nfp+booz_data.n+1;
            m=booz_data.xm(j)+1;
            booz_data.rbc(m,n,k)=booz_data.rmnc(j,k);
            booz_data.zbs(m,n,k)=booz_data.zmns(j,k);
            booz_data.bc(m,n,k)=booz_data.bmnc(j,k);
            booz_data.ps(m,n,k)=booz_data.pmns(j,k);
            booz_data.gc(m,n,k)=booz_data.gmnc(j,k);
        end
    end
end
% Handle asymmetry
if booz_data.lasym
    booz_data.rbs=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
    booz_data.zbc=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
    booz_data.bs=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
    booz_data.pc=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
    booz_data.gs=zeros(booz_data.m,2*booz_data.n+1,booz_data.ns);
    for k=1:booz_data.ns
        if (booz_data.idx(k) == 1)
            for j=1:booz_data.mn
                n=booz_data.xn(j)/booz_data.nfp+booz_data.n+1;
                m=booz_data.xm(j)+1;
                booz_data.rbs(m,n,k)=booz_data.rmns(j,k);
                booz_data.zbc(m,n,k)=booz_data.zmnc(j,k);
                booz_data.bs(m,n,k)=booz_data.bmns(j,k);
                booz_data.pc(m,n,k)=booz_data.pmnc(j,k);
                booz_data.gs(m,n,k)=booz_data.gmns(j,k);
            end
        end
    end
end
% Now clean up the structure
booz_data.ntor=double(booz_data.n);
booz_data.mpol=double(booz_data.m);
booz_data.nu=double(2.*booz_data.m+6);
booz_data.nfp=double(booz_data.nfp);
% Set the dataype
booz_data.datatype=datatype;
return
end

