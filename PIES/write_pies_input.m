function write_pies_input(pies_input,pies_pltflg,exlsta,varargin)
%WRITE_PIES_INPUT Writes the Input namelists for PIES.
%   WRITE_PIES_INPUT writes the input namelist file for PIES.  It requires
%   the user supply the INPUT, PLTFLG, and EXLSTA namelists.
%   
%   WRITE_PIES_INPUT(pies_input,pies_pltflg,exlsta,pies_data,'back')  Will
%   output the background coordinates found in the pies_data structure.
%   
%   WRITE_PIES_INPUT(pies_input,pies_pltflg,exlsta,pies_data,'bfield')
%   Will output the magnetic field found in the pies_data structure.
%   
%   WRITE_PIES_INPUT(pies_input,pies_pltflg,exlsta,pies_data,'modes')  Will
%   output the mode selection matrix found in the pies_data structure.
%   
%   WRITE_PIES_INPUT(pies_input,pies_pltflg,exlsta,pies_data,'write_all')
%   Will output the background coordinates, magnetic field and mode
%   selection matrix found in the pies_data structure.  Note that profile
%   information must be manually copied into this file.
%   
%   WRITE_PIES_INPUT(pies_input,pies_pltflg,exlsta,'restart')  Will cause
%   the 'restart' flag to be utilized.
%   
%   The user may also calculate a new grid.  This is done by passing
%   'mpol', 'ntor' and 'ns' to the function.  In fourier space the
%   matricies are simply expanded with zeros and a new mode selection
%   matrix is calculated.  When a finer radial resolution is desired, the
%   code utilizes the pchip routine (piecewise cubic hermite polynomial
%   interpolation) to calculate values at each new radial gridpoint.  
%   Interpolation is preformed along the radial coordinate for each m and n
%   in the pies_data strcuture.

% Set some defaults
offset=3;
i=1;
restart=0;
piesdata=[];
vmecdata=[];
b1=[];
filename='pies.in';
writex=0;
writeb=0;
writemodes=0;
mpol=[];
ntor=[];
ns=[];
mu0=4*pi*1.e-7;
lfreeb=1;
prof_file=[];


% Handle varargin
while i <= nargin-offset
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch varargin{i}.datatype
                case 'pies_out'
                    piesdata=varargin{i};
                case 'wout'
                    vmecdata=varargin{i};
            end
        end
    elseif ischar(varargin{i})
        switch varargin{i}
            case 'filename'
                i=i+1;
                filename=varargin{i};
            case 'profiles'
                i=i+1;
                prof_file=varargin{i};
            case 'back'
                writex=1;
            case 'bfield'
                writeb=1;
            case 'modes'
                writemodes=1;
            case 'write_all'
                writemodes=1;
                writeb=1;
                writex=1;
            case 'restart'
                restart=1;
            case 'mpol'
                i=i+1;
                mpol=varargin{i};
            case 'ntor'
                i=i+1;
                ntor=varargin{i};
            case 'ns'
                i=i+1;
                ns=varargin{i};
            case 'fixed'
                lfreeb=0;
        end
    end
    i=i+1;
end

% Merge the namelist so everything is defined
pies_input=struct_merge(pies_input,pies_namelist_init('INPUT'));
pies_pltflg=struct_merge(pies_pltflg,pies_namelist_init('PLTFLG'));
exlsta=struct_merge(exlsta,pies_namelist_init('EXLSTA'));

% Check to see if piesdata or vmecdata passed if writing out stuff
if (writex || writeb || writemodes) && (isempty(piesdata) && isempty(vmecdata))
    disp('ERROR:  No PIES or VMEC data structure provided');
    return
elseif (writex || writeb || writemodes) && (isempty(vmecdata))
    x1=piesdata.rbc;
    mid=(size(x1,2)+1)/2;           % Find m=0 components
    x1(1,mid,:)=x1(1,mid,:)-piesdata.rmaj; % Center at x=0
    x2=piesdata.zbs;
    b1=piesdata.bsc;
    b2=piesdata.buc;
    b3=piesdata.bvc;
    mode=piesdata.mode;
    pies_input.k=piesdata.k;
    pies_input.m=piesdata.m;
    pies_input.n=piesdata.n;
    pies_input.mda=piesdata.mda;
    pies_input.nda=piesdata.nda;
    pies_input.vmecf=1;
    exlsta.vmecb=1;
    exlsta.read_brho=1;
    exlsta.vmec_ignore_sym=1;
elseif (writex || writeb || writemodes) && (isempty(piesdata))
    % Remove center of outer boundary
    theta=0:2*pi/359:2*pi;
    zeta=0;
    r=cfunct(theta,zeta,vmecdata.rmnc,vmecdata.xm,vmecdata.xn);
    center=max(max(max(r)))-1.0;
    vmecdata.rbc(1,vmecdata.ntor+1,:) = vmecdata.rbc(1,vmecdata.ntor+1,:) - center;
    pies_input.rmaj= center;
    % Get the spline profiles
    [p_spline j_spline pies_input.beta pies_input.betai pies_input.lp pies_input.lj pies_input.kspln pies_input.iote] = pies_profiles(vmecdata);
    pies_input.ispln=int32(2);
    exlsta.flux_quadrature = int32(1);
    exlsta.use_poly_for_currp_press=int32(0);       % Tell PIES not to look for VMEC AC AM
    exlsta.use_spline_der_for_press=int32(1);
    % We need to expand the X arrays to incoporate the mode padding
    x1=vmecdata.rbc;
    x2=vmecdata.zbs;
    b2=vmecdata.buc;
    b3=vmecdata.bvc;
    pies_input.k=vmecdata.ns;
    pies_input.m=(vmecdata.mpol-1)*2;
    pies_input.n=vmecdata.ntor*2;
    pies_input.mda=(vmecdata.mpol-1)*2;
    pies_input.nda=vmecdata.ntor*2;
    pies_input.vmecf=1;
    if max(vmecdata.iotaf) < 0.0
        disp('ERROR:  negative iota');
        return
    end
    pies_input.iotamx = 1.2*max(vmecdata.iotaf);
    pies_input.iotamn = 0.8*min(vmecdata.iotaf);
    pies_input.beta = mu0;
    pies_input.bsubpi = vmecdata.RBtor;
    pies_input.nper = vmecdata.nfp;
    exlsta.vmecb=1;
    exlsta.read_brho=0;
    % Pad the arrays
    x1_temp=zeros(pies_input.mda+1,2*pies_input.nda+1,vmecdata.ns);
    x2_temp=zeros(pies_input.mda+1,2*pies_input.nda+1,vmecdata.ns);
    b1_temp=zeros(pies_input.mda+1,2*pies_input.nda+1,vmecdata.ns);
    b2_temp=zeros(pies_input.mda+1,2*pies_input.nda+1,vmecdata.ns);
    b3_temp=zeros(pies_input.mda+1,2*pies_input.nda+1,vmecdata.ns);
    mode=zeros(pies_input.mda+1,2*pies_input.nda+1);
    offset = pies_input.nda+1;
    for m=0:vmecdata.mpol-1
        for n = -vmecdata.ntor:vmecdata.ntor
            x1_temp(m+1,n+offset,:)=x1(m+1,n+vmecdata.ntor+1,:);
            x2_temp(m+1,n+offset,:)=x2(m+1,n+vmecdata.ntor+1,:);
            b2_temp(m+1,n+offset,:)=b2(m+1,n+vmecdata.ntor+1,:);
            b3_temp(m+1,n+offset,:)=b3(m+1,n+vmecdata.ntor+1,:);
            mode(m+1,n+offset) = 1;
        end
    end
    x1=x1_temp;
    x2=x2_temp;
    b1=b1_temp;
    b2=b2_temp;
    b3=b3_temp;
    if isempty(ns)
        ns = vmecdata.ns;  % Do this to we spline over rho in the next section
    end
    % We need to go from (mu+nv) to (nv-mu)
    % First we can go from (mu+nv) to (mu-nv) by inverting the n indexes
    % Then we use -(mu-nv) = (nv-mu) so we negate the odd functions
    x1=flipdim(x1,2);
    x2=-flipdim(x2,2);
    b1=-flipdim(b1,2);
    b2=flipdim(b2,2);
    b3=flipdim(b3,2);
end

% Check and interpolate to the new PIES radial grid
if (~isempty(ns) && (~isempty(piesdata) || ~isempty(vmecdata)))
    disp('  - Interpolating to new radial grid');
    if (~isempty(piesdata))
        rho_old=0:1/double(pies_input.k):1;
    elseif (~isempty(vmecdata))
        rho_old=vmecdata.phi./vmecdata.phi(vmecdata.ns);
        rho_old=sqrt(rho_old);
    end
    rho_new=0:1/double(ns):1;
    mmax=size(x1,1);
    nmax=size(x1,2);
    x1_new=zeros(mmax,nmax,ns+1);
    x2_new=zeros(mmax,nmax,ns+1);
    b1_new=zeros(mmax,nmax,ns+1);
    b2_new=zeros(mmax,nmax,ns+1);
    b3_new=zeros(mmax,nmax,ns+1);
    for m=1:mmax
        for n=1:nmax
            if mode(m,n)
                x1_new(m,n,:)=pchip(rho_old,x1(m,n,:),rho_new);
                x2_new(m,n,:)=pchip(rho_old,x2(m,n,:),rho_new);
                b1_new(m,n,:)=pchip(rho_old,b1(m,n,:),rho_new);
                b2_new(m,n,:)=pchip(rho_old,b2(m,n,:),rho_new);
                b3_new(m,n,:)=pchip(rho_old,b3(m,n,:),rho_new);
            end
        end
    end
    x1=x1_new;
    x2=x2_new;
    b1=b1_new;
    b2=b2_new;
    b3=b3_new;
    % HANDLE VMEC DATA
    if (~isempty(vmecdata))
        b1 = [];
    end
    % Now autoset some things
    pies_input.lpinch=ns-1;
    if exlsta.mu_stochf > 0
        pies_input.islres=ns-2;
        exlsta.kstoch=round(ns*exlsta.kstoch/pies_input.k);
        exlsta.islre2=ns-1;
    else
        exlsta.kstoch=0;
        pies_input.islres=round(ns*pies_input.islres/pies_input.k);
        exlsta.islre2=round(ns*exlsta.islre2/pies_input.k);
    end
    pies_input.k=ns;
end

% Check and expand the PIES output arrays
if (~isempty(mpol) || ~isempty(ntor)) && (~isempty(piesdata) || ~isempty(vmecdata))
    % Handle Poloidal Expansion
    if isempty(mpol), mpol=pies_input.m; end
    if isempty(ntor), ntor=pies_input.n; end
    if mpol > pies_input.m
        disp('  - Expanding poloidal modes');
        m1=size(x1,1);
        pies_input.m=mpol;
        pies_input.mda=mpol;
        x1(m1+1:pies_input.m+1,:,:)=0.0;
        x2(m1+1:pies_input.m+1,:,:)=0.0;
        b1(m1+1:pies_input.m+1,:,:)=0.0;
        b2(m1+1:pies_input.m+1,:,:)=0.0;
        b3(m1+1:pies_input.m+1,:,:)=0.0;
        mode(m1+1:pies_input.m+1,:)=0.0;
        m1=size(x1,1);
        nstart=round(pies_input.n/2)+1;
        nend=nstart+pies_input.n;
        mode(1:round(m1/2),nstart:nend)=1;
    elseif (mpol ~= pies_input.m)
        disp(['ERROR: mpol = ' num2str(mpol) ', but m = ' num2str(pies_input.m)]);
        return
    end
    % Handle Toroidal Expansion
    if ntor > pies_input.n
        disp('  - Expanding toroidal modes');
        k1=size(x1,3);
        n1=(size(x1,2)-1)/2;
        m1=size(x1,1);
        temp_arr=zeros(m1,ntor-n1,k1);
        x1=cat(2,temp_arr,x1,temp_arr);
        x2=cat(2,temp_arr,x2,temp_arr);
        b1=cat(2,temp_arr,b1,temp_arr);
        b2=cat(2,temp_arr,b2,temp_arr);
        b3=cat(2,temp_arr,b3,temp_arr);
        temp_arr=zeros(m1,ntor-n1);
        mode=cat(2,temp_arr,mode,temp_arr);
        nstart=(ntor/2)+1;
        nend=nstart+ntor;
        mode(1:round(m1/2),nstart:nend)=1;
        pies_input.n=ntor;
        pies_input.nda=ntor;
    elseif (ntor ~= pies_input.n)
        disp(['ERROR: mpol = ' num2str(ntor) ', but m = ' num2str(pies_input.n)]);
        return
    end
end

% Free boundary settings
if lfreeb
    pies_input.freeb=1;
    pies_input.lpinch=ns-1;
    exlsta.calculate_bvac_once=1;
    exlsta.storvac=1;
    exlsta.virtual_casing=1;
    exlsta.remove_current_in_vacuum_region=2;
    exlsta.setbc_override2=0;
else
    pies_input.freeb=0;
    pies_input.lpinch=0;
    exlsta.calculate_bvac_once=0;
    exlsta.storvac=0;
    exlsta.virtual_casing=0;
    exlsta.remove_current_in_vacuum_region=2;
    exlsta.setbc_override=1;
    exlsta.setbc_override2=1;
end

% Pressure profiles

    

% Open File
fid=fopen(filename,'w+');

% Restart Flag
if restart, fprintf(fid,'''restart''\n'); else fprintf(fid,'''begin''\n'); end

%- INPUT Namelist
fprintf(fid,'&INPUT\n');
fprintf(fid,'  ITER2  = %-4d\n',pies_input.iter2);
fprintf(fid,'  CONVG  = %-9.2e\n',pies_input.convg);
fprintf(fid,'  NUMLST = %-1d\n',pies_input.numlst);
fprintf(fid,'!-----  Grid Control  -----\n');
fprintf(fid,'  K      = %-3d\n',pies_input.k);
fprintf(fid,'  M      = %-3d\n',pies_input.m);
fprintf(fid,'  N      = %-3d\n',pies_input.n);
fprintf(fid,'  MDA    = %-3d\n',pies_input.mda);
fprintf(fid,'  NDA    = %-3d\n',pies_input.nda);
fprintf(fid,'  MDSLCT = %-3d\n',pies_input.mdslct);
fprintf(fid,'  LPINCH = %-4d\n',pies_input.lpinch);
fprintf(fid,'!-----  Field Line Following  -----\n');
fprintf(fid,'  NFOLMX = %-7d\n',pies_input.nfolmx);
fprintf(fid,'  FTPREC = %-9.2e\n',pies_input.ftprec);
fprintf(fid,'  FTFOL  = %-9.2e\n',pies_input.ftfol);
fprintf(fid,'  LININT = %-7d\n',pies_input.linint);
fprintf(fid,'  LINTOL = %-9.2e\n',pies_input.lintol);
fprintf(fid,'  DKLIM  = %-9.2e\n',pies_input.dklim);
fprintf(fid,'  DEVPAR = %-5.3f\n',pies_input.devpar);
fprintf(fid,'!-----  Configuration  -----\n');
fprintf(fid,'  NPER   = %-5.3f\n',pies_input.nper);
fprintf(fid,'  RMAJ   = %-20.10e\n',pies_input.rmaj);
if pies_input.cyl, fprintf(fid,'  CYL    = T\n'); else fprintf(fid,'  CYL    = F\n'); end
if pies_input.analyt, fprintf(fid,'  ANALYT = T\n'); else fprintf(fid,'  ANALYT = F\n'); end
fprintf(fid,'  FREEB  = %-1d\n',pies_input.freeb);
fprintf(fid,'  VMECF  = %-1d\n',pies_input.vmecf);
fprintf(fid,'  BSUBPI = %-20.10e\n',pies_input.bsubpi);
fprintf(fid,'  SETBC  = %-1d\n',pies_input.setbc);
fprintf(fid,'  ISLRES = %-5.3f\n',pies_input.islres);
fprintf(fid,'!-----  Rotational Transform Bounds  -----\n');
fprintf(fid,'  IOTAMX = %-12.10f\n',pies_input.iotamx);
fprintf(fid,'  IOTAMN = %-12.10f\n',pies_input.iotamn);
fprintf(fid,'!-----  Profile Parameters  -----\n');
fprintf(fid,'  ISPLN  = %-1d\n',pies_input.ispln);
fprintf(fid,'  LP     = %-4d\n',pies_input.lp);
fprintf(fid,'  BETA   = %-20.10e\n',pies_input.beta);
fprintf(fid,'  LJ     = %-4d\n',pies_input.lj);
fprintf(fid,'  BETAI  = %-20.10e\n',pies_input.betai);
fprintf(fid,'  RTOKFK = %-20.10e\n',pies_input.rtokfk);
fprintf(fid,'  ADJST  = %-1d\n',pies_input.adjst);
fprintf(fid,'  IOTE   = %12.10f\n',pies_input.iote);
fprintf(fid,'!-----  Other Flags  -----\n');
if pies_input.poinc, fprintf(fid,'  POINC    = T\n'); else fprintf(fid,'  POINC    = F\n'); end
fprintf(fid,'  POSTP   = %1d\n',pies_input.postp);
fprintf(fid,'  FRCFRE  = %1d\n',pies_input.frcfre);
fprintf(fid,'/\n');

%- PLTFLG Namelist
fprintf(fid,'&PLTFLG\n');
fprintf(fid,'!-----  General Plotting Switches  -----\n');
fprintf(fid,'  PLTSF                 = %-1d\n',pies_pltflg.pltsf);
fprintf(fid,'  PLTALF                = %-1d\n',pies_pltflg.pltalf);
fprintf(fid,'  IOTAF                 = %-1d\n',pies_pltflg.iotaf);
fprintf(fid,'  QF                    = %-1d\n',pies_pltflg.qf);
fprintf(fid,'  DELTAF                = %-1d\n',pies_pltflg.deltaf);
fprintf(fid,'  DNTUP1F               = %-1d\n',pies_pltflg.dntup1f);
fprintf(fid,'  DNTUP2F               = %-1d\n',pies_pltflg.dntup2f);
fprintf(fid,'  DNTUP3F               = %-1d\n',pies_pltflg.dntup3f);
fprintf(fid,'  RESDLF                = %-1d\n',pies_pltflg.resdlf);
fprintf(fid,'  EDGEFL                = %-1d\n',pies_pltflg.edgefl);
fprintf(fid,'  EDGF1                 = %-1d\n',pies_pltflg.edgf1);
fprintf(fid,'  EDGF2                 = %-1d\n',pies_pltflg.edgf2);
fprintf(fid,'  EDGF3                 = %-1d\n',pies_pltflg.edgf3);
fprintf(fid,'  ISLPLTF               = %-1d\n',pies_pltflg.islpltf);
fprintf(fid,'  FREEBP                = %-1d\n',pies_pltflg.freebp);
fprintf(fid,'!-----  Poincare  -----\n');
fprintf(fid,'  POINCAF               = %-1d\n',pies_pltflg.poincaf);
fprintf(fid,'  POINCF                = %-1d\n',pies_pltflg.poincf);
fprintf(fid,'  POINCM1               = %-1d\n',pies_pltflg.poincm1);
fprintf(fid,'  POINCM2               = %-1d\n',pies_pltflg.poincm2);
fprintf(fid,'  POINCMG               = %-1d\n',pies_pltflg.poincmg);
fprintf(fid,'  RPOINCF               = %-1d\n',pies_pltflg.rpoincf);
fprintf(fid,'  RPOINC_PLOT_RHO       = %-1d\n',pies_pltflg.rpoinc_plot_rho);
fprintf(fid,'  HUDSON_EDGES_PLT      = %-1d\n',pies_pltflg.hudson_edges_plt);
fprintf(fid,'  WRITE_POINCARE_COORDINATES = %-1d\n',pies_pltflg.write_poincare_coordinates);
fprintf(fid,'!-----  Coordinate/Jacobian Plotting  -----\n');
fprintf(fid,'  RHOMAGF               = %-1d\n',pies_pltflg.rhomagf);
fprintf(fid,'  XPLTF                 = %-1d\n',pies_pltflg.xpltf);
fprintf(fid,'  XIJF                  = %-1d\n',pies_pltflg.xijf);
fprintf(fid,'  UFXF                  = %-1d\n',pies_pltflg.ufxf);
fprintf(fid,'  XMAGF                 = %-1d\n',pies_pltflg.xmagf);
fprintf(fid,'  XMAGPEF               = %-1d\n',pies_pltflg.xmagpef);
fprintf(fid,'  MODAMXF               = %-1d\n',pies_pltflg.modamxf);
fprintf(fid,'  RIJF                  = %-1d\n',pies_pltflg.rijf);
fprintf(fid,'  DXF                   = %-1d\n',pies_pltflg.dxf);
fprintf(fid,'  RHOJAF                = %-1d\n',pies_pltflg.rhojaf);
fprintf(fid,'  BJACF                 = %-1d\n',pies_pltflg.bjacf);
fprintf(fid,'  DPSDNF                = %-1d\n',pies_pltflg.dpsdnf);
fprintf(fid,'  DPSDNIF               = %-1d\n',pies_pltflg.dpsdnif);
fprintf(fid,'  IRREG_GRID_PLOT       = %-1d\n',pies_pltflg.irreg_grid_plot);
fprintf(fid,'!-----  Surface Plotting  -----\n');
fprintf(fid,'  VESSELF               = %-1d\n',pies_pltflg.vesself);
fprintf(fid,'  BACKF                 = %-1d\n',pies_pltflg.backf);
fprintf(fid,'  BGNDF                 = %-1d\n',pies_pltflg.bgndf);
fprintf(fid,'  MAGGNF                = %-1d\n',pies_pltflg.maggnf);
fprintf(fid,'  MAGGNAF               = %-1d\n',pies_pltflg.maggnaf);
fprintf(fid,'  MAGSNFF               = %-1d\n',pies_pltflg.magsnff);
fprintf(fid,'  RHOMAGF               = %-1d\n',pies_pltflg.rhomagf);
fprintf(fid,'!-----  Pressure  -----\n');
fprintf(fid,'  DPDPSI_PLOT           = %-1d\n',pies_pltflg.dpdpsi_plot);
fprintf(fid,'  PMNF                  = %-1d\n',pies_pltflg.pmnf);
fprintf(fid,'  PMMNF                 = %-1d\n',pies_pltflg.pmmnf);
fprintf(fid,'  PMIJF                 = %-1d\n',pies_pltflg.pmijf);
fprintf(fid,'  PMIJBF                = %-1d\n',pies_pltflg.pmijbf);
fprintf(fid,'  PRESSURE_CONTOUR_PLT  = %-1d\n',pies_pltflg.pressure_contour_plt);
fprintf(fid,'  DP1F                  = %-1d\n',pies_pltflg.dp1f);
fprintf(fid,'  DP23F                 = %-1d\n',pies_pltflg.dp23f);
fprintf(fid,'!-----  Current  -----\n');
fprintf(fid,'  MUMNF                 = %-1d\n',pies_pltflg.mumnf);
fprintf(fid,'  MUIJF                 = %-1d\n',pies_pltflg.muijf);
fprintf(fid,'  MUIJBF                = %-1d\n',pies_pltflg.muijbf);
fprintf(fid,'  JJUPF                 = %-1d\n',pies_pltflg.jjupf);
fprintf(fid,'  JPSJUPF               = %-1d\n',pies_pltflg.jpsjupf);
fprintf(fid,'  JPF                   = %-1d\n',pies_pltflg.jpf);
fprintf(fid,'  JJUPMF                = %-1d\n',pies_pltflg.jjupmf);
fprintf(fid,'  JJUPIJF               = %-1d\n',pies_pltflg.jjupijf);
fprintf(fid,'!-----  Magnetic Field  -----\n');
fprintf(fid,'  VMECBFP               = %-1d\n',pies_pltflg.vmecbfp);
fprintf(fid,'  BUPFL                 = %-1d\n',pies_pltflg.bupfl);
fprintf(fid,'  BPHIF                 = %-1d\n',pies_pltflg.bphif);
fprintf(fid,'  BPHIBF                = %-1d\n',pies_pltflg.bphibf);
fprintf(fid,'  MODBPF                = %-1d\n',pies_pltflg.modbpf);
fprintf(fid,'  BXBYFL                = %-1d\n',pies_pltflg.bxbyfl);
fprintf(fid,'  UBXBYF                = %-1d\n',pies_pltflg.ubxbyf);
fprintf(fid,'/\n');


%- EXLSTA Namelist
fprintf(fid,'&EXLSTA\n');
fprintf(fid,'!-----  General Options  -----\n');
fprintf(fid,'  BLEND_B                  = %-5.3f\n',exlsta.blend_b);
fprintf(fid,'  UMINV                    = %-1d\n',exlsta.uminv);
fprintf(fid,'  NSAV                     = %-4d\n',exlsta.nsav);
fprintf(fid,'  USE_VACFLD_KMAG          = %-1d\n',exlsta.use_vacfld_kmag);
fprintf(fid,'  CALCULATE_BVAC_ONCE      = %-1d\n',exlsta.calculate_bvac_once);
fprintf(fid,'  STORVAC                  = %-1d\n',exlsta.storvac);
fprintf(fid,'  LOCAL_J                  = %-1d\n',exlsta.local_j);
fprintf(fid,'  VIRTUAL_CASING           = %-1d\n',exlsta.virtual_casing);
fprintf(fid,'  IFTMTH                   = %-1d\n',exlsta.iftmth);
fprintf(fid,'  FBCF                     = %-1d\n',exlsta.fbcf);
fprintf(fid,'  MBDSF                    = %-1d\n',exlsta.mbdsf);
fprintf(fid,'  KERNBICHLER_WRITE        = %-1d\n',exlsta.kernbichler_write);
fprintf(fid,'  WRITE_EDGE_DATA          = %-1d\n',exlsta.write_edge_data);
fprintf(fid,'  NTOLIT                   = %-1d\n',exlsta.ntolit);
fprintf(fid,'  ITERZ2                   = %-1d\n',exlsta.iterz2);
fprintf(fid,'!-----  Field Line Following  -----\n');
fprintf(fid,'  NAXTOL                   = %-1d\n',exlsta.naxtol);
fprintf(fid,'  LINTOLF                  = %-1d\n',exlsta.lintolf);
fprintf(fid,'  USE_LSODE                = %-1d\n',exlsta.use_lsode);
fprintf(fid,'!-----  VMEC Related  -----\n');
fprintf(fid,'  VMECBF                   = %-1d\n',exlsta.vmecbf);
fprintf(fid,'  READ_BRHO                = %-1d\n',exlsta.read_brho);
fprintf(fid,'  USE_POLY_FOR_CURRP_PRESS = %-1d\n',exlsta.use_poly_for_currp_press);
fprintf(fid,'  VMEC_IGNORE_SYM          = %-1d\n',exlsta.vmec_ignore_sym);
fprintf(fid,'  BLOAT                    = %-20.10e\n',exlsta.bloat);
fprintf(fid,'  HIRSHF                   = %-1d\n',exlsta.hirshf);
fprintf(fid,'  DEV_VMEC_F               = %-1d\n',exlsta.dev_vmec_f);
fprintf(fid,'  REMOVE_CURRENT_IN_VACUUM_REGION = %-1d\n',exlsta.remove_current_in_vacuum_region);
fprintf(fid,'  SETBC_OVERRIDE           = %-1d\n',exlsta.setbc_override);
fprintf(fid,'!-----  STOCHASTIC Options  -----\n');
fprintf(fid,'  MU_STOCHF                = %-1d\n',exlsta.mu_stochf);
fprintf(fid,'  KSTOCH                   = %-4d\n',exlsta.kstoch);
fprintf(fid,'  GRADFL                   = %-1d\n',exlsta.gradfl);
fprintf(fid,'  ISLRE2                   = %-5.3f\n',exlsta.islre2);
fprintf(fid,'  LINEAR_INTERPOLATE_STINE_COORD  = %-1d\n',exlsta.linear_interpolate_stine_coord);
fprintf(fid,'  JDEV_CAL_FROM_DEV_IN_POLAR_COOR = %-1d\n',exlsta.jdev_cal_from_dev_in_polar_coor);
fprintf(fid,'  DEV_NORM                 = %-1d\n',exlsta.dev_norm);
fprintf(fid,'  ISMHMU                   = %-1d\n',exlsta.ismhmu);
fprintf(fid,'  ISMHMU_P                 = %-1d\n',exlsta.ismhmu_p);
fprintf(fid,'  MOREMU                   = %-1d\n',exlsta.moremu);
fprintf(fid,'  DEVRAT                   = %-2d\n',exlsta.devrat);
fprintf(fid,'  OUT_OF_PHASE             = %-1d\n',exlsta.out_of_phase);
fprintf(fid,'  ISLEDGF                  = %-1d\n',exlsta.isledgf);
fprintf(fid,'  IBISEC                   = %-4d\n',exlsta.ibisec);
fprintf(fid,'  BOOT_MODEL_F             = %-1d\n',exlsta.boot_model_f);
fprintf(fid,'  HUDSON_DIAGNOSTIC        = %-1d\n',exlsta.hudson_diagnostic);
fprintf(fid,'  HUDSON_EDGES             = %-1d\n',exlsta.hudson_edges);
fprintf(fid,'  HUDSON_EDGES_IN_PHASE    = %-1d\n',exlsta.hudson_edges_in_phase);
fprintf(fid,'  EDGDVF                   = %-1d\n',exlsta.edgdvf);
fprintf(fid,'  JKLMF                    = %-1d\n',exlsta.jklmf);
fprintf(fid,'!-----  Grid Options  -----\n');
fprintf(fid,'  IRREGULAR_GRID           = %-4d\n',exlsta.irregular_grid);
fprintf(fid,'  DKLIM_IRREG_GRID_INPUT   = %-9.2e\n',exlsta.dklim_irreg_grid_input);
fprintf(fid,'  USE_POLAR_COORD_TO_MAP   = %-1d\n',exlsta.use_polar_coord_to_map);
fprintf(fid,'  FLUX_QUADRATURE          = %-1d\n',exlsta.flux_quadrature);
fprintf(fid,'  USE_INTERPOLATED_GRID    = %-1d\n',exlsta.use_interpolated_grid);
fprintf(fid,'  IAXBIS                   = %-1d\n',exlsta.iaxbis);
fprintf(fid,'  ISMTHM                   = %-1d\n',exlsta.ismthm);
fprintf(fid,'  IRMVXM                   = %-1d\n',exlsta.irmvxm);
fprintf(fid,'  IMAPMG                   = %-1d\n',exlsta.imapmg);
fprintf(fid,'  IDAGTR                   = %-1d\n',exlsta.idagtr);
fprintf(fid,'  SPIMTH                   = %-4d\n',exlsta.spimth);
fprintf(fid,'  IORGFX                   = %-1d\n',exlsta.iorgfx);
fprintf(fid,'  ISTNXM                   = %-1d\n',exlsta.istnxm);
fprintf(fid,'  IWRTMG                   = %-1d\n',exlsta.iwrtmg);
fprintf(fid,'!-----  Spline Options  -----\n');
fprintf(fid,'  USE_EZSPLINE_INTERP               = %-1d\n',exlsta.use_ezspline_interp);
fprintf(fid,'  USE_EZSPLINE_IN_FBPH_AND_CRDINT   = %-1d\n',exlsta.use_ezspline_in_fbph_and_crdint);
fprintf(fid,'  USE_SPLINES_TO_INVERT_COORDINATES = %-1d\n',exlsta.use_splines_to_invert_coordinates);
fprintf(fid,'  B_XSI_B_ETA_TEST                  = %-1d\n',exlsta.b_xsi_b_eta_test);
fprintf(fid,'  USE_SPLINE_DER_FOR_PRESS          = %-1d\n',exlsta.use_spline_der_for_press);
fprintf(fid,'!-----  Coil Options  -----\n');
fprintf(fid,'  FAC_NESCOIL              = %-20.10e\n',exlsta.fac_nescoil);
fprintf(fid,'  NW                       = %-4d\n',exlsta.nw);
fprintf(fid,'  WRITE_BDOTN              = %-1d\n',exlsta.write_bdotn);
fprintf(fid,'  SETBC_OVERRIDE2          = %-1d\n',exlsta.setbc_override2);
fprintf(fid,'  DYNAMICAL_HEALING        = %-1d\n',exlsta.dynamical_healing);
fprintf(fid,'  DYNAMICAL_HEALING_RESTART_ITER = %-4d\n',exlsta.dynamical_healing_restart_iter);
fprintf(fid,'!-----  Error Options  -----\n');
fprintf(fid,'  ISTOP                    = %-1d\n',exlsta.istop);
fprintf(fid,'  CHECK_POINT_FLAG         = %-1d\n',exlsta.check_point_flag);
fprintf(fid,'  CHANGE_JKLIM_IF_OVERLAP_OF_MAG_COORD = %-4d\n',exlsta.change_jklim_if_overlap_of_mag_coord);
fprintf(fid,'!-----  Chebyshev Blending Options  -----\n');
fprintf(fid,'  CHEBYF                   = %-1d\n',exlsta.chebyf);
fprintf(fid,'  NCYC                     = %-1d\n',exlsta.ncyc);
fprintf(fid,'/\n');

% Write Background Coordinates
if writex
    disp('  - Writing background coordinates');
    for k=1:pies_input.k+1
        for n=1:2*pies_input.n+1
            %for m=1:pies_input.m+1
            fprintf(fid,' %20.10e ',x1(:,n,k));
            %end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
        for n=1:2*pies_input.n+1
            %for m=1:pies_input.m+1
            fprintf(fid,' %20.10e ',x2(:,n,k));
            %end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

% Write Magnetic Field
if writeb
    if ~isempty(b1)
        disp('  - Writing B_rho');
        for k=1:pies_input.k+1
            for n=1:2*pies_input.n+1
                fprintf(fid,' %20.10e ',b1(:,n,k));
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    disp('  - Writing B_theta');
    for k=1:pies_input.k+1
        for n=1:2*pies_input.n+1
            fprintf(fid,' %20.10e ',b2(:,n,k));
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    disp('  - Writing B_zeta');
    fprintf(fid,'\n');
    for k=1:pies_input.k+1
        for n=1:2*pies_input.n+1
            fprintf(fid,' %20.10e ',b3(:,n,k));
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

% Write Mode selection
if writemodes
    disp('  - Writing Mode Selection');
    for n=1:2*pies_input.n+1
        fprintf(fid,' %1d ',mode(:,n));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

% Write profile file
if ~isempty(prof_file)
    disp(['  - Writing Profiles from File: ' strtrim(prof_file)]);
    fprintf(fid,'\n');
    pid=fopen(strtrim(prof_file));
    while ~feof(pid)
        line=fgetl(pid);
        fprintf(fid,'%s\n',line);
    end 
    fprintf(fid,'\n');
    fclose(pid);
elseif ~isempty(vmecdata)
    disp('  - Writing Spline Profiles');
    fprintf(fid,'\n');
    for i=1:p_spline.pieces+1
        fprintf(fid,'%20.10e\n',p_spline.breaks(i));
    end
    fprintf(fid,'\n');
    for i=1:j_spline.pieces+1
        fprintf(fid,'%20.10e\n',j_spline.breaks(i));
    end
    fprintf(fid,'\n');
    for i=1:p_spline.pieces
        for j=1:p_spline.order 
            fprintf(fid,' %20.10e ',p_spline.coefs(i,j));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,' %20.10e  %20.10e  %20.10e  %20.10e\n',0,0,0,0);  % Set last value to zero
    fprintf(fid,'\n');
    for i=1:j_spline.pieces
        for j=1:j_spline.order
            fprintf(fid,' %20.10e ',j_spline.coefs(i,j));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,' %20.10e  %20.10e  %20.10e  %20.10e\n',0,0,0,0);  % Set last value to zero
end

%- Close File
fclose(fid);

return

end

