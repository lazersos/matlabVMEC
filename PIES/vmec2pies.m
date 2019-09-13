function vmec2pies(varargin)
%VMEC2PIES(vmec_data,coil_data) Creates a PIES input file from VMEC.
%   EXTENDER_BOUNDARY creates an extender boundary and axis file from a
%   VMEC data structure.  The boundary is created through a bloating of the
%   VMEC boundary in real space.
%
%   Options:
%       VMEC2PIES(data,'filename','test')   Sets the filename of
%       the output boundary and axis file (otherwise taken from
%       data.input_extension)
%
%       VMEC2PIES(data,'ntor',ntor)   Sets the number of toroidal
%       modes to use for the output boundary.  Similar argument exists for
%       'mpol'.
%
%       VMEC2PIES(data,'nu',nu)   Sets the number of gridpoints to
%       use in the real space poloidal direction.  A similar argument 'nv'
%       exists for the toroidal direction.  (suggest 4*mpol)  Note that
%       'nv' is per half period.  (suggest 2*ntor)
%
%       VMEC2PIES(data,'scale',scale)   Sets the scaling factor for
%       the bloated boundary in terms of a %.
%       Default is 20%.
%
%       VMEC2PIES(data,'nocoil')   Do not output the coil data, but
%       output everything else.
%
%       VMEC2PIES(data,'spline')   Calculate the current and
%       pressure splines from the VMEC data on the PIES grid.
%
%   Example:
%       vmec_data=read_vmec('wout.test');
%       vmec2pies(vmec_data);
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           5/8/11

% Set some defaults
i=1;
nu=360;
nv=36;
numin=180;
nvmin=36;
nval=-1;
mval=-1;
ns=99;
lpinch=99;
rtokfk=1.0;
scale=20;
data=[];
coildata=[];
piesinput=[];
pltflg=[];
exlsta=[];
vmecinput=[];
mu0=4*pi*1.e-7;
filename='';
makecoil=1;
iter2=5000;
usespline=0;
run_type='begin';
surfmask=[];

% Handle varargin
while i <= nargin
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch varargin{i}.datatype
                case 'VMEC_input'
                    vmecinput=varargin{i};
                case 'wout'
                    data=varargin{i};
                case 'coil_data'
                    coildata=varargin{i};
                case 'pies_input'
                    piesinput=varargin{i};
                    piesinput=rmfield(piesinput,'datatype');
                case 'pltflg'
                    pltflg=varargin{i};
                    pltflg=rmfield(pltflg,'datatype');
                case 'exlsta'
                    exlsta=varargin{i};
                    exlsta=rmfield(exlsta,'datatype');
            end
        end
    elseif ischar(varargin{i})
        switch varargin{i}
            case 'filename'
                i=i+1;
                filename=varargin{i};
            case 'ntor'
                i=i+1;
                nval=varargin{i};
            case 'mpol'
                i=i+1;
                mval=varargin{i};
            case 'nu'
                i=i+1;
                nu=varargin{i};
            case 'nv'
                i=i+1;
                nv=varargin{i};
            case 'ns'
                i=i+1;
                ns=varargin{i};
            case 'scale'
                i=i+1;
                scale=varargin{i};
            case 'surfmask'
                i=i+1;
                surfmaks=varargin{i};
            case 'nocoil'
                makecoil=0;
            case 'spline'
                usespline=1;
        end
    end
    i=i+1;
end

if makecoil
    if isempty(coildata)
        disp('ERROR:  No Coil structure passed!');
        return
    end
end

% Get the filename
if isempty(filename)
    filename=data.input_extension;
end

% Check nv
if nv < nvmin
    disp(['   - Setting nu < ' num2str(nvmin,'%2d') ' can cause issues, setting nu = ' num2str(nvmin,'%2d')]);
    nv=nvmin;
end

% Check nv
if nu < numin
    disp(['   - Setting nu < ' num2str(numin,'%2d') ' can cause issues, setting nu = ' num2str(nvmin,'%2d')]);
    nu=numin;
end

% Set PIES Input if none passed
if isempty(piesinput)
    piesinput.iter2=int32(iter2);
    piesinput.k=49;
    piesinput.m=int32(mval);
    piesinput.n=int32(nval);
    piesinput.mda=int32(32);
    piesinput.nda=int32(24);
    piesinput.numlst=int32(2);
    piesinput.mreduc=int32(1);
    piesinput.mdslct=int32(1);
    piesinput.bblnd=0;
    piesinput.vmecf=int32(1);
    piesinput.islres=1.0;
    piesinput.lintol=1.0e-7;
    piesinput.nfolmx=int32(90000);
    piesinput.ftprec=1.0e-5;
    piesinput.ftfol=1.0e-5;
    piesinput.dklim=2.0e-3;
    piesinput.linint=int32(2^19);
    piesinput.iotamx=0.16;
    piesinput.iotamn=0.1;
    piesinput.frcfre=int32(0);
    piesinput.ispln=int32(2);
    piesinput.lp=int32(50);
    piesinput.lj=int32(50);
    piesinput.rtokfk=1.0;
    piesinput.adjst=int32(0);
    piesinput.rmaj=3.0;
    piesinput.nper=3.0;
    piesinput.cyl=boolean(0);
    piesinput.island=int32(1);
    piesinput.tokchk=int32(0);
    piesinput.rot1=0.14;
    piesinput.rot2=0.00;
    piesinput.setbc=int32(1);
    piesinput.freeb=int32(1);
    piesinput.convg=1.0e-6;
    piesinput.poinc=boolean(1);
    piesinput.postp=int32(0);
end

% Set PIES Pltflg if none passed
if isempty(pltflg)
    pltflg.dpdpsi_plot=int32(1);
    pltflg.irreg_grid_plot=int32(0);
    pltflg.rpoinc_plot_rho=int32(1);
    pltflg.pltsf=int32(1);
    pltflg.backf=int32(1);              % Background Coordinate Surfaces
    pltflg.poincaf=int32(1);            % Background Coordinates
    pltflg.rpoincf=int32(1);            % Background Coordinates
    pltflg.poincm1=int32(1);            % Poincare in cartesian space at phi=0
    pltflg.poincm2=int32(1);            % Poincare in cartesian space at phi=pi/(2*nper) & phi=pi/nper
    pltflg.rhomagf=int32(0);            % Magnetic Rho coordinates
    pltflg.xmagf=int32(0);              % Magnitde of x(r) Fourier coefficients
    pltflg.dxf=int32(0);                % grad(xs) fourier coefficients
    pltflg.rhojaf=int32(0);             % rhojac in jacobs
    pltflg.bphif=int32(1);              % bphi in jacobs
    pltflg.dpsdnf=int32(0);             
    pltflg.dpsdnif=int32(0);
    pltflg.jjupf=int32(1);
    pltflg.jpsjupf=int32(0);
    pltflg.maggnf=int32(0);
    pltflg.xpltf=int32(0);
    pltflg.maggnaf=int32(0);
    pltflg.magsnff=int32(0);
    pltflg.xmagpef=int32(0);
    pltflg.modamxf=int32(0);
    pltflg.rijf=int32(0);
    pltflg.modbpf=int32(1);
    pltflg.pmmnf=int32(1);
    pltflg.pmijf=int32(0);
    pltflg.mumnf=int32(1);
    pltflg.muijf=int32(0);
    pltflg.xijf=int32(0);
    pltflg.pmijbf=int32(0);
    pltflg.muijbf=int32(0);
    pltflg.bgndf=int32(0);
    pltflg.iotaf=int32(1);
    pltflg.bphibf=int32(0);
    pltflg.ufxf=int32(0);
    pltflg.ubxbyf=int32(0);
    pltflg.deltaf=int32(0);
    pltflg.bjacf=int32(0);
    pltflg.dntup1f=int32(0);
    pltflg.dp1f=int32(0);
    pltflg.dp23f=int32(0);
    pltflg.jpf=int32(0);
    pltflg.jjupmf=int32(0);
    pltflg.jjupijf=int32(0);
    pltflg.dntup2f=int32(0);
    pltflg.dntup3f=int32(0);
    pltflg.resdlf=int32(0);
    pltflg.modbf=int32(1);
    pltflg.bupfl=int32(1);
    pltflg.edgefl=int32(0);
    pltflg.edgf1=int32(0);
    pltflg.edgf2=int32(0);
    pltflg.edgf3=int32(0);
    pltflg.poincf=int32(1);
    pltflg.vmecbfp=int32(1);
end

% Set PIES exlsta if none passed
if isempty(exlsta)
    exlsta.lsodtol=1.e-4;
    exlsta.use_splines_to_invert_coordinates=int32(1);
    exlsta.use_polar_coord_to_map=int32(1);
    exlsta.linear_interpolate_stine_coord=int32(1);
    exlsta.jdev_cal_from_dev_in_polar_coor=int32(1);
    exlsta.dklim_irreg_grid_input=2.0E-003;
    exlsta.irregular_grid=int32(0);
    exlsta.local_j=int32(1);
    exlsta.flux_quadrature=int32(1);
    exlsta.dev_norm=int32(1);
    exlsta.use_interpolated_grid=int32(1);
    exlsta.mu_stochf=int32(0);
    exlsta.calculate_bvac_once=int32(1);
    exlsta.chebyf=int32(0);
    exlsta.ncyc=int32(2);
    exlsta.istop=int32(0);
    exlsta.ismhmu=int32(1);
    exlsta.ismthm=int32(1);
    exlsta.irmvxm=int32(0);
    exlsta.istnxm=int32(1);
    exlsta.moremu=int32(1);
    exlsta.ntolit=int32(-2);
    exlsta.devrat=-1.00000000000000;
    exlsta.imapmg=int32(1);
    exlsta.idagtr=int32(1);
    exlsta.islre2=2.00000000000000;
    exlsta.fbcf=int32(0);
    exlsta.hirshf=int32(0);
    exlsta.mbdsf=int32(0);
    exlsta.gradfl=int32(1);
    exlsta.vmecbf=int32(0);
    exlsta.iorgfx=int32(0);
    exlsta.uminv=int32(5);
    exlsta.spimth=int32(2);
    exlsta.iterz2=int32(0);
    exlsta.lintolf=int32(1);
    exlsta.naxtol=int32(-2);
    exlsta.edgdvf=int32(1);
    exlsta.jklmf=int32(1);
    exlsta.blend_B=0.99;
    exlsta.check_point_flag=int32(0);
    exlsta.out_of_phase=int32(0);
    exlsta.boot_model_f=int32(0);
    exlsta.isledgf=int32(1);
    exlsta.ibisec=int32(7);
    exlsta.iftmth=int32(3);
    exlsta.storvac=int32(1);
    exlsta.nw=int32(0);
    exlsta.virtual_casing=int32(1);
    exlsta.fac_nescoil=-1.0E-007;
    exlsta.use_lsode=int32(0);
    exlsta.dev_vmec_f=int32(0);
    exlsta.b_xsi_b_eta_test=int32(1);
    exlsta.hudson_diagnostic=int32(0);
    exlsta.hudson_edges=int32(0);
    exlsta.hudson_edges_in_phase=int32(0);
    exlsta.dynamical_healing=int32(0);
    exlsta.nsav=int32(50);
    exlsta.bloat=1.000000000000000;
    exlsta.use_ezspline_interp=int32(1);
    exlsta.remove_current_in_vacuum_region=int32(2);
    exlsta.read_brho=int32(1);
    exlsta.vmecbf=int32(1);
    exlsta.use_poly_for_currp_press=int32(1);
    exlsta.vmec_ignore_sym=int32(0);
end

% Now Check for for VMEC data and set namelist values
if isempty(data)
    disp('ERROR:  No VMEC structure passed!');
    return
else
    % We assume that if VMEC_INPUT is passed then we need to add ac and am
    % to the data structure
    if ~isempty(vmecinput)
        data.ac=vmecinput.ac;
        data.am=vmecinput.am;
        if usespline == 1
            disp('WARNING:  Splining profiles, AC and AM ignored!');
        end
    end
    % Check for assymetry
    if data.iasym
        disp('ERROR: Not an up-down symmetric VMEC run');
        return
    end
    mpol=size(data.rbc,1)-1;
    ntor=(size(data.rbc,2)-1)/2;
    if mval < 0
        mval=2*(mpol+1);
    end
    if nval < 0
        nval=2*ntor;
    end
    if ntor == 0
        nval=0;
    end
    piesinput.vmecf=int32(1);
    piesinput.k=int32(ns-1);   % Default value for using grep and sed for replacement
    piesinput.m=int32(mval);
    piesinput.n=int32(nval);
    piesinput.mda=int32(mval);
    piesinput.nda=int32(nval);
    piesinput.iotamx=1.2*max(data.iotaf);
    piesinput.iotamn=0.2*min(data.iotaf);
    piesinput.rot1=data.iotaf(1);
    piesinput.rot2=data.iotaf(data.ns)-data.iotaf(1);
    piesinput.beta=mu0;
    piesinput.bsubpi=data.RBtor;
    piesinput.rmaj=data.Rmajor;
    piesinput.nper=data.nfp;
    piesinput.lp=int32(ns+1);
    piesinput.lj=int32(ns+1);
    phiedge=data.phi(data.ns);
    iint=0.0;
    if isfield(data,'ac')
        for i=1:11
            iint=iint+data.ac(i)/double(i);
        end
    else
        for i=1:11
            iint=iint+input.ac(i)/double(i);
        end
    end
    if (iint*phiedge == 0.0)
        piesinput.betai=0.0;
    else
        piesinput.betai=-(data.Itor*mu0)/(iint*phiedge);
    end
    exlsta.vmec_ignore_sym=int32(1);  % Must be set to one for NMORPH Boundaries
    if isempty(surfmask)
        surfmask=1:data.ns+1;
    elseif (length(surfmask) > data.ns+1) || (max(surfmax) > data.ns+1)
        disp('ERROR: surfmask is too large');
        return
    end
    % Clean up nda
    if ntor == 0
        piesinput.nda=int32(2);
    end
end

if usespline
    % Redefine some stuff for PIES
    piesinput.rtokfk=1.0;                           % Tell PIES not to scale anything
    piesinput.lp=int32(ns);                         % Number of breakpoints for p
    piesinput.lj=int32(ns);                         % Number of breakpoints for j
    piesinput.kspln=int32(4);                       % Order of spline (must be 4)
    exlsta.use_poly_for_currp_press=int32(0);       % Tell PIES not to look for VMEC AC AM
    piesinput.beta=mu0;                             % Scaling factor for pressure spline
    piesinput.betai=-2*pi*mu0/data.phi(data.ns);    % Scaling factor for current spline
else
    piesinput.lp=int32(ns);
    piesinput.lj=int32(ns);
end

% Construct Mode selection Matrix
modes=zeros(mval+1,2*nval+1);
for i=1:mval+1
    m=i-1;
    if (m <= mval/2)
        for j=1:2*nval+1
            n=j-nval-1;
            if ((n >= - nval/2) && (n <= nval/2))
                modes(i,j)=1;
            end
        end
    end
end

% Now we calculate the VMEC surfaces
disp('  - Calculating Surfaces');
theta=0:2*pi/(nu-1):2*pi;
zeta=0:2*pi/(nv-1)/data.nfp:2*pi/data.nfp;
r_v=cfunct(theta,zeta,data.rbc,data.nfp);
z_v=sfunct(theta,zeta,data.zbs,data.nfp);
% Create Boundary deffinition
re=zeros(nu,nv);
ze=re;
% Now calculate re and ze from r's and z's
dr=(r_v(data.ns,:,:)-r_v(1,:,:));
dz=(z_v(data.ns,:,:)-z_v(1,:,:));
dz(logical((z_v(data.ns,:,:)<0).*(dz>0)))=-dz(logical((z_v(data.ns,:,:)<0).*(dz>0)));
dz(logical((z_v(data.ns,:,:)>0).*(dz<0)))=-dz(logical((z_v(data.ns,:,:)>0).*(dz<0)));
drds=dr;
dzds=dz;
re(:,:)=r_v(data.ns,:,:)+dr.*scale/100;
ze(:,:)=z_v(data.ns,:,:)+dz.*scale/100;
% Make a quick plot
fig=figure;
subplot(1,3,1);
i=1;
plot(r_v(data.ns,:,i),z_v(data.ns,:,i),'k');
hold on
plot(r_v(1,:,i),z_v(1,:,i),'k+');
plot(re(:,i),ze(:,i),'r');
axis square; axis equal;
temp1=xlim(gca);
temp2=ylim(gca);
text(temp1(1)+(temp1(2)-temp1(1))/20.,temp2(1)+3*(temp2(2)-temp2(1))/20.,['Scale: ' num2str(scale) '%']);
text(temp1(1)+(temp1(2)-temp1(1))/20.,temp2(1)+2*(temp2(2)-temp2(1))/20.,['nu: ' num2str(nu) '  nv:' num2str(nv)]);
text(temp1(1)+(temp1(2)-temp1(1))/20.,temp2(1)+1*(temp2(2)-temp2(1))/20.,['ns: ' num2str(ns)  '  mpol: ' num2str(mval) '  ntor: ' num2str(nval)]);
hold off
subplot(1,3,2);
i=nv/4;
plot(r_v(data.ns,:,i),z_v(data.ns,:,i),'k');
hold on
plot(r_v(1,:,i),z_v(1,:,i),'k+');
plot(re(:,i),ze(:,i),'r');
axis square; axis equal;
hold off
title('New Surface');
subplot(1,3,3);
i=nv/2;
plot(r_v(data.ns,:,i),z_v(data.ns,:,i),'k');
hold on
plot(r_v(1,:,i),z_v(1,:,i),'k+');
plot(re(:,i),ze(:,i),'r');
axis square; axis equal;
hold off
legend('Boundary (VMEC)','Mag. Axis (VMEC)','Extended Boundary (PIES)');
saveas(gcf,['boundary_' strtrim(filename) '.fig']);
pause(1.5);
close(fig);
% Now add the new surface to the vmec surface
r_v(data.ns+1,:,:)=re;
z_v(data.ns+1,:,:)=ze;
% Create radial grids
r0=data.rbc(1,ntor+1,1);
rvmec=cfunct(0,0,data.rbc,data.nfp);  % VMEC radial locations
rv=rvmec(data.ns);                    % Radial extent of VMEC surfaces
rp=r_v(data.ns+1,1,1);                % PIES boundary
rpies=r0:(rp-r0)/(ns-1):rp;           % PIES radial locations
lpinch=1;
while rpies(lpinch) <= rv
    lpinch=lpinch+1;
end
if ~usespline
    piesinput.lpinch=int32(ns-2);
    rtokfk=(rv-r0)/(rp-r0);
    piesinput.rtokfk=rtokfk;
else
    piesinput.lpinch=int32(ns-2);
    rtokfk=(rv-r0)/(rp-r0);
    piesinput.rtokfk=1.0;
end
fig=figure;
plot(rpies,rpies.*0,'or');
hold on
plot([rv rv],[-1 1],'k');
plot(rvmec,rvmec.*0,'.k');
plot([rpies(piesinput.lpinch+1) rpies(piesinput.lpinch+1)],[0.0 0.0],'^');
text(rpies(1),-0.8,['RTOKFK = ' num2str(piesinput.rtokfk)]);
text(rpies(1),-0.9,['LPINCH = ',num2str(piesinput.lpinch)]);
hold off
title('New Grid');
xlabel('Radial Coordinate [m]');
legend('PIES GRID','VMEC Surface','VMEC GRID','LPINCH');
saveas(gcf,['grid_' strtrim(filename) '.fig']);
pause(1.5);
close(fig);
disp('      Interpolating in Real Space');
% Now calculate rho
rho=sqrt((r_v-r_v(1,1,1)).*(r_v-r_v(1,1,1))+z_v.*z_v);
% For each pies surface we interpolate over real space in rho
r_pies=zeros(ns,nu,nv);
z_pies=zeros(ns,nu,nv);
rho_pies=zeros(ns,nu,nv);
for i=1:nu
    for j=1:nv
        rho_pies(:,i,j)=0:rho(data.ns,i,j)/(ns-1):rho(data.ns,i,j);
    end
end
for i=1:nu
    for j=1:nv
        r_pies(:,i,j)=interp1(rho(surfmask,i,j),r_v(surfmask,i,j),rho_pies(:,i,j),'pchip');
        z_pies(:,i,j)=interp1(rho(surfmask,i,j),z_v(surfmask,i,j),rho_pies(:,i,j),'pchip');
    end
end

% Now Transform to Fourier space
refou=zeros(ns,mpol+1,2*ntor+1);
zefou=zeros(ns,mpol+1,2*ntor+1);
cosu=zeros(nu,mpol+1);
cosv=zeros(nv,2*ntor+1);
sinu=zeros(nu,mpol+1);
sinv=zeros(nv,2*ntor+1);
alu=2*pi/nu;
for i=1:nu
    for j=1:mpol+1
        m=j-1;
        cosu(i,j)=cos(m*(i-1)*alu);
        sinu(i,j)=sin(m*(i-1)*alu);
    end
end
alv=2*pi/nv;
for i=1:nv
    for j=1:2*ntor+1
        n=j-ntor-1;
        cosv(i,j)=cos(n*(i-1)*alv);
        sinv(i,j)=sin(n*(i-1)*alv);
    end
end
disp('      Calculating Surface Fourier Transform');
% Our Amplitude is 1 everywhere since we do a full period
% Our Kernels is
%    read_VMEC arrays
%    cos(mu+nv) = cos(mu)cos(nv)-sin(mu)sin(nv)
%    sin(mu+nv) = sin(mu)cos(nv)+cos(mu)sin(nv)
fnuv=zeros(1,mpol+1);
fnuv(1)=1./(nu*nv);
for i=2:mpol+1;
    fnuv(i)=2*fnuv(1);
end
for m1=1:mpol+1
    for n1=1:2*ntor+1
        for i=1:nu
            for j=1:nv
                refou(:,m1,n1)=refou(:,m1,n1)+r_pies(:,i,j).*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
                zefou(:,m1,n1)=zefou(:,m1,n1)+z_pies(:,i,j).*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
            end
        end
    end
end

% Now add the VMEC nested field to the outter surface
disp('  - Calculating Field');
bu=cfunct(theta,zeta,data.buc,data.nfp);
bv=cfunct(theta,zeta,data.bvc,data.nfp);
bv_vmec=bv;
bu_vmec=bu;
bs_vmec=bu.*0;
bxe=zeros(nu,nv);
bye=zeros(nu,nv);
bze=zeros(nu,nv);
xe=zeros(nu,nv);
ye=zeros(nu,nv);
for i=1:nv
    xe(:,i)=re(:,i).*cos(zeta(i));
    ye(:,i)=re(:,i).*sin(zeta(i));
end
disp('      Calculating Coil Field');
coildata=coil_biot_prep(coildata);
for i=1:nu
    if i==1, tic; end
    for j=1:nv
        [bxe(i,j) bye(i,j) bze(i,j)]=coil_biot(coildata,xe(i,j),ye(i,j),ze(i,j),data.extcur);
    end
    if i==1
        time=toc;
        time=time*(nu-1);
        if time < 60
            disp(['      Est time to Calculate: ' num2str(time) ' [s]']);
        else
            disp(['      Est time to Calculate: ' num2str(time/60) ' [min]']);
        end
    end
end
bre=zeros(nu,nv);
bphie=zeros(nu,nv);
for i=1:nv
    bre(:,i)=bxe(:,i).*cos(zeta(i))+bye(:,i).*sin(zeta(i));
    bphie(:,i)=-bxe(:,i).*sin(zeta(i))+bye(:,i).*cos(zeta(i));
end
bve=bphie./re;
bv_vmec(data.ns+1,:,:)=bve;
% Now get derivatives on the new surface
disp('      Calculating Field Fourier Transform');
ruse=zeros(mpol+1,2*ntor+1);
rvse=zeros(mpol+1,2*ntor+1);
zuce=zeros(mpol+1,2*ntor+1);
zvce=zeros(mpol+1,2*ntor+1);
for i=1:mpol+1
    for j=1:2*ntor+1
        m=i-1;
        n=j-ntor-1;
        ruse(i,j)=-m.*refou(ns,i,j);
        rvse(i,j)=-n.*refou(ns,i,j);
        zuce(i,j)=m.*zefou(ns,i,j);
        zvce(i,j)=n.*zefou(ns,i,j);
    end
end
drdue=squeeze(sfunct(theta,zeta,ruse,data.nfp));
dzdue=squeeze(cfunct(theta,zeta,zuce,data.nfp));
% Now calculate bs and bu
bse=zeros(nu,nv);
bue=zeros(nu,nv);
for i=1:nv
    bse(:,i)=bxe(:,i).*drds(:,i).*cos(zeta(i))+...
        bye(:,i).*drds(:,i).*sin(zeta(i))+...
        bze(:,i).*dzds(:,i);
    bue(:,i)=bxe(:,i).*drdue(:,i).*cos(zeta(i))+...
        bye(:,i).*drdue(:,i).*sin(zeta(i))+...
        bze(:,i).*dzdue(:,i);
end
bu_vmec(data.ns+1,:,:)=bue;
bs_vmec(data.ns+1,:,:)=bse;
disp('      Interpolating in Real Space');
% For each pies surface we interpolate over real space in rho
bs_pies=zeros(ns,nu,nv);
bu_pies=zeros(ns,nu,nv);
bv_pies=zeros(ns,nu,nv);
for i=1:nu
    for j=1:nv
        bs_pies(:,i,j)=interp1(rho(surfmask,i,j),bs_vmec(surfmask,i,j),rho_pies(:,i,j),'pchip');
        bu_pies(:,i,j)=interp1(rho(surfmask,i,j),bu_vmec(surfmask,i,j),rho_pies(:,i,j),'pchip');
        bv_pies(:,i,j)=interp1(rho(surfmask,i,j),bv_vmec(surfmask,i,j),rho_pies(:,i,j),'pchip');
    end
end
% Now we Fourier Transform the Fields
bsfou=zeros(ns,mpol+1,2*ntor+1);
bufou=zeros(ns,mpol+1,2*ntor+1);
bvfou=zeros(ns,mpol+1,2*ntor+1);
for m1=1:mpol+1
    for n1=1:2*ntor+1
        for i=1:nu
            for j=1:nv
                bsfou(:,m1,n1)=bsfou(m1,n1)+bs_pies(:,i,j)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
                bufou(:,m1,n1)=bufou(:,m1,n1)+bu_pies(:,i,j).*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
                bvfou(:,m1,n1)=bvfou(:,m1,n1)+bv_pies(:,i,j).*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
            end
        end
    end
end


disp('  - Transforming to PIES Coordinates');
% Transform to PIES fourier space
rtemp=zeros(ns,mval+1,2*nval+1);
ztemp=zeros(ns,mval+1,2*nval+1);
bstemp=zeros(ns,mval+1,2*nval+1);
butemp=zeros(ns,mval+1,2*nval+1);
bvtemp=zeros(ns,mval+1,2*nval+1);
for i=1:ns
    for j=1:2*ntor+1
        ntemp=j+(nval-ntor);
        rtemp(i,1:mpol+1,ntemp)=refou(i,:,j);
        ztemp(i,1:mpol+1,ntemp)=zefou(i,:,j);
        bstemp(i,1:mpol+1,ntemp)=bsfou(i,:,j);
        butemp(i,1:mpol+1,ntemp)=bufou(i,:,j);
        bvtemp(i,1:mpol+1,ntemp)=bvfou(i,:,j);
    end
end
refou=rtemp;
zefou=ztemp;
bsfou=bstemp;
bufou=butemp;
bvfou=bvtemp;
n0=1+round(nval/2);

% Clean up R
refou= flipdim(refou(:,:,:),3);
zefou=-flipdim(zefou(:,:,:),3);
bsfou=-flipdim(bsfou(:,:,:),3);
bufou= flipdim(bufou(:,:,:),3);
bvfou= flipdim(bvfou(:,:,:),3);

refou(1,2:mval+1,:)=0.0;
zefou(1,2:mval+1,:)=0.0;
zefou(1,1,n0)=0.0;
piesinput.rmaj=refou(1,1,n0);
refou(:,1,n0)=refou(:,1,n0)-piesinput.rmaj; % Center about magnetic axis
xadj=1-sum(sum(refou(ns,:,:)));
refou(:,1,n0)=refou(:,1,n0)+xadj;
piesinput.rmaj=piesinput.rmaj-xadj;

% Create coils file
if makecoil
    disp('  - Creating Coils file');
    coil_pies(coildata,data.extcur);
end

% Output PIES Namelists
disp('  - Creating PIES Input File');
piesfile=strtrim(['pies.' strtrim(filename) '.in']);
fid=fopen(piesfile,'w+');
% Rerun line
fprintf(fid,'''%s''\n',strtrim(run_type));
% Input Namelist
write_fortran_namelist(fid,piesinput,'INPUT');
% PLTFLG
write_fortran_namelist(fid,pltflg,'PLTFLG');
% EXLSTA
write_fortran_namelist(fid,exlsta,'EXLSTA');
% Output VMEC Data
for i=1:ns
    for j=1:2*nval+1
        fprintf(fid,' %20.10e ',refou(i,:,j));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    for j=1:2*nval+1
        fprintf(fid,' %20.10e ',zefou(i,:,j));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');
end
for i=1:ns
    for j=1:2*nval+1
        fprintf(fid,' %20.10e ',bsfou(i,:,j));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
for i=1:ns
    for j=1:2*nval+1
        fprintf(fid,' %20.10e ',bufou(i,:,j));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
for i=1:ns
    for j=1:2*nval+1
        fprintf(fid,' %20.10e ',bvfou(i,:,j));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
% Output Mode selection matrix
for i=1:2*nval+1
    fprintf(fid,' %d',modes(:,i));
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
% Output pressure Profile
if usespline == 1   % Output splined values
    pres_pies=zeros(1,ns);
    jcur_pies=zeros(1,ns);
    % First interpolate p and j to pies grid
    for i=1:ns
        if rpies(i) < rv
            pres_pies(i)=pchip(rvmec,data.presf,rpies(i));
            jcur_pies(i)=pchip(rvmec,data.jcurv,rpies(i));
        else
            pres_pies(i)=0.0;
            jcur_pies(i)=0.0;
        end            
    end
    % Create Splines
    p_spline=spline(0:1/(ns-1):1,pres_pies);
    j_spline=spline(0:1/(ns-1):1,jcur_pies);
    % Output to file
    fprintf(fid,'\n');
    for i=1:p_spline.pieces
        fprintf(fid,'%20.10e\n',p_spline.breaks(i));
    end
    fprintf(fid,'\n');
    for i=1:j_spline.pieces
        fprintf(fid,'%20.10e\n',j_spline.breaks(i));
    end
    fprintf(fid,'\n');
    for i=1:p_spline.pieces
        for j=p_spline.order:-1:1  % Reverse order because matlab does things backwards
            fprintf(fid,' %20.10e ',p_spline.coefs(i,j));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    for i=1:j_spline.pieces
        %for j=1:j_spline.order
        for j=j_spline.order:-1:1  % Reverse order because matlab does things backwards
            fprintf(fid,' %20.10e ',j_spline.coefs(i,j));
        end
        fprintf(fid,'\n');
    end
    % Plot Spline
    fig=figure;
    subplot(1,2,1);
    plot(rvmec,data.presf,'ok');
    hold on
    plot(rpies,ppval(p_spline,0:1/(ns-1):1),'r');
    hold off
    axis tight
    title('Pressure');
    subplot(1,2,2);
    plot(rvmec,data.jcurv,'ok');
    hold on
    plot(rpies,ppval(j_spline,0:1/(ns-1):1),'r');
    hold off
    axis tight
    title('Toroidal Current');
    legend('VMEC','PIES SPLINE');
    pause(1.0);
    saveas(gcf,['splines_' strtrim(filename) '.fig']);
    close(fig);
else  % Output AC and AM Arrays
    if isfield(data,'ac')
        fprintf(fid,'\n');
        fprintf(fid,' %20.10e',data.ac(1:11));
        fprintf(fid,'\n');
        fprintf(fid,' %20.10e',data.am(1:11));
        fprintf(fid,'\n');
    else
        disp('ACAM_MORPHED FILE EMPTY!');
    end
end
fclose(fid);

% Make a total diagnostic plot
zeta=[0 .5*pi/data.nfp pi/data.nfp];
theta=0:2*pi/359:2*pi;
rbc=permute(refou(:,:,:),[2 3 1]);
zbs=permute(zefou(:,:,:),[2 3 1]);
r=cfunct(-theta,zeta,rbc,data.nfp);
z=sfunct(-theta,zeta,zbs,data.nfp);
figure('Position',[1 1 1024 768],'Color','white');
% Plot Flux Surfaces
for i=1:3
    subplot(2,3,i);
    plot(r(1,1,i),z(1,1,i),'k+')
    hold on
    for j=1:piesinput.lpinch
        plot(r(j,:,i),z(j,:,i),'k')
    end
    for j=piesinput.lpinch:piesinput.k
        plot(r(j,:,i),z(j,:,i),'b')
    end
    j=ns;
    plot(r(j,:,i),z(j,:,i),'r')
    %for j=vmecdata.ns+1:piesinput.k
    %    plot(r(j,:,i),z(j,:,i),'r')
    %end
    hold off
    axis tight; axis square; axis equal;
end
% Plot Mode Selection field
subplot(2,3,6);
if nval > 0
    pixplot(double(0:mval),double(-nval:nval),double(modes));
else
    bar(double(0:mval),modes(:,1));
end
xlabel('m');
ylabel('n');
title('Mode Selection Matrix');
% Plot current and pressure profile

return

end

