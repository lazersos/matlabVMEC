function pies_input(varargin)
%PIES_INPUT Creates a pies input file.
%   This function creates a PIES input file.  It can accept a VMEC data
%   structure, VMEC input data structure, coil data structure, pies_input
%   data strucure, pltflg data structure, and exlsta data structure.  The
%   user can also specifiy the filename, Fourier resolution (mpol,ntor),
%   number of surfaces, boundary type, last iteration, and scale factor for
%   the new boundary.  It is suggested that the user play around with a
%   scale factor between 0.5 and 2 for the boundary.
%
%   Example:
%       vmec_data=read_vmec('wout.test');
%       vmec_input=read_vmec_input('input.test');
%       coil_data=read_coil('coils.test_machine');
%       plot_data=read_namelist('default.exa','PLTFLG');
%       plot_data.datatype='pltflg';  % Do this for now
%       pies_input(vmec_data,vmec_input,coil_data,plot_data,...
%           'filename','test.in','mpol',24,'ntor',12,'ns',99,'scale',0.5);
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        0.1
%   Date:           3/22/11

% Defaults
i=1;
iter2=5000;
kval=50;
mval=6;
nval=6;
mu0=4*pi*1.e-7;
vmecdata=[];
coildata=[];
piesinput=[];
pltflg=[];
exlsta=[];
vmecinput=[];
bound_type='expand';
run_type='begin';
filename='pies.in';
scale=1;
novac=0;

% Handle varargin
while i <= nargin
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch varargin{i}.datatype
                case 'VMEC_input'
                    vmecinput=varargin{i};
                case 'wout'
                    vmecdata=varargin{i};
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
            case 'ns'
                i=i+1;
                kval=varargin{i};
            case 'bound_type'
                i=i+1;
                bound_type=varargin{i};
            case 'scale'
                i=i+1;
                scale=varargin{i};
            case 'iter2'
                i=i+1;
                iter2=varargin{i};
            case 'novac'
                novac=1;
        end
    end
    i=i+1;
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
    modes=zeros(piesinput.m+1,2*piesinput.n+1);
    for i=1:piesinput.m+1
        m=i-1;
        if (m <= piesinput.m/2)
            for j=1:2*piesinput.n+1
                n=j-piesinput.n-1;
                if ((n >= - piesinput.n/2) && (n <= piesinput.n/2))
                    modes(i,j)=1;
                end
            end
        end
    end
    piesinput.bblnd=0;
    piesinput.vmecf=int32(0);
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
    piesinput.ispln=int32(1);
    piesinput.lp=int32(6);
    piesinput.lj=int32(6);
    piesinput.rtokfk=1.0;
    piesinput.adjst=int32(0);
    piesinput.rmaj=3.0;
    piesinput.nper=3.0;
    piesinput.cyl=int32(0);
    piesinput.island=int32(1);
    piesinput.tokchk=int32(0);
    piesinput.rot1=0.14;
    piesinput.rot2=0.00;
    piesinput.setbc=int32(0);
    piesinput.freeb=int32(0);
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
    exlsta.blend_B=0.98;
    exlsta.check_point_flag=int32(0);
    exlsta.out_of_phase=int32(0);
    exlsta.boot_model_f=int32(0);
    exlsta.isledgf=int32(1);
    exlsta.ibisec=int32(7);
    exlsta.iftmth=int32(3);
    exlsta.storvac=int32(1);
    exlsta.nw=int32(0);
    exlsta.virtual_casing=int32(1);
    exlsta.fac_nescoil=-1.000000011686097E-007;
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
end

% Check for VMEC data
if ~isempty(vmecdata)
    % Check for assymetry
    if vmecdata.iasym
        disp('ERROR: Not an up-down symmetric VMEC run');
        return
    end
    if isempty(coildata) && vmecdata.lfreeb
        disp('ERROR: Coil datafile required!');
        disp(['       Try file associated with: ' vmecdata.mgrid_file]);
        return
    end
    % Set some defaults
    piesinput.iter2=int32(iter2);
    piesinput.vmecf=int32(1);
    piesinput.k=kval;
    if mval <= 2 * (vmecdata.mpol-1)
        disp('Warning:  mpies < 2 * mvmec');
    end
    if nval <= 2 * (vmecdata.ntor-1)
        disp('Warning:  npies < 2 * nvmec');
    end
    piesinput.m=mval;
    piesinput.n=nval;
    piesinput.nper=vmecdata.nfp;
    piesinput.mda=int32(mval);  %Modeselect = 1/2 mval
    piesinput.nda=int32(nval);  %Modeselect = 1/2 nval
    piesinput.iotamx=1.2*max(vmecdata.iotaf);
    piesinput.iotamn=0.2*min(vmecdata.iotaf);
    piesinput.rot1=vmecdata.iotaf(1);
    piesinput.rot2=vmecdata.iotaf(vmecdata.ns)-vmecdata.iotaf(1);
    % Calculate beta (Pressure scaling factor)
    piesinput.beta=mu0;
    piesinput.lp=int32(kval+1);
    piesinput.lj=int32(kval+1);
    piesinput.ispln=int32(2);
    % bsubpi (toroidal field on axis)
    piesinput.bsubpi=vmecdata.RBtor;
    % Calculate betai (scaling factor for current)
    phiedge=vmecdata.phi(vmecdata.ns);
    iint=0.0;
    if isfield(vmecdata,'ac')
        for i=1:11
            iint=iint+vmecdata.ac(i)/double(i);
        end
    else
        for i=1:11
            iint=iint+vmecinput.ac(i)/double(i);
        end
    end
    if (iint*phiedge == 0.0)
        piesinput.betai=0.0;
    else
        piesinput.betai=-(vmecdata.Itor*mu0)/(iint*phiedge);
    end
    % Set some More Default values
    exlsta.vmecbf=int32(1);
    exlsta.use_poly_for_currp_press=int32(1);
    exlsta.vmec_ignore_sym=int32(0);
    if vmecdata.lfreeb
        piesinput.freeb=int32(1);
        piesinput.lpinch=int32(piesinput.k-1);  %Make the last VMEC flux surface the limiter
        %piesinput.rtokfk=double(piesinput.lpinch)/double(piesinput.k);  %Plasma vacuum interface is the limiter surface
        piesinput.setbc=int32(1);
        exlsta.storvac=int32(1);
        exlsta.virtual_casing=int32(1);
        exlsta.read_brho=int32(1);
    else
        piesinput.freeb=int32(0);
        piesinput.lpinch=int32(0);
        piesinput.setbc=int32(0);
        exlsta.storvac=int32(0);
        exlsta.virtual_casing=int32(0);
        exlsta.setbc_override=int32(1);
    end
    % Now we convert from VMEC to PIES ordering
    rbc=permute(vmecdata.rbc,[3 1 2]);
    zbs=permute(vmecdata.zbs,[3 1 2]);
    bsupuc=permute(vmecdata.buc,[3 1 2]);
    bsupvc=permute(vmecdata.bvc,[3 1 2]);
    % Note that in VMEC m=0:mpol-1
    ntor=vmecdata.ntor;
    mpol=vmecdata.mpol-1;
    % Flush out (0,n)
    for i=1:ntor
        n2=2*ntor+2-i;
        rbc(:,1,i)=rbc(:,1,i)/2.;
        zbs(:,1,i)=zbs(:,1,i)/2.;
        bsupuc(:,1,i)=bsupuc(:,1,i)/2.;
        bsupvc(:,1,i)=bsupvc(:,1,i)/2.;
        rbc(:,1,n2)=rbc(:,1,i);
        zbs(:,1,n2)=zbs(:,1,i);
        bsupuc(:,1,n2)=bsupuc(:,1,i);
        bsupvc(:,1,n2)=bsupvc(:,1,i);
    end
    if vmecdata.lfreeb
        % Convert to centered coordinates
        piesinput.rmaj=rbc(1,1,ntor+1);
        rbc(:,1,ntor+1)=rbc(:,1,ntor+1)-rbc(1,1,ntor+1);
        % Now we extend the VMEC domain
        % Initially assume we have 1 more surface than VMEC
        x=zeros(vmecdata.ns+1,piesinput.m+1,2*piesinput.n+1,2);
        b=zeros(vmecdata.ns+1,piesinput.m+1,2*piesinput.n+1,3);
        for i=1:vmecdata.ns
            for j=1:mpol+1
                for k=1:2*ntor+1
                    i2=i;
                    j2=j;
                    k2=k+(piesinput.n-ntor);
                    x(i2,j2,k2,1)=rbc(i,j,k);
                    x(i2,j2,k2,2)=zbs(i,j,k);
                    b(i2,j2,k2,1)=bsupuc(i,j,k);
                    b(i2,j2,k2,2)=bsupvc(i,j,k);
                end
            end
        end
        % Now we add 1 surface (at s~1+scale*ds)
        [x b]=cartextend(vmecdata,coildata,x,b,vmecdata.ns,scale,bound_type,novac);
        
        % Now we convert back to non-centered coordinates
        x(:,1,piesinput.n+1,1)=x(:,1,piesinput.n+1,1)+piesinput.rmaj;
        
        % Mode selection matrix
        modes=zeros(piesinput.m+1,2*piesinput.n+1);
        for i=1:piesinput.m+1
            m=i-1;
            if (m <= mval/2)
                for j=1:2*piesinput.n+1
                    n=j-piesinput.n-1;
                    if ((n >= - nval/2) && (n <= nval/2))
                        modes(i,j)=1;
                    end
                end
            end
        end
        disp('  - Interpolating to new grid');
        % Now we need to interpolate from s to k space.
        svec=zeros(1,vmecdata.ns+1);
        for i=1:vmecdata.ns+1
            svec(i)=sum(x(i,:,piesinput.n+1,1))-piesinput.rmaj;
        end
        kvec=0:max(svec)/piesinput.k:max(svec);
        %svec=0:1/(vmecdata.ns):1;
        %svec=sqrt(svec);
        %kvec=0:1/(piesinput.k):1;
        %svec_last=1+scale/100;
        %svec=[svec svec_last];
        %svec_max=max(svec);
        %svec=svec./svec_max;
        %kvec=0:svec_last/(piesinput.k):svec_last;
        x2=zeros(piesinput.k+1,piesinput.m+1,2*piesinput.n+1,2);
        b2=zeros(piesinput.k+1,piesinput.m+1,2*piesinput.n+1,3);
        for i=1:piesinput.m+1
            for j=1:2*piesinput.n+1
                if mode(i,j)
                    x2(:,i,j,1)=interp1(svec,x(:,i,j,1),kvec,'pchip');
                    x2(:,i,j,2)=interp1(svec,x(:,i,j,2),kvec,'pchip');
                    b2(:,i,j,1)=interp1(svec,b(:,i,j,1),kvec,'pchip');
                    b2(:,i,j,2)=interp1(svec,b(:,i,j,2),kvec,'pchip');
                    b2(:,i,j,3)=interp1(svec,b(:,i,j,3),kvec,'pchip');
                end
            end
        end
        x=x2;
        b=b2;
        % Fix up rmajor
        piesinput.rmaj=sum(sum(x(piesinput.k+1,:,:,1)))-1.0;
        x(:,1,piesinput.n+1,1)=x(:,1,piesinput.n+1,1)-piesinput.rmaj;
        % Pick lpinch correctly
        x0=zeros(0,piesinput.k+1);
        for i=1:piesinput.k+1
            x0(i)=sum(x(i,:,piesinput.n+1,1))+piesinput.rmaj;
        end
        xvmec=sum(vmecdata.rbc(:,vmecdata.ntor+1,vmecdata.ns));
        for i=1:piesinput.k+1
            if (x0(i) <= xvmec)
                piesinput.lpinch=int32(i-1);
                %piesinput.rtokfk=double(i/piesinput.k+1);
            end
        end
    else
        % Convert to centered coordinates
        rbc(:,1,ntor+1)=rbc(:,1,ntor+1)-vmecdata.rmajor;
        piesinput.rmaj=vmecdata.rmajor;
        % First we put the VMEC data on a pies fourier grid
        x=zeros(vmecdata.ns,piesinput.m+1,2*piesinput.n+1,2);
        b=zeros(vmecdata.ns,piesinput.m+1,2*piesinput.n+1,2);
        for i=1:vmecdata.ns
            for j=1:mpol+1
                for k=1:2*ntor+1
                    i2=i;
                    j2=j;
                    k2=k+(piesinput.n-ntor);
                    x(i2,j2,k2,1)=rbc(i,j,k);
                    x(i2,j2,k2,2)=zbs(i,j,k);
                    b(i2,j2,k2,1)=bsupuc(i,j,k);
                    b(i2,j2,k2,2)=bsupvc(i,j,k);
                end
            end
        end
        % Mode selection matrix
        modes=zeros(piesinput.m+1,2*piesinput.n+1);
        for i=1:piesinput.m+1
            m=i-1;
            if (m <= mval/2)
                for j=1:2*piesinput.n+1
                    n=j-piesinput.n-1;
                    if ((n >= - nval/2) && (n <= nval/2))
                        modes(i,j)=1;
                    end
                end
            end
        end
        disp('  - Interpolating to new grid');
        % Now we need to interpolate from s to k space.
        svec=0:1/(vmecdata.ns-1):1;
        svec=sqrt(svec);
        kvec=0:1/(piesinput.k):1;
        %kvec=0:svec_last*svec_last/(piesinput.k):svec_last*svec_last;
        x2=zeros(piesinput.k+1,piesinput.m+1,2*piesinput.n+1,2);
        b2=zeros(piesinput.k+1,piesinput.m+1,2*piesinput.n+1,2);
        for i=1:piesinput.m+1
            for j=1:2*piesinput.n+1
                if mode(i,j)
                    x2(:,i,j,1)=interp1(svec,x(:,i,j,1),kvec,'spline');
                    x2(:,i,j,2)=interp1(svec,x(:,i,j,2),kvec,'spline');
                    b2(:,i,j,1)=interp1(svec,b(:,i,j,1),kvec,'spline');
                    b2(:,i,j,2)=interp1(svec,b(:,i,j,2),kvec,'spline');
                end
            end
        end
        x=x2;
        b=b2;
    end
    disp('  - Fixing indexing and cleaning up');
    % Transform to VMEC coordinates then invert sin quantities to end up in
    % pies (mu+nv) -> (mu-nv) -> (nv-mu)
    x(:,:,:,1)=flipdim(x(:,:,:,1),3);   % cos(-x)=cos(x)
    x(:,:,:,2)=-flipdim(x(:,:,:,2),3);  % sin(-x)=-sin(x)
    b(:,:,:,1)=flipdim(b(:,:,:,1),3);   % cos(-x)=cos(x)
    b(:,:,:,2)=flipdim(b(:,:,:,2),3);   % cos(-x)=cos(x)
    b(:,:,:,3)=-flipdim(b(:,:,:,3),3);   % cos(-x)=cos(x)
    
    % Now do some cleaning up
    %x(1,2:2*piesinput.m+1,:,:)=0.0;
    %if abs(x(piesinput.k,2,piesinput.n+1,1)) < 1e-10
    %    x(piesinput.k,2,piesinput.n+1,1)=1e-9;
    %end
    %if abs(x(piesinput.k,2,piesinput.n+1,2)) < 1e-10
    %    x(piesinput.k,2,piesinput.n+1,2)=1e-9;
    %end
end


disp('  - Input parameter check');
% This section is to do a check of parameters before writing
piesinput.k=int32(piesinput.k);
piesinput.m=int32(piesinput.m);
piesinput.n=int32(piesinput.n);
piesinput.nper=double(piesinput.nper);
piesinput.cyl=boolean(piesinput.cyl);
piesinput.poinc=boolean(piesinput.poinc);

% Get the pressure profile vector
svec=-1:1/(vmecdata.ns-1):1;
jvec=[fliplr(vmecdata.jcurv) vmecdata.jcurv(2:vmecdata.ns)];
pvec=[fliplr(vmecdata.presf) vmecdata.presf(2:vmecdata.ns)];
if ~novac
    % Open and output the file
    fid=fopen(filename,'w');
    % Rerun line
    fprintf(fid,'''%s''\n',strtrim(run_type));
    % Input Namelist
    write_fortran_namelist(fid,piesinput,'INPUT');
    % PLTFLG
    write_fortran_namelist(fid,pltflg,'PLTFLG');
    % EXLSTA
    write_fortran_namelist(fid,exlsta,'EXLSTA');
    % Output VMEC Data
    if ~isempty(vmecdata)
        for i=1:piesinput.k+1
            for j=1:2*piesinput.n+1
                fprintf(fid,' %20.10e ',x(i,:,j,1));
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
            for j=1:2*piesinput.n+1
                fprintf(fid,' %20.10e ',x(i,:,j,2));
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
            fprintf(fid,'\n');
        end
        if (vmecdata.lfreeb)
            for i=1:piesinput.k+1
                for j=1:2*piesinput.n+1
                    fprintf(fid,' %20.10e ',b(i,:,j,3));
                    fprintf(fid,'\n');
                end
                fprintf(fid,'\n');
            end
        end
        for i=1:piesinput.k+1
            for j=1:2*piesinput.n+1
                fprintf(fid,' %20.10e ',b(i,:,j,1));
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
        end
        for i=1:piesinput.k+1
            for j=1:2*piesinput.n+1
                fprintf(fid,' %20.10e ',b(i,:,j,2));
                fprintf(fid,'\n');
            end
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    % Output Mode selection matrix
    if (piesinput.mdslct >0)
        for i=1:2*piesinput.n+1
            fprintf(fid,' %d',modes(:,i));
            fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    % Output VMEC pressure and current
    if isempty(vmecinput) && ~isfield(vmecdata,'ac') && ~isfield(vmecdata,'am')
        jcoef=fliplr(polyfit(svec,jvec,10));
        j = polyval(fliplr(jcoef),svec);
        fprintf(fid,' %20.10e ',jcoef(:));
        fprintf(fid,'\n');
        pcoef=fliplr(polyfit(svec,pvec,10));
        p = polyval(fliplr(pcoef),svec);
        fprintf(fid,' %20.10e ',pcoef(:));
        fprintf(fid,'\n');
    elseif isempty(vmecinput) && isfield(vmecdata,'ac') && isfield(vmecdata,'am')
        fprintf(fid,' %20.10e ',vmecdata.ac(1:11));
        fprintf(fid,'\n');
        fprintf(fid,' %20.10e ',vmecdata.am(1:11));
        fprintf(fid,'\n');
    elseif ~isempty(vmecinput)
        fprintf(fid,' %20.10e ',vmecinput.ac(1:11).*vmecinput.curtor);
        fprintf(fid,'\n');
        fprintf(fid,' %20.10e ',vmecinput.am(1:11));
        fprintf(fid,'\n');
    else
        disp('ERROR: No Current or Pressure profile specification');
    end
    fclose(fid);
end
%Make some diagnostic plots
if ~isempty(vmecdata)
    zeta=[0 .5*pi/vmecdata.nfp pi/vmecdata.nfp];
    theta=0:2*pi/359:2*pi;
    rbc=permute(x(:,:,:,1),[2 3 1]);
    zbs=permute(x(:,:,:,2),[2 3 1]);
    r=cfunct(-theta,zeta,rbc,vmecdata.nfp);
    z=sfunct(-theta,zeta,zbs,vmecdata.nfp);
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
        j=piesinput.k+1;
        plot(r(j,:,i),z(j,:,i),'r')
        %for j=vmecdata.ns+1:piesinput.k
        %    plot(r(j,:,i),z(j,:,i),'r')
        %end
        hold off
        axis tight; axis square; axis equal;
    end
    % Plot Mode Selection field
    subplot(2,3,6);
    pixplot(double(0:piesinput.m),double(-piesinput.n:piesinput.n),double(modes));
    xlabel('m');
    ylabel('n');
    title('Mode Selection Matrix');
    % Plot current and pressure profile
    if isempty(vmecinput) && ~isfield(vmecdata,'ac') && ~isfield(vmecdata,'am')
        subplot(2,3,4)
        plot(svec,jvec,'k');
        hold on
        plot(svec,j,'ro');
        title('VMEC Current Profile');
        hold off
        subplot(2,3,5)
        plot(svec,pvec,'k');
        hold on
        plot(svec,p,'ro');
        title('VMEC Pressure Profile');
    elseif isempty(vmecinput) && isfield(vmecdata,'ac') && isfield(vmecdata,'am')
        svec=0:1/(vmecdata.ns-1):1;
        j=polyval(fliplr(vmecdata.ac(1:11)),svec)./(2*pi);
        p=polyval(fliplr(vmecdata.am(1:11)),svec);
        svec=-1:1/(vmecdata.ns-1):1;
        j=[fliplr(j) j(2:vmecdata.ns)];
        p=[fliplr(p) p(2:vmecdata.ns)];
        subplot(2,3,4)
        plot(svec,jvec,'k');
        hold on
        plot(svec,j,'ro');
        title('VMEC Current Profile');
        hold off
        subplot(2,3,5)
        plot(svec,pvec,'k');
        hold on
        plot(svec,p,'ro');
        title('VMEC Pressure Profile');
    elseif ~isempty(vmecinput)
        svec=0:1/(vmecdata.ns-1):1;
        j=polyval(vmecinput.curtor*fliplr(vmecinput.ac(1:11)),svec);
        p=polyval(fliplr(vmecinput.am(1:11)),svec);
        svec=-1:1/(vmecdata.ns-1):1;
        j=[fliplr(j) j(2:vmecdata.ns)];
        p=[fliplr(p) p(2:vmecdata.ns)];
        subplot(2,3,4)
        plot(svec,jvec,'k');
        hold on
        plot(svec,j,'ro');
        title('VMEC Current Profile');
        hold off
        subplot(2,3,5)
        plot(svec,pvec,'k');
        hold on
        plot(svec,p,'ro');
        title('VMEC Pressure Profile');
    end    
    subplot(2,3,6)
end
end

function [x b]=cartextend(data,coildata,x,b,ns,scale,type,novac)
% Note that the rbc and zbs arrays are in k,m,n format.
mpol=size(data.rbc,1)-1;
ntor=(size(data.rbc,2)-1)/2;
piesn=(size(b,3)-1)/2;
nu=round(5*mpol+2);
nv=round(5*ntor)+2;
theta=0:2*pi/(nu-1):2*pi;
if (ntor == 0 && data.nfp == 1)  % Axisymmetry
    zeta=0.0;
    nv=1;
    % Get the outter surface (note we're in centered coordinates)
    r0=squeeze(cfunct(theta,zeta,permute(x(ns,:,:,1),[2 3 1]),data.nfp))';
    z0=squeeze(sfunct(theta,zeta,permute(x(ns,:,:,2),[2 3 1]),data.nfp))';
    % Now just get the whole surface
    r=squeeze(cfunct(theta,zeta,data.rbc(:,:,data.ns),data.nfp))';
    z=squeeze(sfunct(theta,zeta,data.zbs(:,:,data.ns),data.nfp))';
    %bu=squeeze(cfunct(theta,zeta,data.buc(:,:,data.ns),data.nfp))';
    %bv=squeeze(cfunct(theta,zeta,data.bvc(:,:,data.ns),data.nfp))';
    %drdu=squeeze(sfunct(theta,zeta,data.rus(:,:,data.ns),data.nfp))';
    %drdv=squeeze(sfunct(theta,zeta,data.rvs(:,:,data.ns),data.nfp))';
    %dzdu=squeeze(cfunct(theta,zeta,data.zuc(:,:,data.ns),data.nfp))';
    %dzdv=squeeze(cfunct(theta,zeta,data.zvc(:,:,data.ns),data.nfp))';
    rax=squeeze(cfunct(theta,zeta,data.rbc(:,:,1),data.nfp))';
    zax=squeeze(sfunct(theta,zeta,data.zbs(:,:,1),data.nfp))';
else
    zeta=0:2*pi/(data.nfp*(nv-1)):2*pi/data.nfp;
    % Get the outter surface (note we're in centered coordinates)
    r0=squeeze(cfunct(theta,zeta,permute(x(ns,:,:,1),[2 3 1]),data.nfp));
    z0=squeeze(sfunct(theta,zeta,permute(x(ns,:,:,2),[2 3 1]),data.nfp));
    % Now just get the whole surface
    r=squeeze(cfunct(theta,zeta,data.rbc(:,:,data.ns),data.nfp));
    z=squeeze(sfunct(theta,zeta,data.zbs(:,:,data.ns),data.nfp));
    %bu=squeeze(cfunct(theta,zeta,data.buc(:,:,data.ns),data.nfp));
    %bv=squeeze(cfunct(theta,zeta,data.bvc(:,:,data.ns),data.nfp));
    %drdu=squeeze(sfunct(theta,zeta,data.rus(:,:,data.ns),data.nfp));
    %drdv=squeeze(sfunct(theta,zeta,data.rvs(:,:,data.ns),data.nfp));
    %dzdu=squeeze(cfunct(theta,zeta,data.zuc(:,:,data.ns),data.nfp));
    %dzdv=squeeze(cfunct(theta,zeta,data.zvc(:,:,data.ns),data.nfp));
    rax=squeeze(cfunct(theta,zeta,data.rbc(:,:,1),data.nfp));
    zax=squeeze(sfunct(theta,zeta,data.zbs(:,:,1),data.nfp));
end
%Now calculate the surface curent (note sign flip on currents)
%br=bu.*drdu+bv.*drdv;
%bphi=bv.*r;
%bz=bu.*dzdu+bv.*dzdv;
%nr=r.*dzdu;
%nphi=drdu.*dzdv-drdv.*dzdu;
%nz=r.*drdu;
%jr=-(bphi.*nz-bz.*nphi);
%jphi=-(bz.*nr-br.*nz);
%jz=-(br.*nphi-bphi.*nr);
% Now allocate the expanded arrays
re=zeros(nu,nv);
ze=zeros(nu,nv);
re0=zeros(nu,nv);
ze0=zeros(nu,nv);
% Now calculate re and ze from r's and z's
if strcmp(type,'D')
    %dr=(r-rax)/(data.ns-1);
    %dz=(z-zax)/(data.ns-1);
    %drds=dr;
    %dzds=dz;
    %rmax=max(max(r));
    %rmin=min(min(r));
    %zmax=max(max(z));
    %zmin=min(min(z));
elseif strcmp(type,'torus')
    dr=(r-rax)/(data.ns-1);
    dz=(z-zax)/(data.ns-1);
    drds=dr;
    dzds=dz;
    rmax=max(max(r));
    rmin=min(min(r));
    zmax=max(max(z));
    zmin=min(min(z));
    r_center=(rmax+rmin)/2;
    radius=max([abs(rmax-r_center) abs(zmax) abs(r_center-rmin) abs(zmin)])*(1+scale/100);
    for i=1:nu
        for j=1:nv
            re0(i,j)=radius*cos(theta(i));
            ze0(i,j)=radius*sin(theta(i));
            re(i,j)=re0(i,j)+r_center;
            ze(i,j)=ze0(i,j);
        end
    end
else
    shift=mean(r0,1);
    r2=zeros(nu,nv);
    for i=1:nv
        r2(:,i)=r0(:,i)-shift(i);
    end
    rhoe0=sqrt(r2.*r2+z0.*z0);
    sinthe0=z0./rhoe0;
    costhe0=r2./rhoe0;
    drhoe0=rhoe0./(data.ns-1);
    dr=drhoe0.*costhe0;
    dz=drhoe0.*sinthe0;
    drds=dr;
    dzds=dz;
    %dr=(r-rax)/(data.ns-1);
    %dz=(z-zax)/(data.ns-1);
    %drds=dr;
    %dzds=dz;
    %dz(logical((z<0).*(dz>0)))=-dz(logical((z<0).*(dz>0)));
    %dz(logical((z>0).*(dz<0)))=-dz(logical((z>0).*(dz<0)));
    re0(:,:)=r0+dr*scale;
    ze0(:,:)=z0+dz*scale;
    re(:,:)=r+dr*scale;
    ze(:,:)=z+dz*scale;
end
% Make a quick plot
fig=figure;
if nv > 1
    subplot(1,3,1);
    plot(r0(:,1),z0(:,1),'k');
    hold on
    plot(re0(:,1),ze0(:,1),'r');
    plot(re0(:,1),ze0(:,1),'or');
    hold off
    axis square; axis equal;
    subplot(1,3,2);
    i=round(nv/4);
    plot(r0(:,i),z0(:,i),'k');
    hold on
    plot(re0(:,i),ze0(:,i),'r');
    plot(re0(:,i),ze0(:,i),'or');
    hold off
    axis square; axis equal;
    title('New Surface');
    subplot(1,3,3);
    i=round(nv/2);
    plot(r0(:,i),z0(:,i),'k');
    hold on
    plot(re0(:,i),ze0(:,i),'r');
    plot(re0(:,i),ze0(:,i),'or');
    hold off
    axis square; axis equal;
else
    plot(r0(:,1),z0(:,1),'k');
    hold on
    plot(re0(:,1),ze0(:,1),'r');
    plot(re0(:,1),ze0(:,1),'or');
    hold off
    axis square; axis equal;
end
pause(1.5);
close(fig);
% Convert to cartesian coords (xe,ye,z,jx,jy,jz)
xe=zeros(nu,nv);
ye=zeros(nu,nv);
xs=zeros(nu,nv);
ys=zeros(nu,nv);
%jx=zeros(nu,nv);
%jy=zeros(nu,nv);
for i=1:nv
    xe(:,i)=re(:,i).*cos(zeta(i));
    ye(:,i)=re(:,i).*sin(zeta(i));
    xs(:,i)=r(:,i).*cos(zeta(i));
    ys(:,i)=r(:,i).*sin(zeta(i));
%    jx(:,i)=jr(:,i).*cos(zeta(i))-jphi(:,i).*sin(zeta(i));
%    jy(:,i)=jr(:,i).*sin(zeta(i))+jphi(:,i).*cos(zeta(i));
end
% Calculate Coil Filed
bxe=zeros(nu,nv);
bye=zeros(nu,nv);
bze=zeros(nu,nv);
if ~novac
    disp('  - Calculating Coil Field');
    coildata=coil_biot_prep(coildata);
    for i=1:nu
        if i==1, tic; end
        for j=1:nv
            [bxe(i,j) bye(i,j) bze(i,j)]=coil_biot(coildata,xe(i,j),ye(i,j),ze(i,j),data.extcur);
        end
        if i==1; time=toc; disp(['      Est time to Calculate: ' num2str(time*(nu-1)) ' [s]']); end
    end
end
% Now extend to entire plasma (This doesn't work)
% if nv > 1
%     disp('  - Calculating Plasma Field');
%     phi=0:2*pi/(data.nfp):2*pi;
%     xe2=xe;
%     ye2=ye;
%     ze2=ze;
%     jx2=jx;
%     jy2=jy;
%     jz2=jz;
%     xs2=xs;
%     ys2=ys;
%     z2=z;
%     bxe2=bxe;
%     bye2=bye;
%     bze2=bze;
%     if (data.nfp > 1)
%         for i=2:data.nfp
%             xe2=[xe2 xe.*cos(phi(i))-ye.*sin(phi(i))];
%             ye2=[ye2 ye.*cos(phi(i))+xe.*sin(phi(i))];
%             ze2=[ze2 ze];
%             xs2=[xs2 xs.*cos(phi(i))-ys.*sin(phi(i))];
%             ys2=[ys2 ys.*cos(phi(i))+xs.*sin(phi(i))];
%             z2=[z2 z];
%             jx2=[jx2 jx.*cos(phi(i))-jy.*sin(phi(i))];
%             jy2=[jy2 jy.*cos(phi(i))+jx.*sin(phi(i))];
%             jz2=[jz2 jz];
%             bxe2=[bxe2 bxe.*cos(phi(i))-bye.*sin(phi(i))];
%             bye2=[bye2 bye.*cos(phi(i))+bxe.*sin(phi(i))];
%             bze2=[bze2 bze];
%         end
%     end
%     xe=xe2;
%     ye=ye2;
%     ze=ze2;
%     xs=xs2;
%     ys=ys2;
%     z=z2;
%     jx=jx2;
%     jy=jy2;
%     jz=jz2;
%     bxe=bxe2;
%     bye=bye2;
%     bze=bze2;
%     % Calculate B due to plasma (not axisymmetric)
%     fac=1/(4.*pi*nu*nv);
%     for i=1:nu
%         for j=1:nv
%             disinv=1./sqrt( (xe(i,j)-xs).*(xe(i,j)-xs)...
%                 +(ye(i,j)-ys).*(ye(i,j)-ys)...
%                 +(ze(i,j)-z).*(ze(i,j)-z));
%             dis3=disinv.*disinv.*disinv;
%             bxe(i,j)=bxe(i,j)+sum(sum(( (jy.*(ze(i,j)-z))...
%                 -(jz.*(ye(i,j)-ys))).*dis3))*fac;
%             bye(i,j)=bye(i,j)+sum(sum(( (jz.*(xe(i,j)-xs))...
%                 -(jx.*(ze(i,j)-z))).*dis3))*fac;
%             bze(i,j)=bze(i,j)+sum(sum(( (jx.*(ye(i,j)-ys))...
%                 -(jy.*(xe(i,j)-xs))).*dis3))*fac;
%         end
%     end
%     % Reducing Arrays and transforming to U,V
%     xe=xe(:,1:nv);
%     ye=ye(:,1:nv);
%     ze=ze(:,1:nv);
%     bxe=bxe(:,1:nv);
%     bye=bye(:,1:nv);
%     bze=bze(:,1:nv);
% end
bre=zeros(nu,nv);
bphie=zeros(nu,nv);
for i=1:nv
    bre(:,i)=bxe(:,i).*cos(zeta(i))+bye(:,i).*sin(zeta(i));
    bphie(:,i)=-bxe(:,i).*sin(zeta(i))+bye(:,i).*cos(zeta(i));
end
bve=bphie./re;
%%%%%%%%%%
% Initilize the Fourier amplitudes
disp('  - Calculating Boundary Fourier Transform');
refou=zeros(mpol+1,2*ntor+1);
zefou=zeros(mpol+1,2*ntor+1);
re0fou=zeros(mpol+1,2*ntor+1);
ze0fou=zeros(mpol+1,2*ntor+1);
bufou=zeros(mpol+1,2*ntor+1);
bvfou=zeros(mpol+1,2*ntor+1);
bsfou=zeros(mpol+1,2*ntor+1);
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
% Transform Position information
for m1=1:mpol+1
    for n1=1:2*ntor+1
        for i=1:nu
            for j=1:nv
                re0fou(m1,n1)=re0fou(m1,n1)+re0(i,j)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
                ze0fou(m1,n1)=ze0fou(m1,n1)+ze0(i,j)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
                refou(m1,n1)=refou(m1,n1)+re(i,j)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
                zefou(m1,n1)=zefou(m1,n1)+ze(i,j)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
            end
        end
    end
end
% Now get derivatives on the new surface
disp('  - Calculating Field Fourier Transform');
ruse=zeros(mpol+1,2*ntor+1);
rvse=zeros(mpol+1,2*ntor+1);
zuce=zeros(mpol+1,2*ntor+1);
zvce=zeros(mpol+1,2*ntor+1);
for i=1:mpol+1
    for j=1:2*ntor+1
        m=i-1;
        n=j-ntor-1;
        ruse(i,j)=-m.*refou(i,j);
        rvse(i,j)=-n.*refou(i,j);
        zuce(i,j)=m.*zefou(i,j);
        zvce(i,j)=n.*zefou(i,j);
    end
end
drdue=squeeze(sfunct(theta,zeta,ruse,data.nfp));
%drdve=squeeze(sfunct(theta,zeta,rvse,data.nfp));
dzdue=squeeze(cfunct(theta,zeta,zuce,data.nfp));
%dzdve=squeeze(cfunct(theta,zeta,zvce,data.nfp));
% Now calculate bu and bs
bue=zeros(nu,nv);
bse=zeros(nu,nv);
%bue=(bre-bze-bve.*(drdve-dzdve))./(drdue-dzdue);
for i=1:nv
    bue(:,i)=bxe(:,i).*drdue(:,i).*cos(zeta(i))+...
        bye(:,i).*drdue(:,i).*sin(zeta(i))+...
        bze(:,i).*dzdue(:,i);
    bse(:,i)=bxe(:,i).*drds(:,i).*cos(zeta(i))+...
        bye(:,i).*drds(:,i).*sin(zeta(i))+...
        bze(:,i).*dzds(:,i);
end
% Transform Field information
for m1=1:mpol+1
    for n1=1:2*ntor+1
        for i=1:nu
            for j=1:nv
                bsfou(m1,n1)=bsfou(m1,n1)+bse(i,j)*...
                    (sinv(j,n1)*cosu(i,m1)+cosv(j,n1)*sinu(i,m1))*fnuv(m1);
                bufou(m1,n1)=bufou(m1,n1)+bue(i,j)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
                bvfou(m1,n1)=bvfou(m1,n1)+bve(i,j)*...
                    (cosv(j,n1)*cosu(i,m1)-sinv(j,n1)*sinu(i,m1))*fnuv(m1);
            end
        end
    end
end
% Now we add these surfaces to the x and b array
for i=1:mpol+1
    for j=1:2*ntor+1
        j2=j+(piesn-ntor);
        x(data.ns+1,i,j2,1)=re0fou(i,j);
        x(data.ns+1,i,j2,2)=ze0fou(i,j);
        b(data.ns+1,i,j2,1)=bufou(i,j);
        b(data.ns+1,i,j2,2)=bvfou(i,j);
        b(data.ns+1,i,j2,3)=bsfou(i,j);
    end
end
disp('  - Symmetrizing Data');
% Do a symmetry check on the fields
for k=1:data.ns+1
        for j=1:piesn % Sweep the negative n
            j2=2*piesn+2-j;
            sum2=0.5*(x(k,1,j,1)+x(k,1,j2,1));
            x(k,1,j,1)=sum2;
            x(k,1,j2,1)=sum2;
            sum2=0.5*(x(k,1,j,2)+x(k,1,j2,2));
            x(k,1,j,2)=sum2;
            x(k,1,j2,2)=sum2;
            sum2=0.5*(b(k,1,j,1)+b(k,1,j2,1));
            b(k,1,j,1)=sum2;
            b(k,1,j2,1)=sum2;
            sum2=0.5*(b(k,1,j,2)+b(k,1,j2,2));
            b(k,1,j,2)=sum2;
            b(k,1,j2,2)=sum2;
            sum2=0.5*(b(k,1,j,3)+b(k,1,j2,3));
            b(k,1,j,3)=sum2;
            b(k,1,j2,3)=sum2;
        end
end
% Also remove (0,0) component of sine Coefficients
x(:,1,piesn+1,2)=0.0;  % Zs
b(:,1,piesn+1,3)=0.0;  % Bs
return
end

