function fieldlines_input(input_file,vmec_file)
%FIELDILNES_INPUT Generates the FIELDLINES input namelist.
%   Detailed explanation goes here

% Set Some Defaults
npass=100;
ntheta=36;
nzeta=36;
poinc_res=0.1;

% Open files
out_id=fopen(input_file,'a+'); %Open input_file for appending
vmec_data=read_vmec(vmec_file);

% Transform the output to real space
theta=0:2*pi/(ntheta-1):2*pi;
zeta=0:2*pi/(ntheta-1):2*pi;
r=cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
z=sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
if (vmec_data.iasym == 1)
    r=r+sfunct(theta,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
    z=z+cfunct(theta,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
end

% Setup the output arrays with radial grid.
r0=r(:,1,1);
z0=z(:,1,1);
phi0=r0.*0.0;

% Now add the surface gridpoints
%r0=[r0 reshape(r(vmec_data.ns,:,:),1,ntheta*nzeta)];
%z0=[z0 reshape(z(vmec_data.ns,:,:),1,ntheta*nzeta)];

fprintf(out_id,'&FIELDLINES_INPUT\n');
write_namelist_int(out_id,'NR',201);
write_namelist_int(out_id,'NZ',201);
write_namelist_int(out_id,'NPHI',vmec_data.ntor*4);
write_namelist_flt(out_id,'RMIN',0.0);
write_namelist_flt(out_id,'RMAX',10.0);
write_namelist_flt(out_id,'ZMIN',-10.0);
write_namelist_flt(out_id,'ZMAX',10.0);
write_namelist_flt(out_id,'PHIMIN',0.0);
write_namelist_flt(out_id,'PHIMAX',2*pi/vmec_data.nfp);
write_namelist_flt(out_id,'MU',0.0);
write_namelist_vec(out_id,'R_START',r0);
write_namelist_vec(out_id,'Z_START',z0);
write_namelist_vec(out_id,'PHI_START',phi0);
write_namelist_vec(out_id,'PHI_END',phi0+npass*2*pi/vmec_data.nfp);
write_namelist_int(out_id,'NPOINC',360/(vmec_data.nfp*poinc_res));
fprintf(out_id,' INT_TYPE = ''NAG'' \n');
write_namelist_flt(out_id,'FOLLOW_TOL',1.0E-9);
write_namelist_flt(out_id,'VC_ADAPT_TOL',1.0E-4);
fprintf(out_id,'&END\n');
fclose(out_id);

end

