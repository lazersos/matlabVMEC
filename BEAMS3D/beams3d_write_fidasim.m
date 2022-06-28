function [f, denf] = beams3d_write_fidasim(data, name)
%BEAMS3D_WRITE FIDASIM produces the FIDASIM input files after a run of
%beams3d. For new versions of BEAMS3D, the same functionality can be
%achieved with the -fidasim flag. The distribution function and all
%quantities are output in the standard BEAMS3D cylindrical grid as in the
%input data structure. The energy range is set to 0-100keV for now. Flow
%velocities is set to 0 and electric field is calculated from gradient of POT_ARR.
%Example usage:
%   data = read_beams3d('test.h5');
%   beams3d_write_fidasim(data);

ec  = 1.60217662E-19; % electron charge [C]
amu = 1.66053906660E-27; % Dalton [kg]

filename_dist = ['fidasim_', name,'_distribution.h5'];
filename_eq = ['fidasim_', name,'_equilibrium.h5'];
time = 0.0;
mass = data.mass(1);

%Get grid points
% Return an RPZ shaped array;
raxis   = data.raxis;
zaxis   = data.zaxis;
paxis   = data.phiaxis;

Emax = 0.5.*mass.*data.partvmax.^2./ec;
Eaxis   = linspace(0,Emax,data.ns_prof4);%%0:10E3:100E3;
pitchaxis = linspace(-1,1,data.ns_prof5);%-1:0.1:1;
[R,P,Z,E,PITCH] = ndgrid(raxis,paxis,zaxis,Eaxis,pitchaxis);
%nsave = size(R);
ntotal = numel(R);
R = reshape(R,[1 ntotal]);
P = reshape(P,[1 ntotal]);
Z = reshape(Z,[1 ntotal]);
E = reshape(E,[1 ntotal]);
PITCH = reshape(PITCH,[1 ntotal]);


% V = sqrt(ec.*2.*E./mass);
% jac = 1.0./(mass.*sqrt(1-PITCH.*PITCH));
% jac = V ./ mass .* ec / 1000;
% f=squeeze(sum(data.dist_prof,1));
% f = f.*jac;
f=beams3d_getdistrpzEpitch(data,R,P,Z,E,PITCH);
f = sum(f,1);%.*ec/1000*1e6;%*1e6/2/pi/100; %keV and cm^-3;
f = reshape(f,[numel(raxis), numel(paxis), numel(zaxis), numel(Eaxis), numel(pitchaxis)]);
f= permute(f,[4, 5, 1, 3, 2]); %Align with FIDASIM axis order

%Quick fix
f(isnan(f))=0;
f(isinf(f))=0;

%Convert to cm
raxis   = data.raxis*100;
zaxis   = data.zaxis*100;
Eaxis = Eaxis/1000;

denf = squeeze(trapz(Eaxis, f,1));
denf = squeeze(trapz(pitchaxis, denf,1));%*100;

write_new_data_with_att(filename_dist,'/f',f,{'description', 'units'}, {'Distribution Function (nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida)', 'part/(cm^3 keV)'})
write_new_data_with_att(filename_dist,'/denf',denf,{'description', 'units'}, {'Fast-ion density (nr_fida,nz_fida,nphi_fida)','part/(cm^3)'})

write_new_data_with_att(filename_dist,'/nenergy',numel(Eaxis),'description','Number of energy values');
write_new_data_with_att(filename_dist,'/npitch',numel(pitchaxis), 'description','Number of pitch values');
write_new_data_with_att(filename_dist,'/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_dist,'/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_dist,'/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_dist,'/time',time,'description','Distribution time');
write_new_data_with_att(filename_dist,'/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_dist,'/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_dist,'/z',zaxis,{'description', 'units'},{'Z','cm'});
write_new_data_with_att(filename_dist,'/energy',Eaxis,{'description', 'units'},{'Energy array','keV'});
write_new_data_with_att(filename_dist,'/pitch',pitchaxis,{'description', 'units'},{'Pitch array','-'});
write_new_data_with_att(filename_dist,'/r2d',repmat(raxis,1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_dist,'/z2d',repmat(zaxis',numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});


% EQUILIBRIUM FILE
%Fields
write_new_data_with_att(filename_eq,'/fields/nenergy',numel(Eaxis),'description','Number of energy values');
write_new_data_with_att(filename_eq,'/fields/npitch',numel(pitchaxis), 'description','Number of pitch values');
write_new_data_with_att(filename_eq,'/fields/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_eq,'/fields/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_eq,'/fields/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_eq,'/fields/time',time,'description','Distribution time');

write_new_data_with_att(filename_eq,'/fields/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_eq,'/fields/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_eq,'/fields/z',zaxis,{'description', 'units'},{'Z','cm'});

write_new_data_with_att(filename_eq,'/fields/r2d',repmat(raxis,1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_eq,'/fields/z2d',repmat(zaxis',numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});

write_new_data_with_att(filename_eq,'/fields/mask',ones(numel(raxis),numel(zaxis)),'description','Boolean mask that indicates where the fields are well defined');

write_new_data_with_att(filename_eq,'/fields/br',permute(data.B_R,[1,3,2]),{'description', 'units'},{'Magnetic field in the r-direction: Br(r,z,phi)','T'});
write_new_data_with_att(filename_eq,'/fields/bt',permute(data.B_PHI,[1,3,2]),{'description', 'units'},{'Magnetic field in the toroidal phi-direction: Bt(r,z,phi)','T'});
write_new_data_with_att(filename_eq,'/fields/bz',permute(data.B_Z,[1,3,2]),{'description', 'units'},{'Magnetic field in the z-direction: Bz(r,z,phi)','T'});

write_new_data_with_att(filename_eq,'/fields/er',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Electric field in the r-direction: Er(r,z,phi)','V/m'});
write_new_data_with_att(filename_eq,'/fields/et',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Electric field in the toroidal phi-direction: Et(r,z,phi)','V/m'});
write_new_data_with_att(filename_eq,'/fields/ez',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Electric field in the toroidal phi-direction: Ez(r,z,phi)','V/m'});

%Plasma
write_new_data_with_att(filename_eq,'/plasma/nenergy',numel(Eaxis),'description','Number of energy values');
write_new_data_with_att(filename_eq,'/plasma/npitch',numel(pitchaxis), 'description','Number of pitch values');
write_new_data_with_att(filename_eq,'/plasma/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_eq,'/plasma/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_eq,'/plasma/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_eq,'/plasma/time',time,'description','Distribution time');

write_new_data_with_att(filename_eq,'/plasma/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_eq,'/plasma/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_eq,'/plasma/z',zaxis,{'description', 'units'},{'Z','cm'});

write_new_data_with_att(filename_eq,'/plasma/r2d',repmat(raxis,1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_eq,'/plasma/z2d',repmat(zaxis',numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});

write_new_data_with_att(filename_eq,'/plasma/mask',ones(numel(raxis),numel(zaxis)),'description','Boolean mask that indicates where the fields are well defined');

[er, et, ez] = gradient(data.POT_ARR,mean(diff(raxis)),mean(diff(paxis)),mean(diff(zaxis))); %first output corresponds to gradient along 2nd dimension???

% write_new_data_with_att(filename_eq,'/plasma/er',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the r-direction: Vr(r,z,phi)','cm/s'});
% write_new_data_with_att(filename_eq,'/plasma/et',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vt(r,z,phi)','cm/s'});
% write_new_data_with_att(filename_eq,'/plasma/ez',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vz(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/er',-er,{'description', 'units'},{'Bulk plasma flow in the r-direction: Vr(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/et',-et,{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vt(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/ez',-ez,{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vz(r,z,phi)','cm/s'});

write_new_data_with_att(filename_eq,'/plasma/dene',data.NE,{'description', 'units'},{'Electron Number Density: Dene(r,z, phi)','cm^-3'});
write_new_data_with_att(filename_eq,'/plasma/te',data.TE,{'description', 'units'},{'Electron Temperature: Te(r,z,phi)','keV'});
write_new_data_with_att(filename_eq,'/plasma/ti',data.TI,{'description', 'units'},{'Ion Temperature: Ti(r,z,phi)','keV'});
write_new_data_with_att(filename_eq,'/plasma/zeff',data.ZEFF_ARR,{'description', 'units'},{'Effective Nuclear Charge: Zeff(r,z,phi)','-'});
end

function write_new_data_with_att(filename,loc,data,attrs,vals)
h5create(filename,loc, size(data));
h5write(filename,loc,data);
if iscell(attrs)
for i = 1:numel(attrs)
    h5writeatt(filename,loc,attrs{i},vals{i})
end
else
    h5writeatt(filename,loc,attrs,vals)
end
end