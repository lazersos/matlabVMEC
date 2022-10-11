function [f, denf] = beams3d_write_fidasim(data, name,varargin)
%BEAMS3D_WRITE FIDASIM produces the FIDASIM input files after a run of
%beams3d. For new versions of BEAMS3D, the same functionality can be
%achieved with the -fidasim flag. The distribution function and all
%quantities are output in the standard BEAMS3D cylindrical grid as standard.
% Alternatively, 'n', nr, nphi, nz, nE, np can be passed as input to change the
% grid spacing, or 'axis',raxis,phiaxis,zaxis can be passed for setting the
% axis precisely. The energy range is set to 0-Emax from the maximum particle
% velocity. Flow velocities are set to 0 and electric field is calculated
% from the gradient of POT_ARR. Optionally, the function outputs the fast
% ion distribution and density as variables.
%
%Example usage:
%   data = read_beams3d('test.h5');
%   beams3d_write_fidasim(data);

ec  = 1.60217662E-19; % electron charge [C]
amu = 1.66053906660E-27; % Dalton [kg]
lmovie = 0;

filename_dist = ['fidasim_', name,'_distribution.h5'];
filename_eq = ['fidasim_', name,'_equilibrium.h5'];
time = 0.0;
mass = data.mass(1);

if isfile(filename_dist)
    disp('Found old distribution file. Overwriting...')
    delete(filename_dist);
end
if isfile(filename_eq)
    disp('Found old equilibrium file. Overwriting...')
    delete(filename_eq);
end

%Get grid points
% Return an RPZ shaped array;
raxis   = data.raxis';%linspace(min(data.raxis), max(data.raxis), 70);
zaxis   = data.zaxis';%linspace(min(data.zaxis), max(data.zaxis), 71);
paxis   = data.phiaxis';%linspace(0, 2*pi, 15);
nE = data.ns_prof4;
np = data.ns_prof5;
% Handle varargin
if ~isempty(varargin)
    i=1;
    while i <= length(varargin)
        switch varargin{i}
            case{'axis'}
                i = i + 1;
                raxis = varargin{i};
                i = i + 1;
                paxis = varargin{i};
                i = i + 1;
                zaxis = varargin{i};
            case{'n'}
                i = i + 1;
                raxis   = linspace(min(data.raxis), max(data.raxis), varargin{i});
                i = i + 1;
                paxis   = linspace(0, 2*pi, varargin{i});
                i = i + 1;
                zaxis   = linspace(min(data.zaxis), max(data.zaxis), varargin{i});
                i = i + 1;
                nE = varargin{i};
                i = i + 1;
                np = varargin{i};
            case{'movie'}
                lmovie = 1;
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end

Emax = ceil(0.5.*mass.*data.partvmax.^2./ec/1e4)*1e4;
Eaxis   = 0.5*max(Emax)/nE:max(Emax)/nE:max(Emax);%linspace(0,Emax,nE);%%0:10E3:100E3;
pitchaxis = -1+1/np:2/np:1;%linspace(-1,1,np);%-1:0.1:1;
[R,P,Z,E,PITCH] = ndgrid(raxis,paxis,zaxis,Eaxis,pitchaxis);
[R3,P3,Z3] = ndgrid(raxis,paxis,zaxis);
[Rd,Pd,Zd] = ndgrid(data.raxis,data.phiaxis,data.zaxis);
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
raxis   = raxis.*100;
%R3 = R3 .*100;
zaxis   = zaxis.*100;
%Z3 = Z3 .*100;
Eaxis = Eaxis./1000;
%E = E/1000;

denf = squeeze(trapz(Eaxis, f,1));
denf = squeeze(trapz(pitchaxis, denf,1));%Integration in velocity space


F = griddedInterpolant(Rd,Pd,Zd, data.B_R, 'spline');
br = F(R3, mod(P3,max(data.phiaxis)), Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.B_PHI, 'spline');
bt = F(R3, mod(P3,max(data.phiaxis)), Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.B_Z, 'spline');
bz = F(R3, mod(P3,max(data.phiaxis)), Z3);

vr = zeros(size(br));
vt = vr;
vz = vr;

[er, et, ez] = gradient(data.POT_ARR,mean(diff(raxis)),mean(diff(paxis)),mean(diff(zaxis))); %first output corresponds to gradient along 2nd dimension???
F = griddedInterpolant(Rd,Pd,Zd, er, 'spline');
er = F(R3, mod(P3,max(data.phiaxis)), Z3);
F = griddedInterpolant(Rd,Pd,Zd, et, 'spline');
et = F(R3, mod(P3,max(data.phiaxis)), Z3);
F = griddedInterpolant(Rd,Pd,Zd, ez, 'spline');
ez = F(R3, mod(P3,max(data.phiaxis)), Z3);

F = griddedInterpolant(Rd,Pd,Zd, data.NE/1e6, 'spline');
ne = F(R3, mod(P3,max(data.phiaxis)), Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.TE/1000, 'spline');
te = F(R3, mod(P3,max(data.phiaxis)), Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.TI/1000, 'spline');
ti = F(R3, mod(P3,max(data.phiaxis)), Z3);
F = griddedInterpolant(Rd,Pd,Zd, data.ZEFF_ARR, 'spline');
zeff = F(R3, mod(P3,max(data.phiaxis)), Z3);
%denn = zeros(size(zeff));


write_new_data_with_att(filename_dist,'/type',1,'description', 'Distribution type: 1="Guiding Center Density Function", 2="Guiding Center Monte Carlo", 3="Full Orbit Monte Carlo"')
write_new_data_with_att(filename_dist,'/f',f,{'description', 'units'}, {'Distribution Function (nenergy_fida,npitch_fida,nr_fida,nz_fida,nphi_fida)', 'part/(cm^3 keV)'})
write_new_data_with_att(filename_dist,'/denf',denf,{'description', 'units'}, {'Fast-ion density (nr_fida,nz_fida,nphi_fida)','part/(cm^3)'})

write_new_data_with_att(filename_dist,'/nenergy',numel(Eaxis),'description','Number of energy values');
write_new_data_with_att(filename_dist,'/npitch',numel(pitchaxis), 'description','Number of pitch values');
write_new_data_with_att(filename_dist,'/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_dist,'/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_dist,'/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_dist,'/time',time,{'description', 'units'},{'Distribution time', 's'});
write_new_data_with_att(filename_dist,'/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_dist,'/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_dist,'/z',zaxis,{'description', 'units'},{'Z','cm'});
write_new_data_with_att(filename_dist,'/energy',Eaxis,{'description', 'units'},{'Energy array','keV'});
write_new_data_with_att(filename_dist,'/pitch',pitchaxis,{'description', 'units'},{'Pitch array','-'});
write_new_data_with_att(filename_dist,'/r2d',repmat(raxis',1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_dist,'/z2d',repmat(zaxis,numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});


% EQUILIBRIUM FILE
%Fields
write_new_data_with_att(filename_eq,'/fields/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_eq,'/fields/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_eq,'/fields/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_eq,'/fields/time',time,'description','Distribution time');

write_new_data_with_att(filename_eq,'/fields/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_eq,'/fields/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_eq,'/fields/z',zaxis,{'description', 'units'},{'Z','cm'});

write_new_data_with_att(filename_eq,'/fields/r2d',repmat(raxis',1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_eq,'/fields/z2d',repmat(zaxis,numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});

write_new_data_with_att(filename_eq,'/fields/mask',int32(ones(numel(raxis),numel(zaxis), numel(paxis))),'description','Boolean mask that indicates where the fields are well defined');



write_new_data_with_att(filename_eq,'/fields/br',permute(br,[1,3,2]),{'description', 'units'},{'Magnetic field in the r-direction: Br(r,z,phi)','T'});
write_new_data_with_att(filename_eq,'/fields/bt',permute(bt,[1,3,2]),{'description', 'units'},{'Magnetic field in the toroidal phi-direction: Bt(r,z,phi)','T'});
write_new_data_with_att(filename_eq,'/fields/bz',permute(bz,[1,3,2]),{'description', 'units'},{'Magnetic field in the z-direction: Bz(r,z,phi)','T'});


write_new_data_with_att(filename_eq,'/fields/er',permute(-er,[1,3,2]),{'description', 'units'},{'Electric field in the r-direction: Er(r,z,phi)','V/m'});
write_new_data_with_att(filename_eq,'/fields/et',permute(-et,[1,3,2]),{'description', 'units'},{'Electric field in the toroidal phi-direction: Et(r,z,phi)','V/m'});
write_new_data_with_att(filename_eq,'/fields/ez',permute(-ez,[1,3,2]),{'description', 'units'},{'Electric field in the toroidal phi-direction: Ez(r,z,phi)','V/m'});

%Plasma
write_new_data_with_att(filename_eq,'/plasma/nr',numel(raxis),'description','Number of R values');
write_new_data_with_att(filename_eq,'/plasma/nphi',numel(paxis),'description','Number of Phi values');
write_new_data_with_att(filename_eq,'/plasma/nz', numel(zaxis),'description','Number of Z values');
write_new_data_with_att(filename_eq,'/plasma/time',time,'description','Distribution time');

write_new_data_with_att(filename_eq,'/plasma/r',raxis,{'description', 'units'},{'Radius','cm'});
write_new_data_with_att(filename_eq,'/plasma/phi',paxis,{'description', 'units'},{'Toroidal angle','rad'});
write_new_data_with_att(filename_eq,'/plasma/z',zaxis,{'description', 'units'},{'Z','cm'});

write_new_data_with_att(filename_eq,'/plasma/r2d',repmat(raxis',1,numel(zaxis)),{'description', 'units'},{'R grid: R(r,z)','cm'});
write_new_data_with_att(filename_eq,'/plasma/z2d',repmat(zaxis,numel(raxis),1),{'description', 'units'},{'Z grid: Z(r,z)','cm'});

write_new_data_with_att(filename_eq,'/plasma/mask',int32(ones(numel(raxis),numel(zaxis), numel(paxis))),'description','Boolean mask that indicates where the fields are well defined');


% write_new_data_with_att(filename_eq,'/plasma/er',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the r-direction: Vr(r,z,phi)','cm/s'});
% write_new_data_with_att(filename_eq,'/plasma/et',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vt(r,z,phi)','cm/s'});
% write_new_data_with_att(filename_eq,'/plasma/ez',zeros(numel(raxis),numel(zaxis)),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vz(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/vr',permute(vr,[1,3,2]),{'description', 'units'},{'Bulk plasma flow in the r-direction: Vr(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/vt',permute(vt,[1,3,2]),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vt(r,z,phi)','cm/s'});
write_new_data_with_att(filename_eq,'/plasma/vz',permute(vz,[1,3,2]),{'description', 'units'},{'Bulk plasma flow in the toroidal phi-direction: Vz(r,z,phi)','cm/s'});



write_new_data_with_att(filename_eq,'/plasma/dene',permute(ne,[1,3,2]),{'description', 'units'},{'Electron Number Density: Dene(r,z, phi)','cm^-3'});
write_new_data_with_att(filename_eq,'/plasma/te',permute(te,[1,3,2]),{'description', 'units'},{'Electron Temperature: Te(r,z,phi)','keV'});
write_new_data_with_att(filename_eq,'/plasma/ti',permute(ti,[1,3,2]),{'description', 'units'},{'Ion Temperature: Ti(r,z,phi)','keV'});
write_new_data_with_att(filename_eq,'/plasma/zeff',permute(zeff,[1,3,2]),{'description', 'units'},{'Effective Nuclear Charge: Zeff(r,z,phi)','-'});
%write_new_data_with_att(filename_eq,'/plasma/denn',permute(denn,[1,3,2]),{'description', 'units'},{'Neutral density: denn(r,z,phi)','cm^-3'});

if lmovie
    v= VideoWriter(['output_' name '.avi']);
    open(v);
    f = figure;
    for i =1:numel(paxis)
        subplot(3,3,1)
        pixplot(raxis,zaxis,squeeze(br(:,i,:)))
        title('BR')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,2)
        pixplot(raxis,zaxis,squeeze(bt(:,i,:)))
        title('BT')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,3)
        pixplot(raxis,zaxis,squeeze(bz(:,i,:)))
        title('BZ')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,4)
        pixplot(raxis,zaxis,squeeze(er(:,i,:)))
        title('ER')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,5)
        pixplot(raxis,zaxis,squeeze(et(:,i,:)))
        title('ET')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,6)
        pixplot(raxis,zaxis,squeeze(ez(:,i,:)))
        title('EZ')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,7)
        pixplot(raxis,zaxis,squeeze(ne(:,i,:)))
        title('NE')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,8)
        pixplot(raxis,zaxis,squeeze(te(:,i,:)))
        title('TE')
        xlabel('R [m]')
        ylabel('Z [m]')
        subplot(3,3,9)
        pixplot(raxis,zaxis,squeeze(ti(:,i,:)))
        title('TI')
        xlabel('R [m]')
        ylabel('Z [m]')

        frame = getframe(f);
        writeVideo(v,frame);
    end
end

end

function write_new_data_with_att(filename,loc,data,attrs,vals)
if sum(size(data)==1)==1
    h5create(filename,loc, length(data)); %Initialize only one dimension for vectors!
else
    h5create(filename,loc, size(data));
end
h5write(filename,loc,squeeze(data));
if iscell(attrs)
    for i = 1:numel(attrs)
        h5writeatt(filename,loc,attrs{i},vals{i})
    end
else
    h5writeatt(filename,loc,attrs,vals)
end
end