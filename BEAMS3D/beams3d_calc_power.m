function [rho, Qe, Qi] = beams3d_calc_power(beam_data,vmec_data)
%BEAMS3D_CALC_POWER Calcluate power deposition from distribution
%   The BEAMS3D_CALC_POWER function calculats the power depoistion from the
%   5D BEAMS3D distribution function.  It takes a BEAMS3D data structure
%   (as read from READ_BEAMS3D) and a VMEC data strucutre (as read from
%   READ_VMEC) as inputs.  The radial grid (rho), power depositied to the
%   electrons (Qe [W/m^3]), and power depositied to the ions (Qi [W/m^3])
%   are returned.  The Qe and Qi variables are [nbeams,nrho] sized.
%
%   Example usage
%      vmec_data = read_vmec('wout_test.nc');
%      beam_data = read_beams3d('beams3d_test.h5');
%      [rho,Qe,Qi]=beams3d_calc_power(beam_data,vmec_data);
%
%   See also read_beams3d, read_vmec, beams3d_calc_pll.
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.00

%Defaults
ec = 1.60217662E-19;
me = 9.10938356E-31;
rho=[];
Qe=[];
Qi=[];

% Check format of inputs
if ~isfield(beam_data,'dist_prof')
    disp('ERROR: No 5D distribution profile in beam_data!');
    return;
end
if ~isfield(beam_data,'plasma_mass')
    disp('ERROR: No plasma_mass variable found in beam_data!');
    return;
end
if ~isfield(vmec_data,'vp')
    disp('ERROR: No Vp volume array found in vmec_data!');
    return;
end
if isfield(beam_data,'dist_prof_fix')
    disp('ERROR: dist_prof is normalized to 5D volume, this isnt tested');
    return;
end

% setup grid
nrho = double(beam_data.ns_prof1);
edges = 0:1/nrho:1;
rho_beam = (edges(1:end-1)+edges(2:end)).*0.5;
rho_vmec = sqrt(0:1/(vmec_data.ns-1):1);
npol = double(beam_data.ns_prof2);
ntor = double(beam_data.ns_prof3);
nvpa = double(beam_data.ns_prof4);
nvpe = double(beam_data.ns_prof5);
edges = 0:1/npol:1;
u_beam = (edges(1:end-1)+edges(2:end)).*pi;
edges = 0:1/ntor:1;
v_beam = (edges(1:end-1)+edges(2:end)).*pi;
nfp_beam = round(2*pi./beam_data.phiaxis(end));
v_beam = v_beam./nfp_beam;
edges = -1:2/nvpa:1;
vpa_beam = 0.5.*(edges(1:end-1)+edges(2:end)).*beam_data.partvmax;
edges = 0:1/nvpe:1;
vpe_beam = 0.5.*(edges(1:end-1)+edges(2:end)).*beam_data.partvmax;

% Create grided data
[u_grid,v_grid]=ndgrid(u_beam,v_beam);
[vpa_grid,vpe_grid]=ndgrid(vpa_beam,vpe_beam);
vto_grid = sqrt(vpa_grid.^2+vpe_grid.^2);

% Now make an array of R and Z
u_lin=reshape(u_grid,[1 numel(u_grid)]);
v_lin=reshape(v_grid,[1 numel(v_grid)]);
R_lin=zeros([length(rho_beam) length(u_lin)]);
Z_lin=zeros([length(rho_beam) length(u_lin)]);
v_hel=zeros([length(rho_beam) length(u_lin)]);
for i=1:length(u_lin)
    R_temp = cfunct(u_lin(i),v_lin(i),vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    Z_temp = sfunct(u_lin(i),v_lin(i),vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    R_lin(:,i) = pchip(rho_vmec,R_temp,rho_beam);
    Z_lin(:,i) = pchip(rho_vmec,Z_temp,rho_beam);
    v_hel(:,i) = v_lin(i);
end
R_lin = reshape(R_lin,[1 numel(R_lin)]);
Z_lin = reshape(Z_lin,[1 numel(Z_lin)]);
v_hel = reshape(v_hel,[1 numel(v_hel)]);

% Calculate the background
TE   = permute(beam_data.TE,[2 1 3]);
NE   = permute(beam_data.NE,[2 1 3]);
%TI   = permute(beam_data.TI,[2 1 3]);
%ZE   = permute(beam_data.ZEFF_ARR,[2 1 3]);
TE_lin = interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    TE,R_lin,v_hel,Z_lin,'cubic');
%TI_lin = interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
%    TI,R_lin,v_hel,Z_lin);
NE_lin = interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    NE,R_lin,v_hel,Z_lin,'cubic');
%ZE_lin = interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
%    ZE,R_lin,v_hel,Z_lin);

% specify some beam parameters as a function of beam
nbeams = double(beam_data.nbeams);
mass   = zeros(1,nbeams);
myz    = zeros(1,nbeams);
for i=1:nbeams
    dex = find(beam_data.Beam==i,1,'first');
    mass(i) = beam_data.mass(dex);
    myz(i)  = beam_data.Zatom(dex);
end

% Since mass and Z can vary we need to divide the problem up by beamline
% After here we should have an [nbeam nrho npol ntor] arrays of C1, C2
TE3 = TE_lin.^3;
C1 = zeros(nbeams,numel(NE_lin));
C2 = zeros(nbeams,numel(NE_lin));
for i=1:nbeams
    coulomb_log=zeros([1 length(TE3)]);
    dex = TE_lin < 10.*myz(i).*myz(i);
    coulomb_log(dex) = 23 - log(myz(i).*sqrt(NE_lin(dex).*1E-6./TE3(dex)));
    dex = ~dex;
    coulomb_log(dex) = 24 - log(sqrt(NE_lin(dex).*1E-6)./TE_lin(dex));
    coulomb_log(coulomb_log <= 1) = 1;
    v_crit = sqrt(2.*ec.*TE_lin./beam_data.plasma_mass).*(0.75.*sqrt(pi.*beam_data.plasma_mass./me)).^(1.0/3.0);
    vcrit_cube = v_crit.^3;
    tau_spit = 3.777183E41.*mass(i).*sqrt(TE3)./(NE_lin.*myz(i).*myz(i).*coulomb_log);
    C1(i,:) = 1./tau_spit;
    C2(i,:) = vcrit_cube./tau_spit;
end
C1 = repmat(reshape(C1,[nbeams nrho npol ntor]),[1 1 1 1 nvpa nvpe]);
C2 = repmat(reshape(C2,[nbeams nrho npol ntor]),[1 1 1 1 nvpa nvpe]);

% Calc the change in velocity
dve = zeros(nbeams,nrho,npol,ntor,nvpa,nvpe);
dvi = zeros(nbeams,nrho,npol,ntor,nvpa,nvpe);
vt2_grid = vto_grid.^2;
for i=1:nvpa
    for j=1:nvpe
        dve(:,:,:,:,i,j) = C1(:,:,:,:,i,j).*vto_grid(i,j);
        dvi(:,:,:,:,i,j) = C2(:,:,:,:,i,j)./vt2_grid(i,j);
    end
end

% Calc Q
Ee = zeros(nbeams,nrho,npol,ntor,nvpa,nvpe);
Ei = zeros(nbeams,nrho,npol,ntor,nvpa,nvpe);
for i=1:nvpa
    for j=1:nvpe
        Ee(:,:,:,:,i,j) = dve(:,:,:,:,i,j).*vto_grid(i,j);
        Ei(:,:,:,:,i,j) = dvi(:,:,:,:,i,j).*vto_grid(i,j);
    end
end
Qe = squeeze(sum(Ee.*beam_data.dist_prof,[3 4 5 6]));
Qi = squeeze(sum(Ei.*beam_data.dist_prof,[3 4 5 6]));
for i=1:nbeams
    Qe(i,:) = Qe(i,:) .*mass(i);
    Qi(i,:) = Qi(i,:) .*mass(i);
end

% ds to dV
vp_beam = repmat(2.*rho_beam.*pchip(rho_vmec,vmec_data.vp.*4.*pi.*pi,rho_beam),[6 1])./nrho;
Qe = Qe./vp_beam;
Qi = Qi./vp_beam;
rho = rho_beam;

end

