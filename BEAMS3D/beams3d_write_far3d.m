function beams3d_write_far3d(beam_data,vmec_data)
%BEAMS3D_WRITE_FAR3D(beam_data,vmec_data) Writes FAR3D Input
%   The BEAMS3D_WRITE_FAR3D routine takes a beams3d and vmec data strcuture
%   an outputs a text file which the FAR3D code can use as input.  The
%   input data structures are those as returned by READ_BEAMS3D and
%   READ_VMEC.
%
%   Example usage
%      vmec_data = read_vmec('wout_test.nc');
%      beam_data = read_beams3d('beams3d_test.h5');
%      beams3d_write_far3d(beam_data,vmec_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.00

% Defaults
ec=1.60217663E-19;
out_filename = 'taefl_input_beams3d.txt';
Raxis = 5.55;
rho_out = linspace(0,1,101);
beam_dex = 1:beam_data.nbeams;

% Calculate values
s_out = rho_out.^2;
Rmajor = beams3d_calc_Rmajor(beam_data);
B0 = beams3d_b0(beam_data);
Aminor = beams3d_calc_aminor(beam_data);
[rho_dist, density, vllb, vperpb, pll, pperp, pcross] = beams3d_calc_moments(beam_data);
[raxis,~] = beams3d_magaxis(beam_data);
[s_prof, ne_prof, te_prof, ti_prof, zeff_prof] = beams3d_profiles(beam_data);

% Sume over beams
density = sum(density(beam_dex,:),1);
vllb = sum(vllb(beam_dex,:),1);
vperpb = sum(vperpb(beam_dex,:),1);
pll = sum(pll(beam_dex,:),1);
pperp = sum(pperp(beam_dex,:),1);
pcross = sum(pcross(beam_dex,:),1);

% Interpolate to single axis
s_vmec = linspace(0,1,vmec_data.ns);
q = pchip(s_vmec,1./vmec_data.iotaf,s_out);
nbeam = pchip(rho_dist,density,rho_out);
ni = pchip(s_prof,ne_prof./zeff_prof,s_out);
ne = pchip(s_prof,ne_prof,s_out);
nimp = rho_out.*0;
ti = pchip(s_prof,ti_prof,s_out);
te = pchip(s_prof,te_prof,s_out);
pbeam = pchip(rho_dist,sqrt(pll.^2+pperp.^2),rho_out);
tbeam = pbeam./(nbeam.*ec);
ptherm = ec.*(ne.*te+ni.*ti);
pequil = pchip(s_vmec,vmec_data.presf,s_out);
zeff = pchip(s_prof,zeff_prof,s_out);
torrot = rho_out.*0;

% Make output array
den_fact = 1E-19; % 10^13 cm^-3
temp_fact = 1E-3; % keV
pres_fact = 1E-3; % kPa
out_array=[rho_out; q; max(nbeam.*den_fact,0); ni.*den_fact; ne.*den_fact; nimp.*den_fact; max(tbeam.*temp_fact,0); ti.*temp_fact; te.*temp_fact; pbeam.*pres_fact; ptherm.*pres_fact; pequil.*pres_fact; zeff; torrot; torrot];

%Write file
fid=fopen(out_filename,'w');
fprintf(fid,'PLASMA GEOMETRY \n');
fprintf(fid,'Vacuum Toroidal magnetic field at R=%7.5fm [Tesla]\n',raxis(1));
fprintf(fid,'    %7.5f\n',B0);
fprintf(fid,'Geometric Center Major radius [m]\n');
fprintf(fid,'    %7.5f\n',Rmajor(1));
fprintf(fid,'Minor radius [m]\n');
fprintf(fid,'    %7.5f\n',Aminor);
fprintf(fid,'Avg. Elongation\n');
fprintf(fid,'    %7.5f\n',1.0);
fprintf(fid,'Avg. Top/Bottom Triangularity\n');
fprintf(fid,'    %7.5f\n',1.0);
fprintf(fid,'Main Contaminant Species\n');
fprintf(fid,'    12C\n');
fprintf(fid,'Main Ion Species mass/proton mass\n');
fprintf(fid,'    %7.5f\n',1.0);
fprintf(fid,'TRYING TO GET TO BETA(0)=%7.5f , Rmax=%7.5f\n',vmec_data.betaxis,vmec_data.rmax_surf);
fprintf(fid,'\nRho(norml. sqrt. toroid. flux), q, Beam Ion Density(10^13 cm^-3), Ion Density(10^13 cm^-3), Elec Density(10^13 cm^-3), Impurity Density(10^13 cm^-3), Beam Ion Effective Temp(keV), Ion Temp(keV), Electron Temp(keV), Beam Pressure(kPa), Thermal Pressure(kPa), Equil.Pressure(kPa), Zeff, Tor Rot(kHz), Tor Rot(10^5 m/s)');
fprintf(fid,'%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f \n',out_array);

fclose(fid);

end