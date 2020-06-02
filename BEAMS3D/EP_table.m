function EP_table(energy,mass,field,ratio)
%EP_table(energy,mass,ratio) Makes a EP table for BEAMS3D
%This function takes a particle total energy (in eV), it's mass in kg,
%a reference magnetic field (in T), and an array of parallel to 
%perpendicular particle velocities.  It then prints to the screen a table
%of velocities and magnetic moments.  It also outputs the necessary 
%namelist parameters for the energetic particle module of STELLOPT.
%   Energy [ev]
%   Mass   [kg]
%   Field  [T]
%   ratio  array (vll/vperp)
%
% Example usage
%      ratio=[0.1 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25];
%      EP_table(60E3,1.6726231E-27,1.5,ratio); 60 [keV] protons @ 1.5 [T]
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00


ec=1.60217733E-19;

energy_j =energy*ec;

vtot = sqrt(2.*energy_j./mass);
vll = vtot./sqrt(1.0+(1.0./ratio).^2);
vperp = vll./ratio;
energy_perp = 0.5.*mass.*vperp.*vperp;
mu = energy_perp./field;

disp('  #     Vpara/Vperp              Vpara                Vperp               Mu');
for i=1:length(ratio)
    disp(num2str([i ratio(i) vll(i) vperp(i) mu(i)],'  %3d       %5f  %20.10e %20.10e %20.10e'));
end
disp(['!  Total Energy: ' num2str(energy./1000,'%5.2f') ' [keV]']);
disp(['!        Ratios: ' num2str(ratio,'%5.2f') ' [Vpara/Vperp]']);
disp(['  NP_ORBIT = ' num2str(length(ratio),'%4d')]);
disp(['  MASS_ORBIT = ' num2str(mass,'%20.10E')]);
disp(['  Z_ORBIT = 1.0']);
disp(['  VLL_ORBIT = ' num2str(vll,' %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E\n')]);
disp(['  VPERP_ORBIT = ' num2str(vperp,' %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E\n')]);
disp(['  MU_ORBIT  = ' num2str(mu,' %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E\n')]);


end

