function [rho,jrad] = beams3d_calc_jrad(beam_data)
%BEAMS3D_CALC_JRAD Calculates radial fast ion current
%   The BEAMS3D_CALC_JRAD routine calculates the radial profile of the fast
%   ion radial current <j_{rad}>.  This is done by differencing the inital 
%   radial and steady state radial distribution functions.  This must be 
%   multiplied by the minor radius of the device to get units of A/m^2.
%
%   Note this whole routine is experimental.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [rho, jrad] = beams3d_calc_jrad(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       0.1

% Defaults
rho=[]; jrad=[];
ec = 1.60217662E-19;

% Determine subset of particles for initial distribution
tdex=1; pdex=true(1,beam_data.nparticles);
if(double(beam_data.lbeam))
    tdex=3;
    pdex = beam_data.neut_lines(3,:)==0;
end

% Calc initial distribution
drho = beam_data.rho(2)-beam_data.rho(1);
drho2 = drho/2;
rho_lines = sqrt(beam_data.S_lines(tdex,pdex));
dist0 = zeros(1,beam_data.ns_prof1);
for i=1:beam_data.ns_prof1
    dex = and(rho_lines>=beam_data.rho(i)-drho2,rho_lines<beam_data.rho(i)+drho2);
    dist0(i) = sum(beam_data.Weight(dex));
end

% Put in volume units
[s,~,Vp]=beams3d_volume(beam_data);
rho_t = sqrt(s);
Vp = pchip(rho_t,Vp,beam_data.rho).*0.5.*rho_t;
dist0=dist0./Vp;

% Radial distribution function
distf=sum(beam_data.dense_prof,[1]);

% Calculate outputs.
rho  = beam_data.rho;
jrad = ec.*(dist0-distf).*rho;


end

