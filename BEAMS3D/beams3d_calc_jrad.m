function [rho,Aminor, jrad] = beams3d_calc_jrad(beam_data)
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
%       [rho, Aminor, jrad] = beams3d_calc_jrad(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       0.1

% Defaults
rho=[]; jrad=[];
ec = 1.60217662E-19;

% Determine subset of particles for initial distribution
tdex=1;
if(double(beam_data.lbeam))
    tdex=3;
end

% Get volume
rho = beam_data.rho;
[s,V,Vp]=beams3d_volume(beam_data);
dVdrho=pchip([0 sqrt(s)],[0 Vp.*2.*sqrt(s)],rho);
Vrho=pchip([0 sqrt(s)],[0 V],rho);

% Calculate thermalized current
Itherm=zeros(beam_data.nbeams,beam_data.ns_prof1);
Ilost=zeros(beam_data.nbeams,beam_data.ns_prof1);
for k = 1: beam_data.nbeams
    thermdex = find(and(beam_data.end_state==1,beam_data.Beam==k)); % just thermalized particles
    lostdex = find(and(beam_data.end_state==2,beam_data.Beam==k)); % just thermalized particles
    rho_edge = 0:1./double(beam_data.ns_prof1):1;
    for i=thermdex'
        slines = beam_data.S_lines(:,i)';
        s1 = slines(tdex);
        j1 = find(beam_data.R_lines(:,i)>0,1,'last');
        s2 = slines(j1);
        j1 = max(sum(rho_edge>s1),1);
        j2 = min(sum(rho_edge>s2),beam_data.ns_prof1);
        Itherm(k,j1:j2) = Itherm(k,j1:j2) + (beam_data.Weight(i).*beam_data.charge(i).*sign(j2-j1));
    end
    
    
    % Calculat the lost current
    for i=lostdex'
        slines = beam_data.S_lines(:,i)';
        s1 = slines(tdex);
        j1 = max(sum(rho_edge>s1),1);
        Ilost(k,j1:end) = Ilost(k,j1:end) + (beam_data.Weight(i).*beam_data.charge(i));
    end
end

% Calc Aminor
Aminor = beams3d_calc_aminor(beam_data);

%V = pi*r*r*L;
%A = 2*pi*r*L;
%  = 2*V/r

% Now calc jrad
area = repmat(2.*Vrho./(Aminor.*rho),[beam_data.nbeams,1]);
jrad = (Ilost+Itherm)./area;


end

