function P = beams3d_powertable(beam_data,varargin)
%BEASM3D_POWERTABLE Calculate total power per beam
%   This subroutine calculates the total power born, going to the port,
%   shining through, deposited in the plasma, and lost to the wall.
%
% Example usage
%      beam_data=read_beams3d('beams3d_test.h5');  % Reads BEAMS3D HDF5 file
%       beams3d_powertable(beam_data,vmec_data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

vmec_data=[];
lverb=1;

if ~isempty(varargin)
    i=1;
    while i<=numel(varargin)
        if isstruct(varargin{i})
            if isfield(varargin{i},'datatype')
                switch varargin{i}.datatype
                    case 'wout'
                        vmec_data=varargin{i};
                end
            end
        else
            switch varargin{i}
                case 'quiet'
                    lverb=0;
            end
        end
        i=i+1;
    end
end

% Handle the calculation of dV/drho from VMEC dV/ds
if ~isempty(vmec_data)
    h_vmec = 1./(vmec_data.ns-1);
    s_vmec = 0:h_vmec:1;
    vp_vmec = vmec_data.vp.*4.*pi.*pi;
else
    [s_vmec, ~, vp_vmec] = beams3d_volume(beam_data);
end
rho_vmec = sqrt(s_vmec);
dVdrho_spl = pchip(rho_vmec,2.*rho_vmec.*vp_vmec);
ns = double(beam_data.ns_prof1);
h = 1./ns;
rho=0:h:1.0;
rho2=rho(1:end-1)+h/2;
dVdrho = ppval(dVdrho_spl,rho2);

nfp = round(2*pi./beam_data.phiaxis(end));

% Calc Vperp
vperp = beams3d_calc_vperp(beam_data);

% Calc Power
vll = beam_data.vll_lines; % Initial energy
v2  = vll.*vll+vperp.*vperp;
P   = zeros(beam_data.npoinc+1,beam_data.nparticles);
for i = 1:beam_data.npoinc+1
    P(i,:) = 0.5.*beam_data.mass'.*v2(i,:).*beam_data.Weight';
end

% Calc wall area
coords = beam_data.wall_vertex;
faces  = beam_data.wall_faces;
v0 = coords(faces(:,1),:);
v1 = coords(faces(:,2),:);
v2 = coords(faces(:,3),:);
d01 = v1-v0;
d02 = v2-v0;
FN  = cross(d01,d02);
A = 0.5*sqrt(sum(FN.*FN,2));

% Shine dexes
port_dex=beam_data.end_state==4;
shine_dex=beam_data.end_state==3;
wall_dex=beam_data.end_state==2;
therm_dex=beam_data.end_state==1;
orbit_dex=beam_data.end_state==0;

% Handle orbiting particles
P_orbit = zeros(1,beam_data.nparticles);
if any(orbit_dex>0)
    if lverb, disp('Circulating Particles found'); end
    orbit_index=beams3d_finddex(beam_data,'orbit_last');
    for i=1:beam_data.nparticles
        if orbit_index(i)==0, continue; end
        P_orbit(i) = P(orbit_index(i),i);
    end
end

% Handle thermalized particles
P_therm = zeros(1,beam_data.nparticles);
if any(therm_dex>0)
    if lverb, disp('Thermalized Particles found'); end
    therm_index=beams3d_finddex(beam_data,'therm_end');
    for i=1:beam_data.nparticles
        if therm_index(i)==0, continue; end
        P_therm(i) = P(therm_index(i),i);
    end
end

for i=1:beam_data.nbeams
    beam_dex=beam_data.Beam==i;
    if any(beam_data.neut_lines>0)
        P_INITIAL(i) = sum(P(1,beam_dex));
        P_PORTS(i)   = sum(P(1,and(beam_dex,port_dex)));
        P_SHINE(i)   = sum(P(1,and(beam_dex,shine_dex)));
    else
        P_PORTS(i)   = 0;
        P_SHINE(i)   = 0;
        P_INITIAL(i) = sum(P(1,and(beam_dex',beam_data.S_lines(1,:)<1)));
    end
    P_IDEPO(i)   = sum(beam_data.ipower_prof(i,:).*dVdrho)./double(beam_data.ns_prof1);
    P_EDEPO(i)   = sum(beam_data.epower_prof(i,:).*dVdrho)./double(beam_data.ns_prof1);
    P_WALL(i)    = sum(beam_data.wall_load(i,:).*A');
    P_ORBIT(i)   = sum(P_orbit(1,beam_dex));
    P_THERM(i)   = sum(P_therm(1,beam_dex));
end

P = [P_INITIAL; P_PORTS; P_SHINE; P_IDEPO; P_EDEPO; P_WALL; P_ORBIT; P_THERM];
if lverb
    disp('POWER TOTAL   PORTS   SHINE    IONS    ELEC    WALL   ORBIT   THERM');
    for i=1:size(P,2)
        disp(['BEAM' num2str(i,'%i') ':  ' num2str(round(P(:,i)./1E3)',' %6i ') ' kW']);
    end
    disp(['TOTAL:  ' num2str(round(sum(P,2)./1E3)',' %6i ') ' kW']);
end
return;
end


