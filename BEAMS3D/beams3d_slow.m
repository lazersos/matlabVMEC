function data=beams3d_slow(varargin)
%BEAMS3D_SLOW(beam_data,vmec_data) Calculates flux surface slowing down
%   The BEAMS3D_SLOW routine calculates the slowing down power deposition
%   assuming that the particles slow down on flux surfaces after being
%   born.
%
%   Optional Arguments
%       'plots'     : Create plots
%       'beams'     : Downselect beamlines considered
%       'mass'      : Plasma mass [kg] (assumes protons otherwise)
%                       data=beams3d_slow(beam_data,vmec_data,'beams',4:6);
%
% Example usage
%      vmec_data = read_vmec('wout_test.nc');
%      beam_data = read_beams3d('beams3d_test.h5');
%      data=beams3d_slow(beam_data,vmec_data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.00

% Helpers
me = 9.10938356D-31;
ec = 1.60217662E-19;
lplot = 0;
beam_dex = []; % Use to downselect beams
vmec_data=[];
beam_data=[];
data=[];
plasma_mass=[];

% Handle varargin
if nargin > 0
    i=1;
    while i <= nargin
        if isstruct(varargin{i})
            if isfield(varargin{i},'datatype')
                if strcmp(varargin{i}.datatype,'wout')
                    vmec_data=varargin{i};
                elseif strcmp(varargin{i}.datatype,'BEAMS3D')
                    beam_data=varargin{i};
                end
            end
        else
            switch varargin{i}
                case 'plots'
                    lplot=1;
                case 'beams'
                    i=i+1;
                    beam_dex=varargin{i};
                case 'mass'
                    i=i+1;
                    plasma_mass=varargin{i};
            end
        end
        i=i+1;
    end
end

if isempty(beam_data)
    disp('  You must provide BEAMS3D structure for Modeling');
    return;
end

if isempty(vmec_data)
    disp('  You must provide VMEC_DATA structure for Vp');
    return;
end

% Handle lack of plasma mass
if isempty(plasma_mass)
    plasma_mass = 1.6726219E-27;
    if isfield(beam_data,'plasma_mass')
        plasma_mass = beam_data.plasma_mass;
    end
end

%Downselect beams
if isempty(beam_dex)
    beam_dex = unique(double(beam_data.Beam)');
end
dex = zeros(1,length(beam_data.Beam));
for i = beam_dex(1:end)
    dex = dex + (beam_data.Beam' == i);
end
dex = dex > 0;

% Extract information from BEAMS3D
NEUT   = beam_data.neut_lines(2,:);
dex    = and(NEUT == 0,dex);
R_BEAM = beam_data.R_lines(2,dex);
P_BEAM = mod(beam_data.PHI_lines(2,dex),max(beam_data.phiaxis));
Z_BEAM = beam_data.Z_lines(2,dex);
S_BEAM = beam_data.S_lines(2,dex);
W_BEAM = beam_data.Weight(dex)';
%BEAM   = double(beam_data.Beam(dex))';
SPEED  = beam_data.vll_lines(1,dex);
MASS   = beam_data.mass(dex)';
CHARGE = beam_data.charge(dex)';
myZ    = CHARGE./ec;
PITCH  = beam_data.vll_lines(3,dex)./SPEED;
E_BEAM = (0.5).*MASS.*SPEED.*SPEED;
NE   = permute(beam_data.NE,[2 1 3]);
TI   = permute(beam_data.TI,[2 1 3]);
TE   = permute(beam_data.TE,[2 1 3]);
NE_BEAM= interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    NE,R_BEAM,P_BEAM,Z_BEAM);
TE_BEAM= interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    TE,R_BEAM,P_BEAM,Z_BEAM);
TI_BEAM= interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    TI,R_BEAM,P_BEAM,Z_BEAM);

% PInj
if size(W_BEAM,2) == 1, W_BEAM=W_BEAM'; end % Flip W
Pinj = sum(E_BEAM.*W_BEAM);
Iinj = sum(CHARGE.*W_BEAM);

% Calculate Values
TE3=TE_BEAM.^3;
coulomb_log=[];
dex = TE_BEAM < 10.*myZ.*myZ;
coulomb_log(dex) = 23 - log(myZ(dex).*sqrt(NE_BEAM(dex).*1E-6./TE3(dex)));
dex = ~dex;
coulomb_log(dex) = 24 - log(sqrt(NE_BEAM(dex).*1E-6)./TE_BEAM(dex));
coulomb_log(coulomb_log <=1) = 1;
%v_crit = ((1.5.*sqrt(pi.*plasma_mass./me)).^(1./3.)).*sqrt(2.*TE_BEAM.*ec./MASS);
v_crit = ((0.75.*sqrt(pi.*plasma_mass./me)).^(1./3.)).*sqrt(2.*TE_BEAM.*ec./MASS);
%v_crit = ((0.75.*sqrt(pi.*MASS./me).*MASS./plasma_mass).^(1./3.)).*sqrt(2.*TE_BEAM.*ec./MASS);
%v_crit = ((0.75*sqrt(pi).*me./MASS).^(1./3)).*sqrt(TE_BEAM).*5.93096892024E5;
vcrit_cube = v_crit.^3;
tau_spit = 3.777183E41.*MASS.*sqrt(TE3)./(NE_BEAM.*myZ.*myZ.*coulomb_log);

% Integrate
C1 = 1./tau_spit;
C2 = vcrit_cube./tau_spit;
%v_sound = 1.5*sqrt(ec.*TI_BEAM./MASS);
v_sound = sqrt(1.5)*sqrt(ec.*TI_BEAM./MASS);
V  = SPEED;
V2 = V;
dt = 1E-4;
Ee = zeros(1,length(W_BEAM));
Ei = zeros(1,length(W_BEAM));
jb = zeros(1,length(W_BEAM));
t=0;
while any(V > v_sound)
    t = t+dt;
    dex = V > v_sound;
    dve = C1(dex).*V(dex);
    dvi = C2(dex)./(V(dex).*V(dex));
    dvt = dve+dvi;
    V2(dex) = V(dex) - dvt.*dt;
    Ee(dex) = Ee(dex) + V(dex).*dve.*dt;
    Ei(dex) = Ei(dex) + V(dex).*dvi.*dt;
    jb(dex) = jb(dex) + V(dex).*PITCH(dex).*dt;
    V = max(V2,v_sound);
    disp([num2str(t) ' ' num2str(V(1:3))]);
end
Pe = MASS.*W_BEAM.*Ee;
Pi = MASS.*W_BEAM.*Ei;
J = CHARGE.*W_BEAM.*jb;

% Define RHO
RHO_BEAM = sqrt(abs(S_BEAM));
[~,RHO] = hist(RHO_BEAM,100);
RHO=[0 RHO];

% Sum by beam and rho
PE_RHO = zeros(1,length(RHO));
PI_RHO = zeros(1,length(RHO));
J_RHO  = zeros(1,length(RHO));
for i = 1:length(RHO)-1
    dex = and(RHO_BEAM<RHO(i+1), RHO_BEAM>= RHO(i));
    %dex = and(dex,BEAM==k);
    PE_RHO(i+1) = sum(Pe(dex));
    PI_RHO(i+1) = sum(Pi(dex));
    J_RHO(i+1) = sum(J(dex));
end 
    

% Add total if not just one beam
%if length(unique(BEAM))>1
%    PE_RHO = [PE_RHO; sum(PE_RHO)];
%    PI_RHO = [PI_RHO; sum(PI_RHO)];
%    J_RHO = [J_RHO; sum(J_RHO)];
%end

% Calculate Vp for new grid
s = 0:1./(vmec_data.ns-1):1;
vp = pchip(s,2.*s.*vmec_data.vp*4*pi*pi./length(RHO),RHO);

if lplot
    if max(PE_RHO) > 1E6 || max(PI_RHO) > 1E6
        factor = 1.0E-6;
        units = '[MW/\Phi]';
        units2 = '[MW/m^3]';
    elseif max(PE_RHO) > 1E3 || max(PI_RHO) > 1E3
        factor = 1.0E-3;
        units = '[kW/\Phi]';
        units2 = '[kW/m^3]';
    else
        factor = 1.0;
        units = '[W/\Phi]';
        units2 = '[W/m^3]';
    end
    if max(J_RHO)>1E6
        factorj = 1E-6;
        unitsj2  = '[MA/m^2]';
    elseif max(J_RHO) > 1E3
        factorj = 1E-3;
        unitsj2  = '[kA/m^2]';
    else
        factorj = 1.0;
        unitsj2  = '[A/m^2]';
    end 
    vp2=vp;
    vp2(1)=vp2(2);
    figure('Position',[1 1 1024 768],'Color','white');
    plot(RHO,factor.*PE_RHO./vp2,'b','LineWidth',4);
    hold on;
    plot(RHO,factor.*PI_RHO./vp2,'r','LineWidth',4);
    set(gca,'FontSize',36);
    xlabel('Effective Radius (\rho/a)');
    ylabel(['Power Density ' units2]);
    title('BEAMS3D Simple Power Deposition');
    legend('P_{electrons}','P_{ions}');
    text(min(xlim)+0.025*diff(xlim),...
        max(ylim)-0.050*diff(ylim),...
        ['P_{depo} = ' num2str(Pinj./1E6,'%5.2f [MW]')],'Color','black','FontSize',36);
    text(min(xlim)+0.025*diff(xlim),...
        max(ylim)-0.110*diff(ylim),...
        ['\tau_{slow} = ' num2str(round(t.*1E3),'%4i [ms]')],'Color','black','FontSize',36);
    figure('Position',[1 1 1024 768],'Color','white');
    plot(RHO,factorj.*J_RHO./vp2,'k','LineWidth',4);
    text(min(xlim)+0.025*diff(xlim),...
        max(ylim)-0.050*diff(ylim),...
        ['I_{inj} = ' num2str(Iinj,'%5.2f [A]')],'Color','black','FontSize',36);
    set(gca,'FontSize',36);
    xlabel('Effective Radius (\rho/a)');
    ylabel(['Beam Current ' unitsj2]);
    title('BEAMS3D Simple Current');
end

data.PE = PE_RHO;
data.PI = PI_RHO;
data.RHO  = RHO;
data.Pinj = Pinj;
data.Iinj = Iinj;
data.tslow = t;
if ~isempty(vp)
    data.VP = vp;
    data.QE = PE_RHO./vp;
    data.QI = PI_RHO./vp;
    data.JB =  J_RHO./vp;
    data.QE(1) = 0;
    data.QI(1) = 0;
    data.JB(1) = 0;
end


return
end