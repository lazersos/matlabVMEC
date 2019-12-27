function data=beams3d_slow(beam_data,varargin)
%BEAMS3D_SLOW(beam_data) Calculates flux surface slowing down
%   The BEAMS3D_SLOW routine calculates the slowing down power deposition
%   assuming that the particles slow down on flux surfaces after being
%   born.
%
%   Optional Arguments
%       'plots'     : Create plots
%       'beams'     : Downselect beamlines considered
%                       data=beams3d_bes(beam_data,file,'beams',4:6);
%
% Example usage
%      beam_data = read_beams3d('beams3d_test.h5');
%      data=beams3d_slow(beam_data);
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

% Helpers
me = 9.10938356D-31;
ec = 1.60217662E-19;
lplot = 0;
beam_dex = 1:6; % Use to downselect beams
vmec_data=[];
vp=[];

% Handle varargin
if nargin > 1
    i=1;
    while i <= nargin-1
        if isstruct(varargin{i})
            if isfield(varargin{i},'datatype')
                if strcmp(varargin{i}.datatype,'wout')
                    vmec_data=varargin{i};
                end
            end
        else
            switch varargin{i}
                case 'plots'
                    lplot=1;
                case 'beams'
                    i=i+1;
                    beam_dex=varargin{i};
                case 'geo_dex'
                    i=i+1;
                    geo_dex=varargin{i};
            end
        end
        i=i+1;
    end
end

%Downselect beams
ntotal = [];
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
W_BEAM = beam_data.Weight(dex);
%SPEED  = ones(1,length(W_BEAM)).*2.03E7;
SPEED  = beam_data.vll_lines(1,dex);
MASS   = beam_data.mass(dex)';
myZ    = beam_data.charge(dex)'./ec;
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
Pinj = sum(E_BEAM.*W_BEAM);

% Calculate Values
TE3=TE_BEAM.^3;
coulomb_log=[];
dex = TE_BEAM < 10.*myZ.*myZ;
coulomb_log(dex) = 23 - log(myZ(dex).*sqrt(NE_BEAM(dex).*1E-6./TE3(dex)));
dex = ~dex;
coulomb_log(dex) = 24 - log(sqrt(NE_BEAM(dex).*1E-6)./TE_BEAM(dex));
v_crit = ((0.75*sqrt(pi).*me./MASS).^(1./3)).*sqrt(TE_BEAM).*5.93096892024E5;
vcrit_cube = v_crit.^3;
tau_spit = 3.777183E41.*MASS.*sqrt(TE3)./(NE_BEAM.*myZ.*myZ.*coulomb_log);

% Integrate
C1 = 1./tau_spit;
C2 = vcrit_cube./tau_spit;
v_sound = 1.5*sqrt(ec.*TI_BEAM./MASS);
V  = SPEED;
V2 = V;
dt = 1E-4;
Ee = zeros(1,length(W_BEAM));
Ei = zeros(1,length(W_BEAM));
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
    V = max(V2,v_sound);
    disp([num2str(t) ' ' num2str(V(1:3))]);
end
Pe = MASS.*W_BEAM.*Ee;
Pi = MASS.*W_BEAM.*Ei;

% Sum up
RHO_BEAM = sqrt(abs(S_BEAM));
[~,RHO] = hist(RHO_BEAM,100);
RHO=[0 RHO];
for i = 1:length(RHO)-1
    dex = and(RHO_BEAM<RHO(i+1), RHO_BEAM>= RHO(i));
    PE_RHO(i) = sum(Pe(dex));
    PI_RHO(i) = sum(Pi(dex));
end
PE_RHO=[0 PE_RHO];
PI_RHO=[0 PI_RHO];

% Calculate Vp for new grid
if ~isempty(vmec_data)
    s = 0:1./(vmec_data.ns-1):1;
    vp = pchip(s,2.*s.*vmec_data.vp,RHO);
end

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
    figure('Position',[1 1 1024 768],'Color','white');
    plot(RHO,PE_RHO.*factor,'b','LineWidth',2);
    hold on;
    plot(RHO,PI_RHO.*factor,'r','LineWidth',2);
    set(gca,'FontSize',24);
    xlabel('Effective Radius (\rho/a)');
    ylabel(['Power Density' units]);
    title('BEAMS3D Simple Power Deposition');
    legend('P_{electrons}','P_{ions}');
    text(min(xlim)+0.025*diff(xlim),...
        max(ylim)-0.025*diff(ylim),...
        ['P_{depo} = ' num2str(Pinj./1E6,'%5.2f [MW]')],'Color','black','FontSize',18);
    text(min(xlim)+0.025*diff(xlim),...
        max(ylim)-0.075*diff(ylim),...
        ['\tau_{slow} = ' num2str(round(t.*1E3),'%4i [ms]')],'Color','black','FontSize',18);
    if ~isempty(vmec_data)
        vp2=vp1;
        vp2(1)=vp2(2);
        figure('Position',[1 1 1024 768],'Color','white');
        plot(RHO,factor.*PE_RHO./vp2,'b','LineWidth',2);
        hold on;
        plot(RHO,factor.*PI_RHO./vp2,'r','LineWidth',2);
        set(gca,'FontSize',24);
        xlabel('Effective Radius (\rho/a)');
        ylabel(['Power Density' units2]);
        title('BEAMS3D Simple Power Deposition');
        legend('P_{electrons}','P_{ions}');
        text(min(xlim)+0.025*diff(xlim),...
            max(ylim)-0.025*diff(ylim),...
            ['P_{depo} = ' num2str(Pinj./1E6,'%5.2f [MW]')],'Color','black','FontSize',18);
        text(min(xlim)+0.025*diff(xlim),...
            max(ylim)-0.075*diff(ylim),...
            ['\tau_{slow} = ' num2str(round(t.*1E3),'%4i [ms]')],'Color','black','FontSize',18);
    end
end

data.PE = PE_RHO;
data.PI = PI_RHO;
data.RHO  = RHO;
data.Pinj = Pinj;
data.tslow = t;
if ~isempty(vp)
    data.VP = vp;
    data.QE = PE_RHO./vp;
    data.QI = PI_RHO./vp;
    data.QE(1) = 0;
    data.QI(1) = 0;
end


return
end