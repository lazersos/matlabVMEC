function w7x_beam_params( source ,varargin)
%W7X_BEAM_PARAMS(SOURCE) Generates the W7-X beamline geometry.
%   The W7X_BEAM_PARAMS(SOURCE) code generates the beamline geometry for
%   the W7-X neutral beams.  The SOURCE parameter is the source number
%   1-4 (NI20 Box) or 5-8 (NI21 Box).  A list of sources may also be
%   specified.  Note that the code does nothing if 'plots' or
%   'write_beams3d' is not specified as an optional argument.
%   Options:
%       'H2':               Hydrogen Beams
%       'D2':               Deterium Beams    
%       'plots':            Generate Geometry Plots
%       'ruidx':            Include Rudix Geometry
%       'write_beams3d':    Generate BEAMS3D Input
%       'grid':             Specify Accelerating Voltage (60 or 100)
%
%   Usage:
%       w7x_beams_params(1:8,'plots','H2','write_beams3d','grid',60);
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Version:    1.2
%   Date:       06/08/16

% Defaults
lplots = 0;
lrudix = 0;
lwrite_beams3d = 0;
grid = 60;
species ='H2';
ec = 1.602176565E-19;
div = deg2rad(4.5);
adist = 1.0;
asize = 0.5;
charge = 1.60217733E-19;

% Handle varargin
if nargin > 1
    for i=1:nargin-1
        switch varargin{i}
            case {'H2','D2'}
                species=varargin{i};
            case 'plots'
                lplots=1;
            case 'write_beams3d'
                lwrite_beams3d=2;
            case 'grid'
                i=i+1;
                grid=varargin{i};
            case {'Rudix','RUDIX','rudix'}
                lrudix=1;
        end
    end
end

% Geometry 
xo_NI20 =  3.68581; yo_NI20 =  5.65498; zo_NI20 = -0.305; %Origin locations
xo_NI21 =  0.34219; yo_NI21 =  6.74132; zo_NI21 =  0.305;
xo_RUDI = -3.81740; yo_RUDI = -7.05885; zo_RUDI = -0.70729; %RUDIX
xt_RUDI = -2.38456; yt_RUDI = -4.72680; zt_RUDI =  0.19521; %RUDIX

alpha_NI20 = 7.4441; %deg
alpha_NI21 = 7.4441; %deg

% Calculate U
nx_NI20 = -xo_NI20; ny_NI20 = -yo_NI20; nz_NI20 = 0.0;
nx_NI21 = -xo_NI21; ny_NI21 = -yo_NI21; nz_NI21 = 0.0;
ux_NI20 = xo_NI20-cosd(alpha_NI20); uy_NI20 = yo_NI20-sind(alpha_NI20);
ux_NI21 = xo_NI21+cosd(alpha_NI21); uy_NI21 = yo_NI21+sind(alpha_NI21);

% Calculate V
vx_NI20 = -uy_NI20; vy_NI20 = ux_NI20;
vx_NI21 = -uy_NI21; vy_NI21 = ux_NI21;

% Normalize
ln_NI20 = sqrt(nx_NI20.*nx_NI20+ny_NI20.*ny_NI20);
ln_NI21 = sqrt(nx_NI21.*nx_NI21+ny_NI21.*ny_NI21);
lu_NI20 = sqrt(ux_NI20.*ux_NI20+uy_NI20.*uy_NI20);
lu_NI21 = sqrt(ux_NI21.*ux_NI21+uy_NI21.*uy_NI21);
lv_NI20 = sqrt(vx_NI20.*vx_NI20+vy_NI20.*vy_NI20);
lv_NI21 = sqrt(vx_NI21.*vx_NI21+vy_NI21.*vy_NI21);
nx_NI20 = nx_NI20./ln_NI20; ny_NI20 = ny_NI20./ln_NI20;
nx_NI21 = nx_NI21./ln_NI21; ny_NI21 = ny_NI21./ln_NI21;
ux_NI20 = ux_NI20./lu_NI20; uy_NI20 = uy_NI20./lu_NI20;
ux_NI21 = ux_NI21./lu_NI21; uy_NI21 = uy_NI21./lu_NI21;
vx_NI20 = vx_NI20./lv_NI20; vy_NI20 = vy_NI20./lv_NI20;
vx_NI21 = vx_NI21./lv_NI21; vy_NI21 = vy_NI21./lv_NI21;

% order [Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8]
su = [6.5   6.5   6.5   6.5   6.5   6.5   6.5   6.5];
sv = [0.47 -0.47 -0.47  0.47 -0.47  0.47  0.47 -0.47];
sw = [0.60  0.60 -0.60 -0.60 -0.60 -0.60  0.60  0.60];

% Targets
tu = [-0.5000 -0.5000 -0.5000 -0.5000 -0.5000 -0.5000 -0.5000 -0.5000];
tv = [-0.0362  0.0362  0.0362 -0.0362  0.0362 -0.0362 -0.0362  0.0362];
tw = [0 0 0 0 0 0 0 0];


if lplots
    plot3([0 xo_NI20*1.25],[0 yo_NI20*1.25],[0 0],'k--');
    xlim([-15 15]);
    ylim([-15 15]);
    zlim([-15 15]);
    hold on;
    plot3([0 xo_NI21*1.25],[0 yo_NI21*1.25], [0 0],'k--');
    quiver3(xo_NI20,yo_NI20,0.0,nx_NI20,ny_NI20,0.0,'r','MaxHeadSize',4)
    quiver3(xo_NI20,yo_NI20,0.0,ux_NI20,uy_NI20,0.0,'m','MaxHeadSize',4)
    quiver3(xo_NI20,yo_NI20,0.0,vx_NI20,vy_NI20,0.0,'m','MaxHeadSize',4)
    quiver3(xo_NI20,yo_NI20,0.0,0,0,1.0,'m','MaxHeadSize',4)
    quiver3(xo_NI21,yo_NI21,0.0,nx_NI21,ny_NI21,0.0,'b','MaxHeadSize',4)
    quiver3(xo_NI21,yo_NI21,0.0,ux_NI21,uy_NI21,0.0,'c','MaxHeadSize',4)
    quiver3(xo_NI21,yo_NI21,0.0,vx_NI21,vy_NI21,0.0,'c','MaxHeadSize',4)
    quiver3(xo_NI21,yo_NI21,0.0,0,0,1.0,'c','MaxHeadSize',4)
    axis square;
    hold off;
end

% Energy
NAME_BEAM = {'Q1' 'Q2' 'Q3' 'Q4' 'Q5' 'Q6' 'Q7' 'Q8' };
P_H2_BEAM = [0 0 1.78 1.64 0 0 1.78 1.64].*1E6;
P_D2_BEAM = [0 0 2.48 2.28 0 0 2.48 2.28].*1E6;
if grid == 60
    E_H2 = 55E3;
    E_D2 = 60E3;
    E_H2_full  = E_H2;    P_H2_full  = 0.546;
    E_H2_half  = E_H2./2; P_H2_half  = 0.309;
    E_H2_third = E_H2./3; P_H2_third = 0.145;
    E_D2_full  = E_D2;    P_D2_full  = 0.742;
    E_D2_half  = E_D2./2; P_D2_half  = 0.208;
    E_D2_third = E_D2./3; P_D2_third = 0.05;
elseif grid == 100
    E_H2 = 72E3;
    E_D2 = 100E3;
    E_H2_full  = E_H2;    P_H2_full  = 0.380;
    E_H2_half  = E_H2./2; P_H2_half  = 0.350;
    E_H2_third = E_H2./3; P_H2_third = 0.270;
    E_D2_full  = E_D2;    P_D2_full  = 0.6207;
    E_D2_half  = E_D2./2; P_D2_half  = 0.2808;
    E_D2_third = E_D2./3; P_D2_third = 0.0985;
end

% Handle non-geometric information
E_str = {'FULL' 'HALF' 'THRID'};
switch species
    case{'H2'}
        Energy(1) = E_H2_full*ec;
        Energy(2) = E_H2_half*ec;
        Energy(3) = E_H2_third*ec;
        Z         = 1.0;
        Mass      = 1.6726231E-27; 
        P_BEAM    = P_H2_BEAM;
        P_FRAC(1) = P_H2_full;
        P_FRAC(2) = P_H2_half;
        P_FRAC(3) = P_H2_third;
    case{'D2'}
        Energy(1) = E_D2_full*ec;
        Energy(2) = E_D2_half*ec;
        Energy(3) = E_D2_third*ec;
        Z         = 1.0;
        Mass      = 3.3435837E-27; 
        P_BEAM    = P_D2_BEAM;
        P_FRAC(1) = P_D2_full;
        P_FRAC(2) = P_D2_half;
        P_FRAC(3) = P_D2_third;
    otherwise
        disp('Unknown Species');
        return;
end

if (lwrite_beams3d)
    disp('!--------Universal Beam Parameters------');
    if lrudix
        beam_str = num2str(length(source)*3+3,'%2.2i');
    else
        beam_str = num2str(length(source)*3+3,'%2.2i');
    end
    disp(['  T_END_IN(1:' beam_str ') = ' beam_str '*0.001']);
    disp(['  DIV_BEAMS(1:' beam_str ') = ' beam_str '*' num2str(div,'%20.10E')]);
    disp(['  ADIST_BEAMS(1:' beam_str ') = ' beam_str '*' num2str(adist,'%20.10E')]);
    disp(['  ASIZE_BEAMS(1:' beam_str ') = ' beam_str '*' num2str(asize,'%20.10E')]);
    disp(['  MASS_BEAMS(1:' beam_str ') = ' beam_str '*' num2str(Mass,'%20.10E')]);
    disp(['  ZATOM_BEAMS(1:' beam_str ') = ' beam_str '*' num2str(Z,'%20.10E')]);
    disp(['  CHARGE_BEAMS(1:' beam_str ') = ' beam_str '*' num2str(charge,'%20.10E')]);
end

for i=1:length(source)
    if (source(i) < 5)
        xo = xo_NI20; yo = yo_NI20; zo = zo_NI20;
        ux = ux_NI20; uy = uy_NI20;
        vx = vx_NI20; vy = vy_NI20;
    else
        xo = xo_NI21; yo = yo_NI21; zo = zo_NI21;
        ux = ux_NI21; uy = uy_NI21;
        vx = vx_NI21; vy = vy_NI21;
    end
    sx = su(source(i))*ux + sv(source(i))*vx + xo;
    sy = su(source(i))*uy + sv(source(i))*vy + yo;
    sz = sw(source(i)) + zo;
    tx = tu(source(i))*ux + tv(source(i))*vx + xo;
    ty = tu(source(i))*uy + tv(source(i))*vy + yo;
    tz = tw(source(i)) + zo;
    if lplots
        hold on;
        plot3(sx,sy,sz,'ok');
        plot3([sx tx],[sy ty],[sz tz],'k');
        plot3(tx,ty,tz,'xk');
        hold off;
        if lrudix
        end
    end
    if lwrite_beams3d
        sr = sqrt(sx*sx+sy*sy);
        sp = atan2(sy,sx);
        tr = sqrt(tx*tx+ty*ty);
        tp = atan2(ty,tx);
        for j = 1:3
            beam_str = num2str((i-1)*3+j,'%2.2i');
            disp(['!----------BEAM ' NAME_BEAM{source(i)} ' (' species ' ' E_str{j} ') ----------']);
            disp(['  E_BEAMS(' beam_str ') = ' num2str(Energy(j),'%20.10E')]);
            disp(['  P_BEAMS(' beam_str ') = ' num2str(P_BEAM(source(i))*P_FRAC(j),'%20.10E')]);
            disp(['  R_BEAMS(' beam_str,',1) = ' num2str(sr,'%20.10E')]);
            disp(['  PHI_BEAMS(' beam_str,',1) = ' num2str(sp,'%20.10E')]);
            disp(['  Z_BEAMS(' beam_str,',1) = ' num2str(sz,'%20.10E')]);
            disp(['  R_BEAMS(' beam_str,',2) = ' num2str(tr,'%20.10E')]);
            disp(['  PHI_BEAMS(' beam_str,',2) = ' num2str(tp,'%20.10E')]);
            disp(['  Z_BEAMS(' beam_str,',2) = ' num2str(tz,'%20.10E')]);
        end
        
    end
end
if lrudix
    i = length(source);
    Energy = [1 0.5 1./3].*55E3;
    P_FRAC = [55 22 22];
    P_RUDI = 10E6;
    p=atan2(yo_RUDI,xo_RUDI);
    r=sqrt(xo_RUDI.^2+yo_RUDI.^2);
    r0=[r p zo_RUDI];
    p=atan2(yt_RUDI,xt_RUDI);
    r=sqrt(xt_RUDI.^2+yt_RUDI.^2);
    r1=[r p zt_RUDI];
    if lplots
        hold on;
        plot3([xo_RUDI xt_RUDI],[yo_RUDI yt_RUDI],[zo_RUDI zt_RUDI],'k');
        hold off;
    end
    if lwrite_beams3d
        for j = 1:3
            beam_str = num2str(i*3+j,'%2.2i');
            disp(['!----------BEAM RUDIX (' species ' ' E_str{j} ') ----------']);
            disp(['  DIV_BEAMS(' beam_str ') = ' num2str(deg2rad(0.7),'%20.10E')]);
            disp(['  E_BEAMS(' beam_str ') = ' num2str(Energy(j),'%20.10E')]);
            disp(['  P_BEAMS(' beam_str ') = ' num2str(P_RUDI*P_FRAC(j),'%20.10E')]);
            disp(['  R_BEAMS(' beam_str,',1) = ' num2str(r0(1),'%20.10E')]);
            disp(['  PHI_BEAMS(' beam_str,',1) = ' num2str(r0(2),'%20.10E')]);
            disp(['  Z_BEAMS(' beam_str,',1) = ' num2str(r0(3),'%20.10E')]);
            disp(['  R_BEAMS(' beam_str,',2) = ' num2str(r1(1),'%20.10E')]);
            disp(['  PHI_BEAMS(' beam_str,',2) = ' num2str(r1(2),'%20.10E')]);
            disp(['  Z_BEAMS(' beam_str,',2) = ' num2str(r1(3),'%20.10E')]);
        end
    end
        
end

return;

end

