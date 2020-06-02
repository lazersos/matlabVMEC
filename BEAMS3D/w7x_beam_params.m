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
%       'Energy':           Specify Accelerating Voltage (explicit) [ev]
%       'Power':            Specify Power (explicit) [W]
%       'Powerfrac':        Specify Power fractions (explicit) [0-1]
%
%           Note: Options marked explicit should have the same shape as the
%           source vector.  Powerfrac should have shape [3:nbeams] where
%           nbeams is the legnth of sources.
%
%   Usage:
%       w7x_beams_params(1:8,'plots','H2','write_beams3d','grid',60);
%
%   Created by: S. Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:    1.5
%   Date:       06/08/16

% Defaults
lplots = 0;
lrudix = 0;
lwrite_beams3d = 0;
grid = 60;
species ='H2';
ec = 1.602176565E-19;
div = deg2rad(1);
div = 0.0125;
adist = 1.0;
asize = 0.5;
charge = 1.60217733E-19;
new_pfrac=[];
new_E=[];
new_power=[];

% Handle varargin
if nargin > 1
    i = 1;
    while i < nargin
        switch varargin{i}
            case {'H2','D2'}
                species=varargin{i};
            case 'plots'
                lplots=1;
            case 'write_beams3d'
                lwrite_beams3d=1;
            case 'grid'
                i=i+1;
                grid=varargin{i};
            case {'Rudix','RUDIX','rudix'}
                lrudix=1;
            case {'Energy','ENERGY','energy'}
                i=i+1;
                new_E = varargin{i};
            case {'Powerfrac','POWERFRAC','powerfrac'}
                i=i+1;
                new_pfrac = varargin{i};
            case {'Power','POWER','power'}
                i=i+1;
                new_power = varargin{i};
        end
        i = i+1;
    end
end

% Geometry 
xo_NI20 =  3.68581; yo_NI20 =  5.65498; zo_NI20 = -0.305; %Origin locations
xo_NI21 =  0.34219; yo_NI21 =  6.74132; zo_NI21 =  0.305;
xo_RUDI = -3.81740; yo_RUDI = -7.05885; zo_RUDI = -0.70729; %RUDIX
xt_RUDI = -2.38456; yt_RUDI = -4.72680; zt_RUDI =  0.19521; %RUDIX

alpha_NI20 = 7.4441; %deg
alpha_NI21 = 7.4441; %deg

% Now define unit vector for u
ro_NI20 = sqrt(xo_NI20*xo_NI20+yo_NI20*yo_NI20);
ro_NI21 = sqrt(xo_NI21*xo_NI21+yo_NI21*yo_NI21);
po_NI20 = atan2d(yo_NI20,xo_NI20);
po_NI21 = atan2d(yo_NI21,xo_NI21);
nx_NI20 = ro_NI20*cosd(po_NI20+alpha_NI20);
nx_NI21 = ro_NI21*cosd(po_NI21-alpha_NI21);
ny_NI20 = ro_NI20*sind(po_NI20+alpha_NI20);
ny_NI21 = ro_NI21*sind(po_NI21-alpha_NI21);
nz_NI20 = 0.0;
nz_NI21 = 0.0;
l_NI20=sqrt(nx_NI20*nx_NI20+ny_NI20*ny_NI20+nz_NI20*nz_NI20);
l_NI21=sqrt(nx_NI21*nx_NI21+ny_NI21*ny_NI21+nz_NI21*nz_NI20);
nx_NI20 = nx_NI20./l_NI20;
nx_NI21 = nx_NI21./l_NI21;
ny_NI20 = ny_NI20./l_NI20;
ny_NI21 = ny_NI21./l_NI21;
nz_NI20 = nz_NI20./l_NI20;
nz_NI21 = nz_NI21./l_NI21;

% Define v and w

ux_NI20 = nx_NI20; uy_NI20 = ny_NI20;
ux_NI21 = nx_NI21; uy_NI21 = ny_NI21;
vx_NI20 = -uy_NI20; vy_NI20 = ux_NI20; vz_NI20 = 0;
vx_NI21 = -uy_NI21; vy_NI21 = ux_NI21; vz_NI21 = 0;

% Now locate the injectors
u_q1 = 6.5; u_q2 = 6.5; u_q3=6.5; u_q4=6.5;
u_q7 = 6.5; u_q8 = 6.5; u_q5=6.5; u_q6=6.5;
v_q1 = 0.47; v_q2 = -0.47; v_q3 = -0.47; v_q4 = 0.47;
v_q7 = 0.47; v_q8 = -0.47; v_q5 = -0.47; v_q6 = 0.47;
w_q1 = 0.6; w_q2 = 0.6; w_q3 = -0.6; w_q4 = -0.6;
w_q7 = 0.6; w_q8 = 0.6; w_q5 = -0.6; w_q6 = -0.6;

% Final injector position
x_q1 = u_q1*ux_NI20+v_q1*vx_NI20+xo_NI20; y_q1 = u_q1*uy_NI20+v_q1*vy_NI20+yo_NI20; z_q1 = w_q1+zo_NI20;
x_q2 = u_q2*ux_NI20+v_q2*vx_NI20+xo_NI20; y_q2 = u_q2*uy_NI20+v_q2*vy_NI20+yo_NI20; z_q2 = w_q2+zo_NI20;
x_q3 = u_q3*ux_NI20+v_q3*vx_NI20+xo_NI20; y_q3 = u_q3*uy_NI20+v_q3*vy_NI20+yo_NI20; z_q3 = w_q3+zo_NI20;
x_q4 = u_q4*ux_NI20+v_q4*vx_NI20+xo_NI20; y_q4 = u_q4*uy_NI20+v_q4*vy_NI20+yo_NI20; z_q4 = w_q4+zo_NI20;
x_q5 = u_q5*ux_NI21+v_q5*vx_NI21+xo_NI21; y_q5 = u_q5*uy_NI21+v_q5*vy_NI21+yo_NI21; z_q5 = w_q5+zo_NI21;
x_q6 = u_q6*ux_NI21+v_q6*vx_NI21+xo_NI21; y_q6 = u_q6*uy_NI21+v_q6*vy_NI21+yo_NI21; z_q6 = w_q6+zo_NI21;
x_q7 = u_q7*ux_NI21+v_q7*vx_NI21+xo_NI21; y_q7 = u_q7*uy_NI21+v_q7*vy_NI21+yo_NI21; z_q7 = w_q7+zo_NI21;
x_q8 = u_q8*ux_NI21+v_q8*vx_NI21+xo_NI21; y_q8 = u_q8*uy_NI21+v_q8*vy_NI21+yo_NI21; z_q8 = w_q8+zo_NI21;

xstart=[x_q1, x_q2, x_q3, x_q4, x_q5, x_q6, x_q7, x_q8];
ystart=[y_q1, y_q2, y_q3, y_q4, y_q5, y_q6, y_q7, y_q8];
zstart=[z_q1, z_q2, z_q3, z_q4, z_q5, z_q6, z_q7, z_q8];
rstart=sqrt(xstart.^2+ystart.^2);
pstart=atan2(ystart,xstart);


% Target points
u_s14 = -0.5; v_s14 = -0.0362; w_s14=0;
u_s23 = -0.5; v_s23 =  0.0362; w_s23=0;
u_s67 = -0.5; v_s67 = -0.0362; w_s67=0;
u_s58 = -0.5; v_s58 =  0.0362; w_s58=0;
x_s14 = u_s14*ux_NI20+v_s14*vx_NI20+xo_NI20; y_s14 = u_s14*uy_NI20+v_s14*vy_NI20+yo_NI20; z_s14=zo_NI20;
x_s23 = u_s23*ux_NI20+v_s23*vx_NI20+xo_NI20; y_s23 = u_s23*uy_NI20+v_s23*vy_NI20+yo_NI20; z_s23=zo_NI20;
x_s58 = u_s58*ux_NI21+v_s58*vx_NI21+xo_NI21; y_s58 = u_s58*uy_NI21+v_s58*vy_NI21+yo_NI21; z_s58=zo_NI21;
x_s67 = u_s67*ux_NI21+v_s67*vx_NI21+xo_NI21; y_s67 = u_s67*uy_NI21+v_s67*vy_NI21+yo_NI21; z_s67=zo_NI21;

xtarget=[x_s14, x_s23, x_s23, x_s14, x_s58, x_s67, x_s67, x_s58];
ytarget=[y_s14, y_s23, y_s23, y_s14, y_s58, y_s67, y_s67, y_s58];
ztarget=[z_s14, z_s23, z_s23, z_s14, z_s58, z_s67, z_s67, z_s58];
rtarget=sqrt(xtarget.^2+ytarget.^2);
ptarget=atan2(ytarget,xtarget);

if lplots
    % plot box geometry
    hold on;
    if (any(source < 5))
        plot3([0 xo_NI20],[0 yo_NI20],[0 zo_NI20],'-k')
        plot3(xo_NI20,yo_NI20,zo_NI20,'ok')
        quiver3(xo_NI20,yo_NI20,zo_NI20,ux_NI20,uy_NI20,0,6.5,'r')
        quiver3(xo_NI20,yo_NI20,zo_NI20,vx_NI20,vy_NI20,0,1,'r')
        quiver3(xo_NI20,yo_NI20,zo_NI20,0,0,1,1,'r')
        plot3(x_s14,y_s14,z_s14,'+k')
        plot3(x_s23,y_s23,z_s23,'+k')
    end
    if (any(source > 4))
        plot3([0 xo_NI21],[0 yo_NI21],[0 zo_NI21],'-k')
        plot3(xo_NI21,yo_NI21,zo_NI21,'ok')
        quiver3(xo_NI21,yo_NI21,zo_NI21,ux_NI21,uy_NI21,0,6.5,'r')
        quiver3(xo_NI21,yo_NI21,zo_NI21,vx_NI21,vy_NI21,0,1,'r')
        quiver3(xo_NI21,yo_NI21,zo_NI21,0,0,1,1,'r')
        plot3(x_s58,y_s58,z_s58,'+k')
        plot3(x_s67,y_s67,z_s67,'+k')
    end
    d = 0.15;
    if (find(source == 1))
        plot3(x_q1,y_q1,z_q1,'ok'); text(x_q1,y_q1,z_q1+d,'Q1');
        quiver3(x_q1,y_q1,z_q1,[x_s14-x_q1],[y_s14-y_q1],[z_s14-z_q1],'b')
    end
    if (find(source == 2))
        plot3(x_q2,y_q2,z_q2,'ok'); text(x_q2,y_q2,z_q2+d,'Q2');
        quiver3(x_q2,y_q2,z_q2,[x_s23-x_q2],[y_s23-y_q2],[z_s23-z_q2],'b')
    end
    if (find(source == 3))
        plot3(x_q3,y_q3,z_q3,'ok'); text(x_q3,y_q3,z_q3+d,'Q3');
        quiver3(x_q3,y_q3,z_q3,[x_s23-x_q3],[y_s23-y_q3],[z_s23-z_q3],'b')
    end
    if (find(source == 4))
        plot3(x_q4,y_q4,z_q4,'ok'); text(x_q4,y_q4,z_q4+d,'Q4');
        quiver3(x_q4,y_q4,z_q4,[x_s14-x_q4],[y_s14-y_q4],[z_s14-z_q4],'b')
    end
    if (find(source == 5))
        plot3(x_q5,y_q5,z_q5,'ok'); text(x_q5,y_q5,z_q5+d,'Q5');
        quiver3(x_q5,y_q5,z_q5,[x_s58-x_q5],[y_s58-y_q5],[z_s58-z_q5],'b')
    end
    if (find(source == 6))
        plot3(x_q6,y_q6,z_q6,'ok'); text(x_q6,y_q6,z_q6+d,'Q6');
        quiver3(x_q6,y_q6,z_q6,[x_s67-x_q6],[y_s67-y_q6],[z_s67-z_q6],'b')
    end
    if (find(source == 7))
        plot3(x_q7,y_q7,z_q7,'ok'); text(x_q7,y_q7,z_q7+d,'Q7');
        quiver3(x_q7,y_q7,z_q7,[x_s67-x_q7],[y_s67-y_q7],[z_s67-z_q7],'b')
    end
    if (find(source == 8))
        plot3(x_q8,y_q8,z_q8,'ok'); text(x_q8,y_q8,z_q8+d,'Q8');
        quiver3(x_q8,y_q8,z_q8,[x_s58-x_q8],[y_s58-y_q8],[z_s58-z_q8],'b')
    end
    axis tight; axis equal;
end

% Energy
NAME_BEAM = {'Q1' 'Q2' 'Q3' 'Q4' 'Q5' 'Q6' 'Q7' 'Q8' };
P_full=[]; P_half=[]; P_third=[];
E_full=[]; E_half=[]; E_third=[];
P_H2_BEAM = [1.78 1.64 1.78 1.64 1.78 1.64 1.78 1.64].*1E6;
P_D2_BEAM = [2.48 2.28 2.48 2.28 2.48 2.28 2.48 2.28].*1E6;
if grid ==60
    switch species
        case{'H2'}
            E_full = [1 1 1 1 1 1 1 1].*55E3;
            P_BEAM = [1.78 1.64 1.78 1.64 1.78 1.64 1.78 1.64].*1E6;
            P_frac_full  = [1 1 1 1 1 1 1 1].*0.546;
            P_frac_half  = [1 1 1 1 1 1 1 1].*0.309;
            P_frac_third = [1 1 1 1 1 1 1 1].*0.145;
            Mass      = 1.6726231E-27;
        case{'D2'}
            E_full = [1 1 1 1 1 1 1 1].* 60E3;
            P_BEAM = [2.48 2.28 2.48 2.28 2.48 2.28 2.48 2.28].*1E6;
            P_frac_full  = [1 1 1 1 1 1 1 1].*0.742;
            P_frac_half  = [1 1 1 1 1 1 1 1].*0.208;
            P_frac_third = [1 1 1 1 1 1 1 1].*0.050;
            Mass      = 3.3435837E-27;
    end
elseif grid == 100
    switch species
        case{'H2'}
            E_full = [1 1 1 1 1 1 1 1].*72E3;
            P_BEAM = [1.78 1.64 1.78 1.64 1.78 1.64 1.78 1.64].*1E6; %????
            P_frac_full  = [1 1 1 1 1 1 1 1].*0.380;
            P_frac_half  = [1 1 1 1 1 1 1 1].*0.350;
            P_frac_third = [1 1 1 1 1 1 1 1].*0.270;
            Mass      = 1.6726231E-27;
        case{'D2'}
            E_full = [1 1 1 1 1 1 1 1].*100E3;
            P_BEAM = [2.48 2.28 2.48 2.28 2.48 2.28 2.48 2.28].*1E6; %????
            P_frac_full  = [1 1 1 1 1 1 1 1].*0.6207;
            P_frac_half  = [1 1 1 1 1 1 1 1].*0.2808;
            P_frac_third = [1 1 1 1 1 1 1 1].*0.0985;
            Mass      = 3.3435837E-27;
    end
end

% Handel new energy
if ~isempty(new_E)
    E_full(source) = new_E;
end
E_half  = E_full./2;
E_third = E_full./3;


% Handle non-geometric information
E_str = {'FULL' 'HALF' 'THRID'};
Energy(1,:) = E_full*ec;
Energy(2,:) = E_half*ec;
Energy(3,:) = E_third*ec;
Z = 1.0;
P_FRAC(1,:) = P_frac_full;
P_FRAC(2,:) = P_frac_half;
P_FRAC(3,:) = P_frac_third;

% Handel new Power
if ~isempty(new_power)
    P_BEAM(source) = new_power;
end

% Handel new PowerFracp
if ~isempty(new_pfrac)
    P_FRAC(1:3,source) = new_pfrac;
end

if (lwrite_beams3d)
    disp('!--------Universal Beam Parameters------');
    if lrudix
        beam_str = num2str(length(source)*3+3,'%2.2i');
    else
        beam_str = num2str(length(source)*3,'%2.2i');
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
    dex_beam = source(i);
    if lwrite_beams3d
        for j = 1:3
            beam_str = num2str((i-1)*3+j,'%2.2i');
            disp(['!----------BEAM ' NAME_BEAM{source(i)} ' (' species ' ' E_str{j} ') ----------']);
            disp(['  DEX_BEAMS(' beam_str ') = ' num2str(source(i),'%2i')]);
            disp(['  E_BEAMS(' beam_str ') = ' num2str(Energy(j,dex_beam),'%20.10E')]);
            disp(['  P_BEAMS(' beam_str ') = ' num2str(P_BEAM(dex_beam)*P_FRAC(j,dex_beam),'%20.10E')]);
            disp(['  R_BEAMS(' beam_str,',1) = ' num2str(rstart(dex_beam),'%20.10E')]);
            disp(['  PHI_BEAMS(' beam_str,',1) = ' num2str(pstart(dex_beam),'%20.10E')]);
            disp(['  Z_BEAMS(' beam_str,',1) = ' num2str(zstart(dex_beam),'%20.10E')]);
            disp(['  R_BEAMS(' beam_str,',2) = ' num2str(rtarget(dex_beam),'%20.10E')]);
            disp(['  PHI_BEAMS(' beam_str,',2) = ' num2str(ptarget(dex_beam),'%20.10E')]);
            disp(['  Z_BEAMS(' beam_str,',2) = ' num2str(ztarget(dex_beam),'%20.10E')]);
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
            %disp(['  DEX_BEAMS(' beam_str ') = ' num2str(deg2rad(0.7),'%20.10E')]);
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

