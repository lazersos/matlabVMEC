function [R,PHI,Z,VLL,MU] = beams3d_part2gc(data,Rpart,PHIpart,Zpart,VRpart,VPHIpart,VZpart,MASS,CHARGE)
%BEAMS3D_PART2GC Converts a particle position to a gyrocenter position
%   The BEAMS3D_PART2GC takes particle positions and transforms them to
%   gyrocenter positions.  The code takes a beams3d data structure as
%   returned by READ_BEAMS3D along with particle arrays.  The particles
%   arrays define the cylindrical positions (Rpart,PHIpart,Zpart), the 
%   cylindrical velocity vectors (VRpart,VPHIpart,VZpart), the particle 
%   mass (MASS), and the particle charge (CHARGE).  These arrays can be any
%   shape.  The code returns the cylindrical coordinates (R,PHI,Z) the
%   parallel velocity (VLL) and magnetic moment (MU) in the same shape as
%   the input arrays.
%
%   Example usage
%       beam_data = read_beams3d('beams3d_test.h5');
%       MASS  = repmat(data.mass',[data.npoinc+1 1]);
%       CHARGE  = repmat(data.charge',[data.npoinc+1 1]);
%       [R,PHI,Z,VLL,MU] = beams3d_part2gc(beam_data,beam_data.R_lines, ...
%                               beam_data.PHI_lines,beam_data.Z_lines, ...
%                               beam_data.VR_lines,beam_data.VPHI_lines,...
%                               beam_data.VZ_lines,MASS,CHARGE);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

lplot=1;

% Save the dimension and vectorize
nsave = size(Rpart);
ntotal = prod(nsave);
R = reshape(Rpart,[1 ntotal]);
P = reshape(PHIpart,[1 ntotal]);
Z = reshape(Zpart,[1 ntotal]);
VR = reshape(VRpart,[1 ntotal]);
VP = reshape(VPHIpart,[1 ntotal]);
VZ = reshape(VZpart,[1 ntotal]);
M = reshape(MASS,[1 ntotal]);
C = reshape(CHARGE,[1 ntotal]);

% Fix PHI
pgrid = mod(P,data.phiaxis(end));

% Get the field components
BR   = permute(data.B_R,[2 1 3]);
BP   = permute(data.B_PHI,[2 1 3]);
BZ   = permute(data.B_Z,[2 1 3]);
BRval= interp3(data.raxis,data.phiaxis,data.zaxis,...
    BR,R,pgrid,Z);
BPval= interp3(data.raxis,data.phiaxis,data.zaxis,...
    BP,R,pgrid,Z);
BZval= interp3(data.raxis,data.phiaxis,data.zaxis,...
    BZ,R,pgrid,Z);

% Go Cartesian
X  = R.*cos(P);
Y  = R.*sin(P);
BX = BRval.*cos(P)-BPval.*sin(P);
BY = BRval.*sin(P)+BPval.*cos(P);
BZ = BZval;
VX = VR.*cos(P)-VP.*sin(P);
VY = VR.*sin(P)+VP.*cos(P);


% Create inverse B
Binv = 1.0./sqrt(BX.*BX+BY.*BY+BZ.*BZ);
Binv(R<0) = 0;

% Create V parallel
VLL = (VX.*BX+VY.*BY+VZ.*BZ).*Binv;

% Calc Vperp
VPERP = sqrt(abs(VX.*VX+VY.*VY+VZ.*VZ-VLL.*VLL));

% Calc Gyroradius
rhox = VY.*BZ - VZ.*BY;
rhoy = VZ.*BX - VX.*BZ;
rhoz = VX.*BY - VY.*BX;
rhonorm = sign(C).*sqrt(rhox.*rhox+rhoy.*rhoy+rhoz.*rhoz);
rhox = rhox./rhonorm;
rhoy = rhoy./rhonorm;
rhoz = rhoz./rhonorm;

% Recalc gyroradius
rhox = rhox.*M.*VPERP.*Binv./C;
rhoy = rhoy.*M.*VPERP.*Binv./C;
rhoz = rhoz.*M.*VPERP.*Binv./C;

% Step to gyrocenter
Z = Z + rhoz;
X = X + rhox;
Y = Y + rhoy;
R = sqrt(X.*X+Y.*Y);
P = atan2(Y,X);
MU = 0.5.*M.*VPERP.*VPERP.*Binv;


if lplot
    quiver3(X,Y,Z,rhox,rhoy,rhoz,1,"filled",'ok');
    hold on;
    plot3(X,Y,Z,'om');
end

% Reshape
R = reshape(R,nsave);
PHI = reshape(P,nsave);
Z = reshape(Z,nsave);
VLL = reshape(VLL,nsave);
MU = reshape(MU,nsave);



end