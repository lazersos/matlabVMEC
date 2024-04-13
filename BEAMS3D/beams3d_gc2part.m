function [R,PHI,Z,VR,VPHI,VZ] = beams3d_gc2part(data,Rpart,PHIpart,Zpart,VLLpart,MUpart,MASS,CHARGE)
%BEAMS3D_GC2PART Converts a gyrocenter position to a particle position
%   The BEAMS3D_PART2GC takes particle positions and transforms them to
%   gyrocenter positions.  The code takes a beams3d data structure as
%   returned by READ_BEAMS3D along with particle arrays.  The particles
%   arrays define the cylindrical positions (Rpart,PHIpart,Zpart), the
%   parallel velocity (VLLpart), the magnetic moment (MUpart), the particle 
%   mass (MASS), and the particle charge (CHARGE).  These arrays can be any
%   shape.  The code returns the cylindrical coordinates (R,PHI,Z) and the
%   cylindcrical velocity (VR,VPHI,VZ) in the same shape as the input
%   arrays.
%
%   Example usage
%       beam_data = read_beams3d('beams3d_test.h5');
%       MASS  = repmat(data.mass',[data.npoinc+1 1]);
%       CHARGE  = repmat(data.charge',[data.npoinc+1 1]);
%       [R,PHI,Z,VR,VPHI,VZ] = beams3d_gc2part(beam_data,beam_data.R_lines, ...
%                               beam_data.PHI_lines,beam_data.Z_lines, ...
%                               beam_data.vll_lines,beam_data.MOMENT_lines,...
%                               MASS,CHARGE);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% For plotting
lplot = 0;

% Save the dimension and vectorize
nsave = size(Rpart);
ntotal = prod(nsave);
R = reshape(Rpart,[1 ntotal]);
P = reshape(PHIpart,[1 ntotal]);
Z = reshape(Zpart,[1 ntotal]);
VLL = reshape(VLLpart,[1 ntotal]);
MU = reshape(MUpart,[1 ntotal]);
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

% Create inverse B
B = sqrt(BX.*BX+BY.*BY+BZ.*BZ);
Binv = 1.0./B;
Binv(R<0) = 0;
BX = BX.*Binv;
BY = BY.*Binv;
BZ = BZ.*Binv;

% Calculate Gyroradius
VPERP = sqrt(2.*B.*MU./M);
RGYRO = M.*VPERP.*Binv./C;

% Calculate the perpendicular vector from (BxZ)xB
XGYRO = -BX.*BZ.*RGYRO;
YGYRO = -BZ.*BY.*RGYRO;
ZGYRO = (BY.*BY+BX.*BX).*RGYRO;

if lplot
    quiver3(X(1),Y(1),Z(1),BX(1),BY(1),BZ(1),0.01,'filled','ro');
    hold on;
    quiver3(X(1),Y(1),Z(1),XGYRO(1),YGYRO(1),ZGYRO(1),1,'filled','go');
end

% Generate a random gyrophase
PGYRO = (rand([1,ntotal])-0.5).*2.*pi;

% Now Compute a rotation matrix
% https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
XG = zeros([1,ntotal]);
YG = zeros([1,ntotal]);
ZG = zeros([1,ntotal]);
ROT_MATRIX = zeros(3,3,ntotal);
COSPGRYO = cos(PGYRO);
SINPGRYO = sin(PGYRO);
OMCOSPTYRO = 1-COSPGRYO;
ROT_MATRIX(1,1,:) = COSPGRYO+BX.*BX.*OMCOSPTYRO;
ROT_MATRIX(1,2,:) = BX.*BY.*OMCOSPTYRO-BZ.*SINPGRYO;
ROT_MATRIX(1,3,:) = BX.*BZ.*OMCOSPTYRO+BY.*SINPGRYO;
ROT_MATRIX(2,1,:) = BY.*BX.*OMCOSPTYRO+BZ.*SINPGRYO;
ROT_MATRIX(2,2,:) = COSPGRYO+BY.*BY.*OMCOSPTYRO;
ROT_MATRIX(2,3,:) = BY.*BZ.*OMCOSPTYRO-BX.*SINPGRYO;
ROT_MATRIX(3,1,:) = BZ.*BX.*OMCOSPTYRO-BY.*SINPGRYO;
ROT_MATRIX(3,2,:) = BZ.*BY.*OMCOSPTYRO+BX.*SINPGRYO;
ROT_MATRIX(3,3,:) = COSPGRYO+BZ.*BZ.*OMCOSPTYRO;
for i = 1:ntotal
    A = ROT_MATRIX(:,:,i);
    XTEMP = A*[XGYRO(i);YGYRO(i);ZGYRO(i)];
    XG(i) = XTEMP(1);
    YG(i) = XTEMP(2);
    ZG(i) = XTEMP(3);
end

% Calc the perp vector
XPERP = BY.*ZG - BZ.*YG;
YPERP = BZ.*XG - BX.*ZG;
ZPERP = BX.*YG - BY.*XG;
RHOPERP = sqrt(XPERP.*XPERP+YPERP.*YPERP+ZPERP.*ZPERP);
RHOPERPinv = 1./RHOPERP;


% Step to particle position
X = X + XG;
Y = Y + YG;
Z = Z + ZG;

%Convert to Cylindrical
R = sqrt(X.*X+Y.*Y);
PHI = atan2(Y,X);
VX = VLL.*BX + VPERP.*XPERP.*RHOPERPinv;
VY = VLL.*BY + VPERP.*YPERP.*RHOPERPinv;
VZ = VLL.*BZ + VPERP.*ZPERP.*RHOPERPinv;
VR   =  VX.*cos(PHI)+VY.*sin(PHI);
VPHI = -VX.*sin(PHI)+VY.*cos(PHI);

if lplot
    quiver3(X,Y,Z,VX-VLL.*BX,VY-VLL.*BY,VZ-VLL.*BZ,1,'filled','ob');
    quiver3(X(1),Y(1),Z(1),XPERP(1),YPERP(1),ZPERP(1),1,'filled','go');
    quiver3(X,Y,Z,VLL.*BX,VLL.*BY,VLL.*BZ,1','filled','oc');
    axis equal;
end

% Reshape
R = reshape(R,nsave);
PHI = reshape(PHI,nsave);
Z = reshape(Z,nsave);
VR = reshape(VR,nsave);
VPHI = reshape(VPHI,nsave);
VZ = reshape(VZ,nsave);



end