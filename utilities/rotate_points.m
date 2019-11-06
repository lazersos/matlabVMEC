function pts_out = rotate_points(pts,vec,angle)
%ROTATE_POINTS(pts,vec,angle) Rotates points about vector through angle.
%   The ROTATE_POINTS funciton rotates an array of 3D points of size (N,3)
%   about a vector through and angle (rad).
%
% Example usage
%      pts = rand(100,3);
%      vec = [0, 0, 1];
%      pts_out = rotate_points(pts,vec,0.1);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.00

pts_temp=pts;
ct = cos(angle); st = sin(angle);
ct1 = (1-ct);
n  = 1./norm(vec);
ux = vec(1).*n;
uy = vec(2).*n;
uz = vec(3).*n;
m = [ct + ux*ux*ct1, ux*uy*ct1-uz*st, ux*uz*ct1+uy*st;...
    uy*ux*ct1+uz*st, ct+uy*uy*ct1,    uy*uz*ct1-ux*st;...
    uz*ux*ct1-uy*st, uz*uy*ct1+ux*st, ct+uz*uz*ct1];
if size(pts,1) == 3
    pts_out=pts'*m;
elseif size(pts,2) == 3
    pts_out=pts*m;
else
    pts_out=[];
    disp('Error: size(pts)/=(N,3) or (3,N)');
end
return
    
end

