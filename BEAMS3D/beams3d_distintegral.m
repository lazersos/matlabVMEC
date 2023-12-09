function dist = beams3d_distintegral(beam_data,dims)
%BEAMS3D_DISTINTEGRAL Computes integrals over the BEAMS3D distribution
%   This subroutine computes integrals over the beams3d distribution
%   funciton stored as [beam,s,u,phi,vll,vperp] with units of s^3/m^6.  The
%   functions BEAMS3D_CALC_AMINOR and BEAMS3D_CALC_RMAJOR are used to get
%   the physical units.  There only the Rmajor on axis is used to perform
%   the R*dphi integral.
% Example usage
%   beam_data = read_beams3d('beams3d_filename.h5');
%   dist = beams3d_distintegral(beam_data,[1 4 5 6]) %rho vs U.
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.20

dist = beam_data.dist_prof;

% Do phi before rho because R(rho) dependence

if any(dims==1) % beam index
        dist = sum(dist,1);
end
if any(dims==4) % phi
    Rmajor = beams3d_calc_Rmajor(beam_data);
    x_arr = linspace(0,2.*pi,beam_data.ns_prof3+1).*Rmajor(1); % axis value
    dx    = diff(x_arr(1:2));
    %x_arr = 0.5.*(x_arr(1:end-1)+x_arr(2:end));
    for j = 1:beam_data.ns_prof(1)
        dist(:,j,:,:,:,:) = dist(:,j,:,:,:,:).*Rmajor(j);
    end
    dist = sum(dist,4).*dx;
end
if any(dims==2) % rho
        Aminor = beams3d_calc_aminor(beam_data);
        x_arr = linspace(0,1,beam_data.ns_prof1+1).*Aminor;
        dx    = diff(x_arr(1:2));
        x_arr = 0.5.*(x_arr(1:end-1)+x_arr(2:end));
        temp = size(dist);
        temp(2) = 1;
        temp = zeros(temp);
        for j = 1:beam_data.ns_prof1
            temp(:,1,:,:,:,:) = temp(:,1,:,:,:,:) + dist(:,j,:,:,:,:).*dx.*x_arr(j);
        end
        dist = temp;
end
if any(dims==3)% U
        x_arr = linspace(0,2.*pi,beam_data.ns_prof2+1);
        dx    = diff(x_arr(1:2));
        %x_arr = 0.5.*(x_arr(1:end-1)+x_arr(2:end));
        dist = sum(dist,3).*dx;
end
if any(dims==5)% Vll
    x_arr = linspace(-1,1,beam_data.ns_prof4+1).*beam_data.partvmax;
    dx    = diff(x_arr(1:2));
    %x_arr = 0.5.*(x_arr(1:end-1)+x_arr(2:end));
    dist = sum(dist,5).*dx;
end
if any(dims==6) % Vll
    x_arr = linspace(0,1,beam_data.ns_prof5+1).*beam_data.partvmax;
    dx    = diff(x_arr(1:2));
    x_arr = 0.5.*(x_arr(1:end-1)+x_arr(2:end));
    temp = size(dist);
    temp(6) = 1;
    temp = zeros(temp);
    for j = 1:beam_data.ns_prof5
        temp(:,:,:,:,:,1) = temp(:,:,:,:,:,1) + dist(:,:,:,:,:,j).*dx.*x_arr(j);
    end
    dist = temp;
end

dist = squeeze(dist);
return;


end