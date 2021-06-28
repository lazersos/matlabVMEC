function jout = johmic(s,te,ne,I0,t,tau)
%JOHMIC Calculates an approximate Ohmic plasma response
%   The JOHMIC routine calculates the approximate ohmic plasma response
%   given an electron temperature and density profiles.  It takes the
%   radial coordinate (s), electron temperature (te), electron density
%   (ne), total non-drive plasma current (I0), the time (t), and the
%   resistive diffusion time L/R (tau) as input.  It outputs the values on
%   the input radial grid.  Values are normalized to ne and te so any units
%   of density and temperature may be used.
%
% Example usage
%      s    = 0:0.01:1;
%      te   = 1000.*polyval([-1 1],s);
%      ne   = 1E19.*polyval([-1 2 1],s);
%      t    = 5; tau = 30; I0=15E3;
%      jout = johmic(s,te,ne,I0,t,tau);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0


jnorm = cumtrapz(s,ne.*te.^1.5);
joh = -I0.*ne.*(te.^1.5).*exp(-t./tau);
jout = joh./jnorm(end);


end

