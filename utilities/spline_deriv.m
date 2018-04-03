function pp2 = spline_deriv( pp,jderiv )
%SPLINE_DERIV Calculates the jth deriviative of a MATLAB spline
%   SPLINE_DERIV(pp,jderiv) calculates the jth deriviative of a polynomial
%   spline as produced by mkpp.  The polynomial spline is provided through
%   pp and the order of a derivative is supplied by jderiv (=0,1,or 2).
%   Inspired by: http://www.mathworks.com/support/solutions/en/data/1-15NSB/index.html?product=ML
%
%   Example;
%       x=0:1/50:1;
%       y=x.^2;
%       pp=spline(x,y);   % Spline the data
%       pp1=spline_deriv(pp,1);   % First Derivative
%       pp2=spline_deriv(pp,2);   % Second Derivative
%
%   See also:  spline, mkpp and unmkpp.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           8/01/11

if (jderiv < 0 || jderiv > 2)
    disp('ERROR:  jderiv must be 0, 1, or 2.');
    return;
end
pp2=pp;
for i=1:jderiv
    [breaks,coefs,l,k,d] = unmkpp(pp2);
    % make the polynomial that describes the derivative
    pp2 = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
end
return
end

