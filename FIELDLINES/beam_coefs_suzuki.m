function output = beam_coefs_suzuki(E,Ne,Te)
%BEAM_COEFS_SUZUKI Calculates the beam stopping coefs. from Suzuki
%   The BEAM_COEFS_SUZUKI routine calcualtes the beam stopping coefficients
%   based on the model by Suzuki:
%   Suzuki, S., Shirai, T., Nemoto, M., Tobita, K., Kubo, H., Sugie, 
%       T., et al. (1999). Attenuation of high-energy neutral hydrogen 
%       beams in high-density plasmas. PPCF, 40(12), 2097?2111. 
%       http://doi.org/10.1088/0741-3335/40/12/009
%   The code take the beam energy E [keV/amu], electron density Ne
%   [m^-3], and the electron temperature Te [eV].


output=[];
eps = log(E.*1E-3);
U   = log(Te.*1E-3);
N   = Ne.*1E-19;
isotope = 'H';

coefs=[];
switch isotope
    case('H')
        coefs = [1.27E1, 1.25E0, 4.52E-1, 1.05E-2, 5.47E-1, -1.02E-1, 3.60E-1, -2.98E-2, -9.59E-2, 4.21E-3];
    case('D')
        coefs = [1.41E1, 1.11E0, 4.08E-1, 1.05E-2, 5.47E-1, -4.03E-2, 3.45E-1, -2.88E-2, -9.71E-2, 4.74E-3];
        % @1T
        %coefs = [-2.47E1, -1.30E0, -1.54E-1, 8.02E-3, 4.81E-1, -1.49E-1, 3.92E-1, -2.99E-2, -9.76E-2, 4.79E-3];
    case('T')
        coefs = [1.27E1, 1.26E0, 4.49E-1, 1.05E-2, 5.47E-1, -5.77E-3, 3.36E-1, -2.82E-2, -9.74E-2, 4.87E-3];
    otherwise
        disp(['Unknown Isotope selected: ' isotope]);
        return;
end

fact1  = 1.0+coefs(2).*eps+coefs(3).*eps.*eps;
fact2  = (1-exp(-coefs(4).*N)).^coefs(5);
fact3  = coefs(6)+coefs(7).*eps+coefs(8).*eps.*eps;
fact4  = 1.0+coefs(9).*U+coefs(10).*U.*U;
output = coefs(1).*1.0E-16.*fact1.*(1+fact2.*fact3).*fact4;

return;

end

