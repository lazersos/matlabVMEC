function val = vmec_profile(type,aa_in,s)
%val = VMEC_PROFILE(type,aa,s) Evaluate VMEC profiles at s.
%   Detailed explanation goes here

aa(1,1:21) = 0.0;
aa(1:size(aa_in,2)) = aa_in;
switch lower(type)
    case{'power_series'}
        val = polyval(fliplr(aa),s);
    case{'two_power'}
        val = aa(1).*(1.0-s.^aa(2)).^aa(3);
    case{'gauss_trunc'}
        glx = [0.01304673574141414, 0.06746831665550774, 0.1602952158504878,...
            0.2833023029353764, 0.4255628305091844, 0.5744371694908156,...
            0.7166976970646236, 0.8397047841495122, 0.9325316833444923,...
            0.9869532642585859];
        glw = [0.03333567215434407, 0.0747256745752903, 0.1095431812579910,...
            0.1346333596549982, 0.1477621123573764, 0.1477621123573764,...
            0.1346333596549982, 0.1095431812579910, 0.0747256745752903,...
            0.03333567215434407 ];
        for i=1:10
            xp = s.*glx(i);
            val = val + glw(i).*aa(1).*(exp(-(xp/aa(2)).^2) - exp(-(1./aa(2)).^2));
        end
        val = val.*s;
    case{'gauss_trunc_offset'}
        val = aa(1)+(aa(2)./(1. - exp(-(1./aa(4)).^2)))*...
                  (exp(-((s-aa(3))./aa(2)).^2)-exp(-(1./aa(4)).^2));
    case{'sum_atan'}
        val = aa(1)+(2./pi).*( aa(2).*atan(aa(3).*s.^aa(4)./(1-s).^aa(5))...
                              +aa(6).*atan(aa(7).*s.^aa(8)./(1-s).^aa(9))...
                              +aa(10).*atan(aa(11).*s.^aa(12)./(1-s).^aa(13))...
                              +aa(14).*atan(aa(15).*s.^aa(16)./(1-s).^aa(17))...
                              +aa(18).*atan(aa(19).*s.^aa(20)./(1-s).^aa(21)));
        val(s>=1) = aa(1) + aa(2) + aa(6) + aa(10) + aa(14) + aa(17);
    case{'nice_quadratic'}
        val  = aa(1).*(1.0-s)+aa(2).*s+4.*aa(3).*s.*(1.0-s);
    case{'two_lorentz'}
        val  = aa(1).*(aa(2).*(1./(1+(s./aa(3).^2).^aa(4)).^aa(5)...
                              -1./(1+(1./aa(3).^2).^aa(4)).^aa(5))./...
                      (1-1./(1+(1./aa(3).^2).^aa(4)).^aa(5))+...
                      (1-aa(2)).*(1./(1+(s./aa(6).^2).^aa(7)).^aa(8)...
                                 -1./(1+(1./aa(6).^2).^aa(7)).^aa(8))./...
                      (1-1./(1+(1./aa(6).^2).^aa(7)).^aa(8)));
    case{'pedestal'}
        val = polyval(fliplr(aa(1:16),s));
        i   = 17;
        if (aa(i+3) <= 0.0)
            aa(i:i+4) = 0.0;
            aa(i+3)   = 1.0E30;
        else
            aa(i+4) = 1.0 / (tanh(2*aa(i+2)/aa(i+3))-tanh(2*(aa(i+2)-1)/aa(i+3)));
        end
        val = val + aa(i+4).*aa(i+1).*( tanh(2.*(aa(i+2)-sqrt(s))./aa(i+3))...
                                       -tanh(2.*(aa(i+2)-1.0)    ./aa(i+3)));
    case{'bump'}
        x0 = aa(1);
        x1 = 1.0-2*(1.0-x0);
        x2 = 1.0;
        x3 = aa(3);
        h  = aa(2)/((x0-x1)*(x0-x2));
        val= 0.0*s;
        val(s>x1) = h.*(s(s>x1)-x1).*(s(s>x1)-x2);
        val = val + x3.*s.*(s-1)./(-0.25);
        
                                 
end
val(abs(s) > 1.0) = 0.0;
return;
end

