function plot_profs(type,array)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

s=0:0.001:1.0;
switch lower(type)
    case 'power_series'
        subplot(2,1,1);
        plot(s,polyval(fliplr(array),s),'k');
        subplot(2,1,2);
        bar(array);
    case 'two_power'
        subplot(2,1,1);
        plot(s,array(1).*(1-s.^array(2)).^array(3),'k');
        subplot(2,1,2);
        bar(array);
    case 'pedestal'
        val = polyval(fliplr(array(1:15)),s);
        i=17;
        if (array(i+3) <= 0.0)
            i1 = 0.0;
            i2 = 0.0;
            i3 = 0.0;
            i4 = 0.0;
            i3 = 1.0E30;
        else
            i1 = array(i+1);
            i2 = array(i+2);
            i3 = array(i+3);
            i4 = 1.0 / (tanh(2*i2/i3)-tanh(2*(i2-1)/i3));
        end
        val = val + i4*i1*( tanh(2*(i2-sqrt(s))./i3) - tanh(2*(i2-1.0)./i3));
        plot(s,val,'k');
    otherwise
        disp(['!!!!!  Error unknown plot type: ' type ' !!!!!']);
end
    


end

