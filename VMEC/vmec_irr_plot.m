function vmec_irr_plot( vmec_data, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pdata=[];
for i=1:nargin-1
    if isstruct(varargin{i})
        if isfield(varargin{i},'datatype')
            switch upper(varargin{i}.datatype)
                case 'STELLOPT_PRESSURE'
                    pdata=varargin{i};
            end
        end
    end
end


subplot(1,3,1);
if isempty(pdata)
    phin=vmec_data.phi./vmec_data.phi(vmec_data.ns);
    p_spl=pchip(phin,vmec_data.presf);
else
    dex=pdata.data(:,4)<=1.0;
    s=pdata.data(dex,4);
    p=pdata.data(dex,7);
    temp=length(s);
    s(temp+1) = 0.0;
    s(temp+2) = 1.0;
    p(temp+1) = vmec_data.presf(1);
    p(temp+2) = vmec_data.presf(vmec_data.ns);   
    eps=1.0E-4;
    for i=1:length(s)
        for j=1:length(s)
            if (s(i) == s(j) )&&(i~=j)
                s(j) = s(j) + eps;
            end
        end
    end
    [s, sdex] = sort(s);
    for i=1:length(s)
        p(i)=p(sdex(i));
    end
    datap_spl=pchip(s,p);
    phin=vmec_data.phi./vmec_data.phi(vmec_data.ns);
    vmec_data.presf=ppval(datap_spl,phin);
    p_spl=pchip(phin,vmec_data.presf);
end 
dpds_spl=spline_deriv(p_spl,1);
dpds=ppval(dpds_spl,phin);
plot(vmec_data.iotaf,-dpds,'k');
title('Pressure Gradient vs. Iota');
ylabel('Pressure Gradient')
xlabel('Rotational Transform (\iota)');
axis tight
hold on
%iota=1/2
plot([1 1]*(1./2.),[min(-dpds) max(-dpds)],'r');
plot([1 1]*(2./2.),[min(-dpds) max(-dpds)],'r');
plot([1 1]*(3./2.),[min(-dpds) max(-dpds)],'r');
plot([1 1]*(4./2.),[min(-dpds) max(-dpds)],'r');
plot([1 1]*(5./2.),[min(-dpds) max(-dpds)],'r');
plot([1 1]*(6./2.),[min(-dpds) max(-dpds)],'r');
%iota=1/3 2/3
plot([1 1]*(1./3.),[min(-dpds) max(-dpds)],'b');
plot([1 1]*(2./3.),[min(-dpds) max(-dpds)],'b');
plot([1 1]*(4./3.),[min(-dpds) max(-dpds)],'b');
plot([1 1]*(5./3.),[min(-dpds) max(-dpds)],'b');
plot([1 1]*(7./3.),[min(-dpds) max(-dpds)],'b');
plot([1 1]*(8./3.),[min(-dpds) max(-dpds)],'b');
%iota=1/4
plot([1 1]*(1./4.),[min(-dpds) max(-dpds)],'g');
plot([1 1]*(3./4.),[min(-dpds) max(-dpds)],'g');
plot([1 1]*(5./4.),[min(-dpds) max(-dpds)],'g');
plot([1 1]*(7./4.),[min(-dpds) max(-dpds)],'g');
plot([1 1]*(9./4.),[min(-dpds) max(-dpds)],'g');
plot([1 1]*(11./4.),[min(-dpds) max(-dpds)],'g');
%iota=1/5
plot([1 1]*(1./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(2./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(3./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(4./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(6./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(7./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(8./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(9./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(11./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(12./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(13./5.),[min(-dpds) max(-dpds)],'c');
plot([1 1]*(14./5.),[min(-dpds) max(-dpds)],'c');

subplot(1,3,2)
plot(phin,-dpds,'k');
title('Pressure Gradient');
ylabel('Pressure Gradient')
xlabel('Normalized Toroidal Flux');
axis tight

subplot(1,3,3)
plot(phin,vmec_data.iotaf,'k');
title('Rotational Transform (\iota)');
ylabel('Rotational Transform')
xlabel('Normalized Toroidal Flux');
axis tight

end