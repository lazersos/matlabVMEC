function [ output_args ] = stellopt_prof_refit( filename, s_vals )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

minf=0.8;
maxf=1.2;

switch filename(1:5)
    case{'tprof'}
        data = importdata(filename,' ',1);
        s= data.data(:,1);
        ne= data.data(:,2);
        te= data.data(:,3);
        ti= data.data(:,4);
        zeff= data.data(:,5);
        disp(['  NE_AUX_S = ' num2str(s_vals,'%20.10E')]);
        disp(['  NE_AUX_F = ' num2str(pchip(s,ne,s_vals)./1E18,'%20.10E')]);
        disp(['  TE_AUX_S = ' num2str(s_vals,'%20.10E')]);
        disp(['  TE_AUX_F = ' num2str(pchip(s,te,s_vals),'%20.10E')]);
        disp(['  TI_AUX_S = ' num2str(s_vals,'%20.10E')]);
        disp(['  TI_AUX_F = ' num2str(pchip(s,ti,s_vals),'%20.10E')]);
        disp(['  NE_F_MIN = ' num2str(minf.*pchip(s,ne,s_vals)./1E18,'%20.10E')]);
        disp(['  NE_F_MAX = ' num2str(maxf.*pchip(s,ne,s_vals)./1E18,'%20.10E')]);
        disp(['  TE_F_MIN = ' num2str(minf.*pchip(s,te,s_vals),'%20.10E')]);
        disp(['  TE_F_MAX = ' num2str(maxf.*pchip(s,te,s_vals),'%20.10E')]);
        disp(['  TI_F_MIN = ' num2str(minf.*pchip(s,ti,s_vals),'%20.10E')]);
        disp(['  TI_F_MAX = ' num2str(maxf.*pchip(s,ti,s_vals),'%20.10E')]);
    case{'dprof'}
        data = importdata(filename,' ',1);
        s= data.data(:,1);
        ex= data.data(:,2);
        disp(['  EMIS_XICS_S = ' num2str(s_vals,'%20.10E')]);
        disp(['  EMIS_XICS_F = ' num2str(pchip(s,ex,s_vals),'%20.10E')]);
        disp(['  EMIS_XICS_F_MIN = ' num2str(minf.*pchip(s,ex,s_vals),'%20.10E')]);
        disp(['  EMIS_XICS_F_MAX = ' num2str(maxf.*pchip(s,ex,s_vals),'%20.10E')]);
    case{'stell'}
        data=read_stellopt(filename);
        niter = length(data.iter);
        if isfield(data,'NE_S')
            s = data.NE_S(niter,:);
            f = data.NE_target(niter,:);
            f = f(s<=1.0);
            s = s(s<=1.0);
            [s, dex] = sort(s);
            f = f(dex);
            disp(['  NE_AUX_S = ' num2str(s_vals,'%20.10E')]);
            disp(['  NE_AUX_F = ' num2str(pchip(s,f,s_vals),'%20.10E')]);
        end
        if isfield(data,'TE_S')
            s = data.TE_S(niter,:);
            f = data.TE_target(niter,:);
            f = f(s<=1.0);
            s = s(s<=1.0);
            [s, dex] = sort(s);
            f = f(dex);
            disp(['  TE_AUX_S = ' num2str(s_vals,'%20.10E')]);
            disp(['  TE_AUX_F = ' num2str(pchip(s,f,s_vals),'%20.10E')]);
        end
end

end

