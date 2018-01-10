function data = read_bootsj( filename )
%READ_BOOTSJ(filename) The function reads the BOOTSJ output files
%   The routine reads the answers, answers_plot, and jBbs files produced by
%   the BOOTSJ code.
%
% Example usage
%      boot_data=read_bootsj('answers.test');     % Reads BOOTSJ  file
%      boot_data=read_bootsj('answers_plot.test');     % Reads BOOTSJ  file
%      boot_data=read_bootsj('jBbs.test');     % Reads BOOTSJ  file
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

% Defaults and constants
ec = 1.6021766208e-19;
me = 9.10938356e-31;
e0 = 8.854187817e-12;
fit_poly=7;

% Number of Arguments
if (nargin == 0)
    disp('read_vmec requires a filename.');
    data=1;
    return
end

% Check to see the file exists
fid=fopen(filename,'r+');
if (fid < 0)
    disp(['ERROR: Could not find file ' filename]);
    data=fid;
    return
end

if strfind(filename,'answers_plot.')
    fclose(fid);
    temp = importdata(filename)';
    data.s = temp(1,:);
    data.gnorm = temp(2,:);
    data.amain = temp(3,:);
    data.lam_int = temp(4,:);
    data.other = temp(5,:);
    data.dIds = temp(6,:);
    data.j_grad_ne = temp(7,:);
    data.j_grad_ni = temp(8,:);
    data.j_grad_Te = temp(9,:);
    data.j_grad_Ti = temp(10,:);
    data.Q = temp(11,:);
    data.ftrapped = temp(12,:);
    data.jbsnorm = temp(13,:);
    data.te = temp(14,:);
    data.ti = temp(15,:);
    data.ne = temp(16,:);
    data.ni = temp(17,:);
    data.beta = temp(18,:);
    data.jB = temp(19,:);
    data.I = cumsum(data.dIds);
    data.ne = data.ne.*1E20;
    data.ni = data.ni.*1E20;
    data.curtor = data.I(end);
    data.ac_poly = fliplr(polyfit(data.s,data.dIds,fit_poly));
    data.ac_poly = [data.ac_poly -sum(data.ac_poly)]; % assume J(1) = 0
    data.j_fit   = polyval(fliplr(data.ac_poly),0:1/(length(data.j_grad_Ti)-1):1);
elseif strfind(filename,'answers.')
    data = read_namelist(fid,'BOOTIN');
    fgetl(fid);
    head = fgetl(fid);
    i=1;
    while strfind(head,'nh ')
         temp=fscanf(fid,'%d %e %e %e %e %e %e %e %e %e %e %e',[12 inf]);
         data.array(1:11,1:data.nbuse+1,i) = temp(2:12,:);
         i=i+1;
         head=fgetl(fid);
    end
    data.xn = temp(1,:);
    data.xm = -5:5;
    temp = fscanf(fid,'%f %e%e%e%e%e%e%e%e',[9 inf]);
    data.s = temp(1,:);
    data.R = temp(2,:);
    data.bigS = temp(3,:);
    data.ftrapped = temp(4,:);
    data.H2 = temp(5,:);
    data.amain = temp(6,:);
    data.lam_int = temp(7,:);
    data.other = temp(8,:);
    data.gnorm = temp(9,:);
    fgetl(fid);
    temp = fscanf(fid,'%f %e%e%e%e%e%e%e',[8 inf]);
    data.Q = temp(2,:);
    data.thetamax = temp(3,:);
    data.zetamax = temp(4,:);
    data.Bmax = temp(5,:);
    data.ft_fttok = temp(6,:);
    data.fptok_fp = temp(7,:);
    data.cpu_secs = temp(8,:);
    fgetl(fid);
    temp = fscanf(fid,'%f %e%e%e%e%e%e%e%e',[9 inf]);
    data.gnorm = temp(2,:);
    data.jbsnorm = temp(3,:);
    data.dIds = temp(4,:);
    data.I = temp(5,:);
    data.j_grad_ne = temp(6,:);
    data.j_grad_ni = temp(7,:);
    data.j_grad_Te = temp(8,:);
    data.j_grad_Ti = temp(9,:);
    fgetl(fid);
    temp = fscanf(fid,'%d %e%e%e%e%e%e',[7 inf]);
    data.te = temp(2,:);
    data.ti = temp(3,:);
    data.ne = temp(4,:);
    data.ni = temp(5,:);
    data.beta = temp(6,:);
    data.jB = temp(7,:);
    % Add some routines
    data.ne = data.ne.*1E20;
    data.ni = data.ni.*1E20;
    data.curtor = data.I(end);
    f            = data.I;
    pfit         = polyfit(data.s,f,fit_poly);
    pfit         = polyder(pfit);
    data.ac_poly = fliplr(pfit);
    %data.ac_poly = fliplr(polyfit(data.s,data.dIds,fit_poly));
    data.ac_poly(1) = 0;
    data.ac_poly = [data.ac_poly -sum(data.ac_poly)]; % assume J(1) = 0
    data.ac_poly = data.I(end).*data.ac_poly./polyval(polyint(fliplr(data.ac_poly)),1);
    data.j_fit   = polyval(fliplr(data.ac_poly),0:1/(length(data.j_grad_Ti)-1):1);
    data.omega_pe = sqrt(data.ne.*ec.*ec./(me.*e0));
    data.lambda   = sqrt(e0.*ec.*data.te./(data.ne.*ec.*ec));
    data.lambda   = data.ne.*data.lambda.^3;
    data.nu_ei    = data.omega_pe.*log10(data.lambda)./(16.*2.*pi.*data.lambda.*pi.*2);
    data.E        = sqrt(data.s);
    data.nu_star  = data.nu_ei.*sqrt(me./(ec.*data.te))./(data.E.^3./data.Q);
    data.j_nu_star = data.dIds./(1+data.nu_star);
    fclose(fid);
elseif strfind(filename,'jBbs.')
    fgetl(fid);
    temp = fscanf(fid,'%d %e %e',[3 inf]);
    data.jbsnorm = temp(3,:);
    data.jB = temp(2,:);
    fclose(fid);
else
    disp('Unknown file type or version.');
    data=2;
    return
end

return;

end

