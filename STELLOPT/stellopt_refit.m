function stellopt_refit(filename,opt_data)
%STELLOPT_REFIT Refits polynomial profiles for a new run.
%   STELLOPT_REFIT(filename,opt_data)   This function calculates the
%   polynomial fit of the pressure data to the normalized flux stored in
%   opt_data.  It then writes a stellopt input file which initializes all
%   values to those in filename except it uses the new polynomial fits.
%   Note that the new file is renamed from .min to _refit where filename is
%   the full name of the .min file.
%
%   Example:
%       pdata=read_stellopt('p_prof.test.min');
%       stellopt_refit('input.test.min',pdata);
%
%   See also read_stellopt, read_vmec_input, and read_namelist.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           6/06/11


% Defaults
numargs=2;
polyorder=10;

% Check opt_data
if isstruct(opt_data) && isfield(opt_data,'datatype')
    if ~strcmp(opt_data.datatype,'stellopt_pressure')
        disp('ERROR: opt_data is not of type stellopt_pressure');
        return
    end
else
    disp('ERROR: Check opt_data for datatype field');
    return
end

% Calculate poly fit for pressure
s=opt_data.data(opt_data.data(:,4) <= 1.0,4);
press=opt_data.data(opt_data.data(:,4) <= 1.0,5)*opt_data.norm;
poly=polyfit(s,press,polyorder); %Fit over s space
am=fliplr(poly); % Flip to get am in VMEC order
poly_data=polyval(poly,s);

% Make a plot
figure;
plot(opt_data.data(:,4),opt_data.data(:,5)*opt_data.norm,'ok');
hold on
plot(s,poly_data,'k');
xlabel('Normalized Flux (s)');
ylabel('Pressure [Pa]');

% Read STELLOPT file
vmec_input=read_vmec_input(filename);
stellopt_input=read_namelist(filename,'OPTIMUM');

% Plot old fit
poly_data_old=polyval(fliplr(vmec_input.am),s);
plot(s,poly_data_old,'r');
hold off
legend('Thompson','New Fit','Old fit');
pause(3.0);

% Modify VMEC input
vmec_input.am=am;

% Move up a directory
cd('..');
disp('    - Moving up a directory to output!');

% Write STELLOPT File
newfilename=[filename(1:length(filename)-4) '_refit'];
fid=fopen(newfilename,'w+');
% VMEC INDATA
fprintf(fid,'%s\n','&INDATA');
fprintf(fid,'%s\n','!----- Runtime Parameters -----');
write_namelist_flt(fid,'delt',vmec_input.delt);
write_namelist_int(fid,'niter',vmec_input.niter);
write_namelist_int(fid,'nstep',vmec_input.nstep);
write_namelist_flt(fid,'tcon0',vmec_input.tcon0);
write_namelist_vec(fid,'ns_array',vmec_input.ns_array,'int');
write_namelist_vec(fid,'ftol_array',vmec_input.ftol_array);
fprintf(fid,'%s\n','!----- Grid Parameters -----');
fprintf(fid,'%s\n','  LASYM = F');
write_namelist_int(fid,'nfp',vmec_input.nfp);
write_namelist_int(fid,'mpol',vmec_input.mpol);
write_namelist_int(fid,'ntor',vmec_input.ntor);
if isfield(vmec_input,'ntheta')
    write_namelist_int(fid,'ntheta',vmec_input.ntheta);
end
if isfield(vmec_input,'nzeta')
    write_namelist_int(fid,'nzeta',vmec_input.nzeta);
end
write_namelist_flt(fid,'phiedge',vmec_input.phiedge);
fprintf(fid,'%s\n','!----- Free Boundary Parameters -----');
if vmec_input.lfreeb
    write_namelist_boo(fid,'LFREEB',1);
    write_namelist_str(fid,'MGRID_FILE',vmec_input.mgrid_file);
    write_namelist_vec(fid,'EXTCUR',vmec_input.extcur);
    write_namelist_int(fid,'NVACSKIP',vmec_input.nvacskip);
else
    write_namelist_boo(fid,'LFREEB',0);
end
fprintf(fid,'%s\n','!----- Pressure Parameters -----');
write_namelist_flt(fid,'gamma',vmec_input.gamma);
write_namelist_flt(fid,'bloat',vmec_input.bloat);
write_namelist_flt(fid,'spres_ped',vmec_input.spres_ped);
write_namelist_vec(fid,'am',vmec_input.am);
fprintf(fid,'%s\n','!----- Current/Iota Parameters -----');
write_namelist_flt(fid,'curtor',vmec_input.curtor);
write_namelist_int(fid,'ncurr',vmec_input.ncurr);
write_namelist_int(fid,'AC_FORM',vmec_input.ac_form);
write_namelist_vec(fid,'AI',vmec_input.ai);
write_namelist_vec(fid,'AC',vmec_input.ac);
fprintf(fid,'%s\n','!----- Axis Parameters -----');
write_namelist_vec(fid,'RAXIS',vmec_input.raxis);
write_namelist_vec(fid,'ZAXIS',vmec_input.zaxis);
fprintf(fid,'%s\n','!----- Boundary Parameters -----');
write_namelist_arr(fid,'RBC',vmec_input.rbc,-vmec_input.ntor,vmec_input.ntor,0,vmec_input.mpol);
write_namelist_arr(fid,'ZBS',vmec_input.zbs,-vmec_input.ntor,vmec_input.ntor,0,vmec_input.mpol);
fprintf(fid,'%s\n','/');
% STELLOPT OPTIMUM
fprintf(fid,'&OPTIMUM\n');
fprintf(fid,'!-----------------------------------------------------------------------\n');
fprintf(fid,'!          OPTIMIZER RUN CONTROL PARAMETERS\n');
fprintf(fid,'!-----------------------------------------------------------------------\n');
write_namelist_flt(fid,'EPSFCN',stellopt_input.epsfcn);
write_namelist_int(fid,'NITER_OPT',stellopt_input.niter_opt);
write_namelist_int(fid,'NUM_PROCESSORS',stellopt_input.num_processors);
write_namelist_int(fid,'NUM_LEVMAR_PARAMS',stellopt_input.num_levmar_params);
write_namelist_int(fid,'NOPT_ALG',stellopt_input.nopt_alg);
write_namelist_int(fid,'NOPT_BOUNDARY',stellopt_input.nopt_boundary);
write_namelist_boo(fid,'LRESET_OPT',stellopt_input.lreset_opt);
write_namelist_boo(fid,'LDIAG_OPT',stellopt_input.ldiag_opt);
write_namelist_boo(fid,'LKEEP_MINS',stellopt_input.lkeep_mins);
fprintf(fid,'!-----------------------------------------------------------------------\n');
fprintf(fid,'!          PHYSICS MODULES\n');
fprintf(fid,'!-----------------------------------------------------------------------\n');
write_namelist_boo(fid,'LBETA_MIN',stellopt_input.lbeta_min);
write_namelist_boo(fid,'LCUR_PROF_OPT',stellopt_input.lcur_prof_opt);
write_namelist_boo(fid,'LCUR_OPT_EDGE0',stellopt_input.lcur_opt_edge0);
fprintf(fid,'  AC_MASK =  1  0  1  0  1  0  0  0  0  0  0\n');
write_namelist_boo(fid,'LDIAGNO_OPT',stellopt_input.ldiagno_opt);
write_namelist_boo(fid,'LPRES_PROF_OPT',stellopt_input.lpres_prof_opt);
write_namelist_boo(fid,'LPRES_OPT_EDGE0',stellopt_input.lpres_opt_edge0);
write_namelist_boo(fid,'LPRES_OPT_EDGEGR0',stellopt_input.lpres_opt_edgegr0);
write_namelist_boo(fid,'LP_PROF_INCL_EDGE',stellopt_input.lp_prof_incl_edge);
fprintf(fid,'!-----------------------------------------------------------------------\n');
fprintf(fid,'!          EQUILIBRIUM AND GEOMETRY OPTIMIZATION PARAMETERS\n');
fprintf(fid,'!-----------------------------------------------------------------------\n');
write_namelist_flt(fid,'TARGET_BETA',stellopt_input.target_beta);
write_namelist_flt(fid,'SIGMA_BETA',stellopt_input.sigma_beta);
write_namelist_flt(fid,'TARGET_CURTOR',stellopt_input.target_curtor);
write_namelist_flt(fid,'SIGMA_CURTOR',stellopt_input.sigma_curtor);
write_namelist_flt(fid,'TARGET_EPLASMA',stellopt_input.target_eplasma);
write_namelist_flt(fid,'SIGMA_EPLASMA',stellopt_input.sigma_eplasma);
fprintf(fid,'!-----------------------------------------------------------------------\n');
fprintf(fid,'!          PRESSURE PROFILE OPTIMIZATION PARAMETERS\n');
fprintf(fid,'!-----------------------------------------------------------------------\n');
write_namelist_int(fid,'PRES_OPT_NMAX',11);
write_namelist_int(fid,'NP_PROF',stellopt_input.np_prof);
write_namelist_flt(fid,'FACTOR_P_PROF',stellopt_input.factor_p_prof);
write_namelist_vec(fid,'P_PROF',stellopt_input.p_prof);
write_namelist_vec(fid,'SIGMA_P_PROF',stellopt_input.sigma_p_prof);
write_namelist_vec(fid,'R_P_PROF',stellopt_input.r_p_prof);
write_namelist_vec(fid,'Z_P_PROF',stellopt_input.z_p_prof);
write_namelist_vec(fid,'PHI_P_PROF',stellopt_input.phi_p_prof);
if isfield(stellopt_input,'diagno_control');
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    fprintf(fid,'!          DIAGNO OPTIMIZATION PARAMETERS\n');
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    write_namelist_str(fid,'DIAGNO_CONTROL',stellopt_input.diagno_control);
    write_namelist_vec(fid,'TARGET_DIAGNO_FLX',stellopt_input.target_diagno_flx);
    write_namelist_vec(fid,'SIGMA_DIAGNO_FLX',stellopt_input.sigma_diagno_flx);
end
fprintf(fid,'/\n'); % End Optimizer namelist

end

