function f = write_vmec_input(filename,data)
%WRITE_VMEC_INPUT(filename,data) Writes a VMEC input structure to a file.
%   This function writes a VMEC input structure to a VMEC input namelist
%   file.  It returns 1 if successful otherwise it returns a error number.
%   
%   Example:
%       data=read_vmec_input('input.test');
%       write_vmec_input('input.test_b');
%
%   See also read_vmec_input,VMECedit.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.1
%   Date:           3/28/11

%   01/12/2013 - SAL   Modified to produce user readable NS and FTOL Arrays

if ~isfield(data,'datatype')
    disp('Error:  Unknown data structure type!');
    f=-1;
    return
elseif ~strcmp(data.datatype,'VMEC_input');
    disp('Error:  Not VMEC input data !');
    disp(['    data.datatype = ' data.datatype]);
    f=-1;
    return
end

fid=fopen(filename,'w');
[message, f] = ferror(fid);
if ~(f==0)
    disp(message);
    return
end
fprintf(fid,'%s\n','&INDATA');
fprintf(fid,'%s\n','!----- Runtime Parameters -----');
write_namelist_flt(fid,'delt',data.delt);
write_namelist_int(fid,'niter',data.niter);
write_namelist_int(fid,'nstep',data.nstep);
write_namelist_flt(fid,'tcon0',data.tcon0);
ns=length(data.ns_array);
fprintf(fid,'  NS_ARRAY    =');
for i=1:ns
    fprintf(fid,' %8d ',data.ns_array(i));
end
if isfield(data,'niter_array')
    fprintf(fid,'\n  NITER_ARRAY =');
    for i=1:ns
        fprintf(fid,' %8d ',data.niter_array(i));
    end
end
fprintf(fid,'\n  FTOL_ARRAY  =');
for i=1:ns
    fprintf(fid,' %8.2E ',data.ftol_array(i));
end
fprintf(fid,'\n');
%write_namelist_vec(fid,'ns_array',data.ns_array,'int');
%write_namelist_vec(fid,'ftol_array',data.ftol_array);
fprintf(fid,'%s\n','!----- Grid Parameters -----');
write_namelist_boo(fid,'lasym',data.lasym);
%fprintf(fid,'%s\n','  LASYM = F');
write_namelist_int(fid,'nfp',data.nfp);
write_namelist_int(fid,'mpol',data.mpol);
write_namelist_int(fid,'ntor',data.ntor);
if isfield(data,'ntheta')
    write_namelist_int(fid,'ntheta',data.ntheta);
end
if isfield(data,'nzeta')
    write_namelist_int(fid,'nzeta',data.nzeta);
end
write_namelist_flt(fid,'phiedge',data.phiedge);
fprintf(fid,'%s\n','!----- Free Boundary Parameters -----');
if data.lfreeb
    write_namelist_boo(fid,'LFREEB',1);
    write_namelist_str(fid,'MGRID_FILE',data.mgrid_file);
    write_namelist_vec(fid,'EXTCUR',data.extcur);
    write_namelist_int(fid,'NVACSKIP',data.nvacskip);
else
    write_namelist_boo(fid,'LFREEB',0);
    if any(data.extcur ~= 0)
        write_namelist_vec(fid,'EXTCUR',data.extcur);
    end
    if ~isempty(data.mgrid_file)
        write_namelist_str(fid,'MGRID_FILE',data.mgrid_file);
    end
end
fprintf(fid,'%s\n','!----- Pressure Parameters -----');
write_namelist_flt(fid,'gamma',data.gamma);
write_namelist_flt(fid,'bloat',data.bloat);
write_namelist_flt(fid,'spres_ped',data.spres_ped);
write_namelist_str(fid,'pmass_type',data.pmass_type);
if isfield(data,'pres_scale')
        write_namelist_vec(fid,'pres_scale',data.pres_scale);
end
if isfield(data,'am')
        write_namelist_vec(fid,'am',data.am);
end
if isfield(data,'am_aux_s') && isfield(data,'am_aux_f')
        dex=find(data.am_aux_s == max(data.am_aux_s(4:length(data.am_aux_s))),1,'first');
        write_namelist_vec(fid,'AM_AUX_S',data.am_aux_s(1:dex));
        write_namelist_vec(fid,'AM_AUX_F',data.am_aux_f(1:dex));
end
if isfield(data,'bcrit')
    fprintf(fid,'%s\n','!----- FLOW/ANISOTROPY Parameters -----');
    write_namelist_flt(fid,'BCRIT',data.bcrit);
    write_namelist_str(fid,'ph_type',data.ph_type);
    if isfield(data,'ah')
        write_namelist_vec(fid,'ah',data.ah);
    end
    if isfield(data,'ah_aux_s') && isfield(data,'ah_aux_f')
        dex=find(data.ah_aux_s == max(data.ah_aux_s(4:length(data.ah_aux_s))),1,'first');
        write_namelist_vec(fid,'AH_AUX_S',data.ah_aux_s(1:dex));
        write_namelist_vec(fid,'AH_AUX_F',data.ah_aux_f(1:dex));
    end
    write_namelist_str(fid,'pt_type',data.pt_type);
    if isfield(data,'at')
        write_namelist_vec(fid,'at',data.ah);
    end
    if isfield(data,'at_aux_s') && isfield(data,'at_aux_f')
        dex=find(data.at_aux_s == max(data.at_aux_s(4:length(data.at_aux_s))),1,'first');
        write_namelist_vec(fid,'AT_AUX_S',data.at_aux_s(1:dex));
        write_namelist_vec(fid,'AT_AUX_F',data.at_aux_f(1:dex));
    end
end
fprintf(fid,'%s\n','!----- Current/Iota Parameters -----');
write_namelist_flt(fid,'curtor',data.curtor);
write_namelist_int(fid,'ncurr',data.ncurr);
if isfield(data,'ac_form')
    write_namelist_int(fid,'AC_FORM',data.ac_form);
end
write_namelist_str(fid,'piota_type',data.piota_type);
if isfield(data,'ai')
        write_namelist_vec(fid,'ai',data.ai);
end
if isfield(data,'ai_aux_s') && isfield(data,'ai_aux_f')
        dex=find(data.ai_aux_s == max(data.ai_aux_s(4:length(data.ai_aux_s))),1,'first');
        write_namelist_vec(fid,'AI_AUX_S',data.ai_aux_s(1:dex));
        write_namelist_vec(fid,'AI_AUX_F',data.ai_aux_f(1:dex));
end
write_namelist_str(fid,'pcurr_type',data.pcurr_type);
if isfield(data,'ac')
        write_namelist_vec(fid,'ac',data.ac);
end
if isfield(data,'ac_aux_s') && isfield(data,'ac_aux_f')
        dex=find(data.ac_aux_s == max(data.ac_aux_s(4:length(data.ac_aux_s))),1,'first');
        write_namelist_vec(fid,'AC_AUX_S',data.ac_aux_s(1:dex));
        write_namelist_vec(fid,'AC_AUX_F',data.ac_aux_f(1:dex));
end
fprintf(fid,'%s\n','!----- Axis Parameters -----');
if (isfield(data,'raxis'))
    write_namelist_vec(fid,'RAXIS',data.raxis);
elseif (isfield(data,'raxis_cc'))
    write_namelist_vec(fid,'RAXIS_cc',data.raxis_cc);
    if (data.lasym)
        write_namelist_vec(fid,'RAXIS_cs',data.raxis_cs);
    end
end
if (isfield(data,'zaxis'))
    write_namelist_vec(fid,'ZAXIS',data.zaxis);
elseif (isfield(data,'zaxis_cc'))
    write_namelist_vec(fid,'ZAXIS_cs',data.zaxis_cs);
    if (data.lasym)
        write_namelist_vec(fid,'ZAXIS_cc',data.zaxis_cc);
    end
end
fprintf(fid,'%s\n','!----- Boundary Parameters -----');
for i=1:data.mpol
    for j=1:2*data.ntor+1
        if (data.rbc(j,i)~=0.0) || (data.zbs(j,i)~=0.0) 
            fprintf(fid,'RBC(%3i,%3i) = %17.10e  ZBS(%3i,%3i) = %17.10e\n',j-1-data.ntor,i-1,data.rbc(j,i),j-1-data.ntor,i-1,data.zbs(j,i));
            if (data.lasym)
                fprintf(fid,'   RBS(%3i,%3i) = %17.10e  ZBC(%3i,%3i) = %17.10e\n',j-1-data.ntor,i-1,data.rbs(j,i),j-1-data.ntor,i-1,data.zbc(j,i));
            end
        end
    end
end
fprintf(fid,'%s\n',['!----- Created by write_vmec ' datestr(now) ' -----']);
fprintf(fid,'%s\n','/');
fclose(fid);
f=1;
return
end

