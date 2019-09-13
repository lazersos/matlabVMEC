function f = read_ipec(filename)
%READ_IPEC(filename)
%   Read IPEC output files.

fid = fopen(filename,'r+');
if (fid < 0)
    disp(['ERROR: Could not find file ' filename]);
    f=fid;
    return
else
    fclose(fid);
end
f=-1;
% Determine routine to call
sdex=strfind(filename,'.');
if isempty(sdex)
    disp(['ERROR: Filename ' filename '!']);
    return
end
if ~isempty(strfind(filename,'ipecopt.'))
        f=read_stellopt(filename);
elseif ~isempty(strfind(filename,'idcon_equil.'))
        f=read_dcon_equil(filename);
elseif ~isempty(strfind(filename,'ipec_pmodb_fun'))
        f=read_ipec_pmodb_fun(filename);
elseif ~isempty(strfind(filename,'ipec_pmodb'))
        f=read_ipec_pmodb(filename);
elseif ~isempty(strfind(filename,'ipec_xbnormal_fun'))
        f=read_ipec_xbnormal_fun(filename);
elseif ~isempty(strfind(filename,'ipec_xbnormal'))
        f=read_ipec_xbnormal(filename);
elseif ~isempty(strfind(filename,'pent'))
        f=read_ipec_pent(filename);
else
        disp(['ERROR: Unsuported filetype' filetype '!']);
end
return
end

function f=read_dcon_equil(filename)
fid=fopen(filename,'r+');
line=fgetl(fid); % header
line=fgetl(fid);
line_type=fgetl(fid);
line_size=fgetl(fid);
line_psi=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
dex=strfind(line_type,'=');
f.jac_type = strtrim(line_type(dex+1:end));
temp=sscanf(line_size,'%*9c %d %*9c %d',2);
f.mpsi = temp(1);
f.mtheta=temp(2);
temp=sscanf(line_psi,'%*16c %e %*16c %e',2);
f.psi_edge = temp(1);
f.psitor_edge = temp(2);
temp=fscanf(fid,'%e',[6 f.mpsi+1]);
f.psi = temp(1,:);
f.psitor = temp(2,:);
f.press = temp(3,:);
f.q = temp(4,:);
f.iota = 1./f.q;
f.g = temp(5,:);
f.I = temp(6,:);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
temp = fscanf(fid,'%e',[8 (f.mpsi+1)*(f.mtheta+1)]);
f.psi2d = reshape(temp(1,:),[f.mtheta+1 f.mpsi+1]);
f.theta = reshape(temp(2,:),[f.mtheta+1 f.mpsi+1]);
f.r = reshape(temp(3,:),[f.mtheta+1 f.mpsi+1]);
f.z = reshape(temp(4,:),[f.mtheta+1 f.mpsi+1]);
f.eta = reshape(temp(5,:),[f.mtheta+1 f.mpsi+1]);
f.dphi = reshape(temp(6,:),[f.mtheta+1 f.mpsi+1]);
f.jac = reshape(temp(7,:),[f.mtheta+1 f.mpsi+1]);
f.b0 = reshape(temp(8,:),[f.mtheta+1 f.mpsi+1]);
f.type = 'DCON_EQUIL';
fclose(fid);
return;
end

function f=read_ipec_pmodb(filename)
fid=fopen(filename,'r+');
line=fgetl(fid);
line=fgetl(fid);
line_jac=fgetl(fid);
line_size=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
dex=strfind(line_jac,'=');
f.jac_out = strtrim(line_jac(dex+1:end));
temp=sscanf(line_size,'%*13c %d %*7c %d %*9c %d',3);
f.mstep = temp(1);
f.mpert = temp(2);
f.mthsurf = temp(3);
temp=fscanf(fid,'%e %d %e %e %e %e %e %e %e %e',[10 f.mpert*f.mstep]);
f.psi2d   = reshape(temp(1,:),[f.mpert f.mstep]);
f.xm      = temp(2,1:f.mpert);
f.eulbmnc = reshape(temp(3,:),[f.mpert f.mstep]);
f.eulbmns = reshape(temp(4,:),[f.mpert f.mstep]);
f.lagbmnc = reshape(temp(5,:),[f.mpert f.mstep]);
f.lagbmns = reshape(temp(6,:),[f.mpert f.mstep]);
f.Bdivxprpmnc = reshape(temp(7,:),[f.mpert f.mstep]);
f.Bdivxprpmns = reshape(temp(8,:),[f.mpert f.mstep]);
f.Bkxprpmnc = reshape(temp(9,:),[f.mpert f.mstep]);
f.Bkxprpmns = reshape(temp(10,:),[f.mpert f.mstep]);
f.xn      = str2double(filename(13)).*ones(1,f.mpert);
f.type = 'IPEC_PMODB';
fclose(fid);
return;
end

function f=read_ipec_pmodb_fun(filename)
fid=fopen(filename,'r+');
line=fgetl(fid);
line=fgetl(fid);
line_type=fgetl(fid);
line_size=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
temp=sscanf(line_type,'%*13c %s',1);
f.jac_out = temp;
temp=sscanf(line_type,'%*37c %e',1);
f.R0 = temp;
temp=sscanf(line_size,'%*13c %d %*7c %d %*9c %d',3);
f.mstep = temp(1);
f.mpert = temp(2);
f.mthsurf = temp(3);
temp=fscanf(fid,'%e',[15 f.mstep*f.mthsurf]);
f.psi2d   = reshape(temp(1,:),[f.mthsurf f.mstep]);
f.theta   = reshape(temp(2,:),[f.mthsurf f.mstep]);
f.r   = reshape(temp(3,:),[f.mthsurf f.mstep]);
f.z   = reshape(temp(4,:),[f.mthsurf f.mstep]);
f.eulb_r = reshape(temp(5,:),[f.mthsurf f.mstep]);
f.eulb_i = reshape(temp(6,:),[f.mthsurf f.mstep]);
f.lagb_r = reshape(temp(7,:),[f.mthsurf f.mstep]);
f.lagb_i = reshape(temp(8,:),[f.mthsurf f.mstep]);
f.Bdivxprp_r = reshape(temp(9,:),[f.mthsurf f.mstep]);
f.Bdivxprp_i = reshape(temp(10,:),[f.mthsurf f.mstep]);
f.Bkxprp_r = reshape(temp(11,:),[f.mthsurf f.mstep]);
f.Bkxprp_i = reshape(temp(12,:),[f.mthsurf f.mstep]);
f.equilb = reshape(temp(13,:),[f.mthsurf f.mstep]);
f.dequilbdpsi = reshape(temp(14,:),[f.mthsurf f.mstep]);
f.dequilbdtheta = reshape(temp(15,:),[f.mthsurf f.mstep]);
f.type = 'IPEC_PMODB_FUN';
fclose(fid);
return;
end

function f=read_ipec_pent(filename)
fid=fopen(filename,'r+');
line=fgetl(fid);
line=fgetl(fid);
line_R0=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
line_totals=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
temp=sscanf(line_R0,'%*5c %e %*5c %e %*6c %e %*3c %d',4);
f.R0=temp(1);
f.B0=temp(2);
f.chi1=temp(3);
f.n=temp(4);
temp=sscanf(line_totals,'%*30c %e %*19c %e',2);
f.T_phi_total=temp(1);
f.total_2nddeltaW=temp(2);
temp=fscanf(fid,'%e',[6 inf]);
f.psi=temp(1,:);
f.T_phi=temp(2,:);
f.secdeltaW=temp(3,:);
f.int_T_phi=temp(4,:);
f.int_secdeltaW=temp(5,:);
f.dvdpsi=temp(6,:);
f.type = 'PENT';
fclose(fid);
return;
end

function f=read_ipec_xbnormal(filename)
fid=fopen(filename,'r+');
line=fgetl(fid);
line=fgetl(fid);
line_jac=fgetl(fid);
line_size=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
dex=strfind(line_jac,'=');
f.jac_out = strtrim(line_jac(dex+1:end));
temp=sscanf(line_size,'%*13c %d %*7c %d %*9c %d',3);
f.mstep = temp(1);
f.mpert = temp(2);
f.mthsurf = temp(3);
temp=fscanf(fid,'%e %e %d %e %e %e %e %e %e',[9 f.mpert*f.mstep]);
f.psi2d   = reshape(temp(1,:),[f.mpert f.mstep]);
f.q2d     = reshape(temp(2,:),[f.mpert f.mstep]);
f.xm      = temp(3,1:f.mpert);
f.dmnc     = reshape(temp(4,:),[f.mpert f.mstep]);
f.dmns     = reshape(temp(5,:),[f.mpert f.mstep]);
f.bnmnc     = reshape(temp(6,:),[f.mpert f.mstep]);
f.bnmns     = reshape(temp(7,:),[f.mpert f.mstep]);
f.bwpmnc     = reshape(temp(8,:),[f.mpert f.mstep]);
f.bwpmns     = reshape(temp(8,:),[f.mpert f.mstep]);
f.xn      = str2double(filename(16)).*ones(1,f.mpert);
f.type = 'IPEC_XBNORMAL';
fclose(fid);
return;
end

function f=read_ipec_xbnormal_fun(filename)
fid=fopen(filename,'r+');
line=fgetl(fid);
line=fgetl(fid);
line_jac=fgetl(fid);
line_size=fgetl(fid);
line=fgetl(fid);
line=fgetl(fid);
dex=strfind(line_jac,'=');
f.jac_out = strtrim(line_jac(dex+1:end));
temp=sscanf(line_size,'%*13c %d %*9c %d',2);
f.mstep = temp(1);
f.mthsurf = temp(2);
temp=fscanf(fid,'%e',[10 f.mstep*f.mthsurf]);
f.psi2d   = reshape(temp(1,:),[f.mthsurf f.mstep]);
f.theta   = reshape(temp(2,:),[f.mthsurf f.mstep]);
f.r   = reshape(temp(3,:),[f.mthsurf f.mstep]);
f.z   = reshape(temp(4,:),[f.mthsurf f.mstep]);
f.rvec = reshape(temp(5,:),[f.mthsurf f.mstep]);
f.zvec = reshape(temp(6,:),[f.mthsurf f.mstep]);
f.dxr = reshape(temp(7,:),[f.mthsurf f.mstep]);
f.dxi = reshape(temp(8,:),[f.mthsurf f.mstep]);
f.bnr = reshape(temp(9,:),[f.mthsurf f.mstep]);
f.bni = reshape(temp(10,:),[f.mthsurf f.mstep]);
f.type = 'IPEC_XBNORMAL_FUN';
fclose(fid);
return;
end


