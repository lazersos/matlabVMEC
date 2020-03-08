function f = read_siesta(filename)
%READ_SIESTA This function reads the SIESTA netCDF file.
%   This subroutine read the netCDF file produced by SIESTA.  It outputs
%   arrays in into a similar format to that of read_vmec.
%
% Example usage
%      data=read_siesta('siesta_test.nc');     % Reads VMEC wout file
%
% Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
% Version:       1.00


% Defaults
f=[];
% Number of Arguments
if (nargin == 0)
    disp('read_siesta requires a filename.');
    return
end
f=read_netcdf(filename);

% Create helpers
ds = 1./double(f.nrad-1);

% Reformat to VMEC 
f.nfp = double(f.nfp);
f.mpol = double(f.mpol);
f.ntor = double(f.ntor);
f.ns = double(f.nrad);
f.input_extension = filename;
f.Itor = f.curtor;
f.Rmajor = f.rmajor;
f.iotaf = f.chipf_r_'./f.phipf_r_';
f.phipf = f.phipf_r_';
f.chipf = f.chipf_r_';
f.phif  = cumsum(f.phipf);
f.phi   = f.phif;

% Define xm and xn
% Siesta uses the NESCOIL kernel (mu+nv)
mn = 1;
f.xm=[];
f.xn=[];
for m=0:f.mpol
   for n = -f.ntor:f.ntor
       f.xm(mn) = m;
       f.xn(mn) = n.*f.nfp;
       mn = mn+1;
   end
end
f.mnmax = mn-1;

% Now transform array quantities
shape = [f.ns (2*f.ntor+1)*(f.mpol+1)];
if isfield(f,'rmnc_m_n_r_')
    f.rmnc=reshape(f.rmnc_m_n_r_,shape)';
    temp=[zeros(f.mnmax,1) diff(f.rmnc,1,2)]./ds;
    f.rsmnc=half2full(temp);
    f.rsmnc(:,f.ns) = f.rsmnc(:,f.ns-1);
    f.rsmnc(:,1) = temp(:,2).*2.*ds;
    f.rsmnc(f.xm~=1,1)=0;
    f.rumns=-f.rmnc.*repmat(f.xm',[1 f.ns]);
    f.rvmns=-f.rmnc.*repmat(f.xn',[1 f.ns]);
end
if isfield(f,'zmns_m_n_r_')
    f.zmns=reshape(f.zmns_m_n_r_,shape)';
    temp=[zeros(f.mnmax,1) diff(f.zmns,1,2)]./ds;
    f.zsmns=half2full(temp);
    f.zsmns(:,f.ns) = f.zsmns(:,f.ns-1);
    f.zsmns(:,1) = temp(:,2).*2.*ds;
    f.zsmns(f.xm~=1,1)=0;
    f.zumnc=f.zmns.*repmat(f.xm',[1 f.ns]);
    f.zvmnc=f.zmns.*repmat(f.xn',[1 f.ns]);
end
if isfield(f,'pmnch_m_n_r_')
    temp=reshape(f.pmnch_m_n_r_,shape)';
    f.pmnc=half2full(temp);
end
% Not sure if handled yet
%if isfield(f,'rmns_m_n_r_')
%    f.rmns=reshape(f.rmns_m_n_r_,shape)';
%end
%if isfield(f,'zmnc_m_n_r_')
%    f.zmnc=reshape(f.zmnc_m_n_r_,shape)';
%end
if isfield(f,'bsupsmnsh_m_n_r_')
    temp=reshape(f.bsupsmnsh_m_n_r_,shape)';
    f.bsupsmns=half2full(temp);
    f.bsupsmns(f.xm~=1,1) = 0; % Filter a_s_full(1) = a_s_half(2)   For m = 1    and   0. For m ? 1
end
if isfield(f,'bsupumnch_m_n_r_')
    temp=reshape(f.bsupumnch_m_n_r_,shape)';
    f.bsupumnc=half2full(temp);
    f.bsupumnc(:,1) = 0; % Filter a_u_full(1) = 0
end
if isfield(f,'bsupvmnch_m_n_r_')
    temp=reshape(f.bsupvmnch_m_n_r_,shape)';
    f.bsupvmnc=half2full(temp);
    f.bsupvmnc(f.xm~=0,1) = 0; % Filter a_v_full(1) = a_v_half(2)   For m = 0    and  0 for m > 0
end
if isfield(f,'bsubsmnsh_m_n_r_')
    temp=reshape(f.bsubsmnsh_m_n_r_,shape)';
    f.bsubsmns=half2full(temp);
    f.bsubsmns(f.xm~=1,1) = 0; % Filter a_s_full(1) = a_s_half(2)   For m = 1    and   0. For m ? 1
end
if isfield(f,'bsubumnch_m_n_r_')
    temp=reshape(f.bsubumnch_m_n_r_,shape)';
    f.bsubumnc=half2full(temp);
    f.bsubumnc(:,1) = 0; % Filter a_u_full(1) = 0
end
if isfield(f,'bsubvmnch_m_n_r_')
    temp=reshape(f.bsubvmnch_m_n_r_,shape)';
    f.bsubvmnc=half2full(temp);
    f.bsubvmnc(f.xm~=0,1) = 0; % Filter a_v_full(1) = a_v_half(2)   For m = 0    and  0 for m > 0
end
if isfield(f,'jksupsmnsf_m_n_r_')
    temp=reshape(f.jksupsmnsf_m_n_r_,shape)';
    f.currsmns=half2full(temp);
    f.currsmns(f.xm~=1,1) = 0; % Filter a_s_full(1) = a_s_half(2)   For m = 1    and   0. For m ? 1
end
if isfield(f,'jksupumncf_m_n_r_')
    temp=reshape(f.jksupumncf_m_n_r_,shape)';
    f.currumnc=half2full(temp);
    f.currumnc(:,1) = 0; % Filter a_u_full(1) = 0
end
if isfield(f,'jksupvmncf_m_n_r_')
    temp=reshape(f.jksupvmncf_m_n_r_,shape)';
    f.currvmnc=half2full(temp);
    f.currvmnc(f.xm~=0,1) = 0; % Filter a_v_full(1) = a_v_half(2)   For m = 0    and  0 for m > 0
end

f.datatype='siesta';
return;

end

function f = half2full(a)
mn = size(a,1);
ns = size(a,2);
f = zeros(mn,ns);
f(:,1:ns-1) = 0.5.*(a(:,2:ns)+a(:,1:ns-1));
f(:,ns) = 0;
f(:,1)  = a(:,2);
f(:,ns) = a(:,ns);
return;
end

