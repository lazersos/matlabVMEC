function data =  vmec_namelist_init(namelist,varargin)
%VMEC_NAMELIST_INIT(namelist) Return initialized namelist structure.
%   VMEC_NAMELIST_INIT(namelist) Returns and initizlied namelst structure.
%   The namelist input is the name of the initizlied namelist you'd like to
%   like.  Available options include:
%   'indata'
%
%   Example:
%       piesinput=pies_namelist_init('input');
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           8/03/11

switch strtrim(upper(namelist))
    case 'INDATA'
        data.datatype = 'VMEC_input';
        data.nfp=0;
        data.ncurr=0;
        data.nsin=0;
        data.niter=0;
        data.nstep=0;
        data.nvacskip=0;
        data.mpol=0;
        data.ntor=0;
        data.ntheta=0;
        data.nzeta=0;
        data.mfilter_fbdy=0;
        data.nfilter_fbdy=0;
        data.ns_array=zeros(1,100);
        data.niter_array=zeros(1,100);
        data.imse=0;
        data.itse=0;
        data.ipnodes=0;
        data.iopt_raxis=0;
        data.imatch_phiedge=0;
        data.nflxs=0;
        data.nbfld=zeros(1,5);
        data.indexflx=zeros(1,5);
        data.rbc=zeros(123,61);
        data.zbs=zeros(123,61);
        data.rbs=zeros(123,61);
        data.zbc=zeros(123,61);
        data.time_slice=0.0;
        data.curtor=0.0;
        data.delt=0.0;
        data.ftol=0.0;
        data.tcon0=0.0;
        data.gamma=0.0;
        data.phiedge=0.0;
        data.spres_ped=0.0;
        data.bloat=0.0;
        data.pres_scale=0.0;
        data.am=zeros(1,21);
        data.ai=zeros(1,21);
        data.ac=zeros(1,21);
        data.aphi=zeros(1,20);
        data.pcurr_type='power_series';
        data.piota_type='power_series';
        data.pmass_type='power_series';
        data.ph_type = 'power_series';
        data.pt_type = 'power_series';
        data.ai_aux_s=zeros(1,100);
        data.ai_aux_f=zeros(1,100);
        data.am_aux_s=zeros(1,100);
        data.am_aux_f=zeros(1,100);
        data.ac_aux_s=zeros(1,100);
        data.ac_aux_f=zeros(1,100);
        data.ah=zeros(1,21);
        data.at=zeros(1,21);
        data.bcrit=0.0;
        data.raxis=zeros(1,102);
        data.zaxis=zeros(1,102);
        data.raxis_cc=zeros(1,102);
        data.raxis_cs=zeros(1,102);
        data.zaxis_cc=zeros(1,102);
        data.zaxis_cs=zeros(1,102);
        data.folt_array=zeros(1,100);
        data.extcur=zeros(1,100);
        data.lpofr=0;
        data.lmac=0;
        data.lfreeb=0;
        data.lrecon=0;
        data.loldout=0;
        data.ledge_dump=0;
        data.lasym=0;
        data.lforbal=0;
        data.lrfp=0;
        data.l_v3fit=0;
        data.lspectrum_dump=0;
        data.loptim=0;
        data.mgrid_file='';
        data.precon_type='none';
        data.arg1='';
    otherwise
        data=-1;
        disp(['  ERROR:  Unknown namelist, ' strtrim(namelist)]);
        return
end
return
end