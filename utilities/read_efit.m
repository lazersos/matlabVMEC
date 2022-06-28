function efit_data = read_efit(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

filetype='gfile';
rotate = 1;

temp=strfind(filename,'/');
if ~isempty(temp)
    dex=temp(length(temp));
    filename_short=strtrim(filename(dex+1:length(filename)));
else
    filename_short=strtrim(filename);
end

switch filename_short(1)
    case 'g'
        filetype='gfile';
    case 'a'
        filetype='afile';
    case 'q'
        filetype='qfile'; %For trgui exported eqdsk files (CLISTE)
    otherwise
        if strcmp(filename_short(end-5:end-1),'eqdsk') || strcmp(filename_short(end-6:end-2),'eqdsk')
            filetype='gfile';
        else
        disp(['Unknown filetype: ' filename_short]);
        efit_data=-1;
        return
        end
end

switch filetype
    case 'gfile'
        fid=fopen(filename,'r');
        header=fgetl(fid);
        efit_data.type=sscanf(header,'%10c %*s',1);
        efit_data.date=sscanf(header,'%*10c %10c %*s',1);
        efit_data.shot=sscanf(header,'%*10c %*10c %9c %*s',1);
        efit_data.time=sscanf(header,'%*10c %*10c %*9c %10c %*s',1);
        temp=sscanf(header,'%*49c %d %d %d');
        efit_data.ipest=temp(1);
        nx=temp(2);
        nz=temp(3);
        efit_data.nx=nx;
        efit_data.nz=nz;
        efit_data.xdim=fscanf(fid,'%e',1);
        efit_data.zdim=fscanf(fid,'%e',1);
        efit_data.zc=fscanf(fid,'%e',1);
        efit_data.redge=fscanf(fid,'%e',1);
        efit_data.zmid=fscanf(fid,'%e',1);
        efit_data.xaxis=fscanf(fid,'%e',1);
        efit_data.zaxis=fscanf(fid,'%e',1);
        efit_data.psiaxis=fscanf(fid,'%e',1);
        efit_data.psilim=fscanf(fid,'%e',1);
        efit_data.btor=fscanf(fid,'%e',1);
        efit_data.totcur=fscanf(fid,'%e',1);
        efit_data.psimx(1)=fscanf(fid,'%e',1);
        efit_data.psimx(2)=fscanf(fid,'%e',1);
        efit_data.xax(1)=fscanf(fid,'%e',1);
        efit_data.xax(2)=fscanf(fid,'%e',1);
        efit_data.zax(1)=fscanf(fid,'%e',1);
        efit_data.zax(2)=fscanf(fid,'%e',1);
        efit_data.psisep=fscanf(fid,'%e',1);
        efit_data.xsep=fscanf(fid,'%e',1);
        efit_data.zsep=fscanf(fid,'%e',1);
        efit_data.sf=fscanf(fid,'%e',nx);
        efit_data.sp=fscanf(fid,'%e',nx);
        efit_data.sffp=fscanf(fid,'%e',nx);
        efit_data.spp=fscanf(fid,'%e',nx);
        efit_data.psixz=fscanf(fid,'%e',[nx nz]);
        efit_data.qpsi=fscanf(fid,'%e',nx);
        efit_data.nbndry=fscanf(fid,'%e',1);
        efit_data.nlim=fscanf(fid,'%e',1);
        bdry_temp=fscanf(fid,'%e',efit_data.nbndry*2);
        efit_data.xbndry=bdry_temp(1:2:efit_data.nbndry*2);
        efit_data.zbndry=bdry_temp(2:2:efit_data.nbndry*2);
        lim_temp=fscanf(fid,'%e',efit_data.nlim*2);
        efit_data.xlim=lim_temp(1:2:efit_data.nlim*2);
        efit_data.zlim=lim_temp(2:2:efit_data.nlim*2);
        if rotate
            efit_data.kvtor=fscanf(fid,'%d',1);
            efit_data.rvtor=fscanf(fid,'%e',1);
            efit_data.nmass=fscanf(fid,'%d',1);
            if efit_data.kvtor > 0
                efit_data.pressw=fscanf(fid,'%e',nx);
                efit_data.pwprim=fscanf(fid,'%e',nx);
            end
            if efit_data.nmass > 0
                efit_data.rho0=fscanf(fid,'%e',nx);
            end
        end
        fclose(fid);
        efit_data.out=read_namelist(filename,'OUT1');
        efit_data.basis=read_namelist(filename,'BASIS');
        efit_data.chiout=read_namelist(filename,'CHIOUT');
        efit_data.datatype='EFIT_G';
        efit_data.xgrid=efit_data.redge+efit_data.xdim.*(0:(nx-1))./(nx-1);
        efit_data.zgrid=efit_data.zdim.*(0:(nz-1))./(nz-1)-efit_data.zdim/2.;
    case 'afile'
    case 'qfile'
        fid=fopen(filename,'r');
        header=fgetl(fid);
%        header2 = textscan(header,'%s %s %s %s %s %d %d %d',1);
%         efit_data.type=header2{1}{1};
%         efit_data.date=datestr(now);
%         efit_data.shot=header2{2}{1};
%         efit_data.time=header2{3}{1};
%         efit_data.ipest=header2{4}{1};
%         nx=header2{6};
%         nz=header2{7};
        efit_data.type=sscanf(header,'%10c %*s',1);
        efit_data.date=datestr(now);
        efit_data.shot=sscanf(header,'%*10c %8c %*s',1);
        efit_data.time=sscanf(header,'%*10c %*8c %6c %*s',1);
        temp=sscanf(header,'%*35c %d %d %d');
        efit_data.ipest=temp(1);
        nx=temp(2);
        nz=temp(3);
        efit_data.nx=nx;
        efit_data.nz=nz;
        efit_data.xdim=fscanf(fid,'%e',1);
        efit_data.zdim=fscanf(fid,'%e',1);
        efit_data.zc=fscanf(fid,'%e',1);
        efit_data.redge=fscanf(fid,'%e',1);
        efit_data.zmid=fscanf(fid,'%e',1);
        efit_data.xaxis=fscanf(fid,'%e',1);
        efit_data.zaxis=fscanf(fid,'%e',1);
        efit_data.psiaxis=fscanf(fid,'%e',1);
        efit_data.psilim=fscanf(fid,'%e',1);
        efit_data.btor=fscanf(fid,'%e',1);
        efit_data.totcur=fscanf(fid,'%e',1);
        efit_data.psimx(1)=fscanf(fid,'%e',1);
        efit_data.psimx(2)=fscanf(fid,'%e',1);
        efit_data.xax(1)=fscanf(fid,'%e',1);
        efit_data.xax(2)=fscanf(fid,'%e',1);
        efit_data.zax(1)=fscanf(fid,'%e',1);
        efit_data.zax(2)=fscanf(fid,'%e',1);
        efit_data.psisep=fscanf(fid,'%e',1);
        efit_data.xsep=fscanf(fid,'%e',1);
        efit_data.zsep=fscanf(fid,'%e',1);
        efit_data.sf=fscanf(fid,'%e',nx);
        efit_data.sp=fscanf(fid,'%e',nx);
        efit_data.sffp=fscanf(fid,'%e',nx);
        efit_data.spp=fscanf(fid,'%e',nx);
        efit_data.psixz=fscanf(fid,'%e',[nx nz]);
        efit_data.qpsi=fscanf(fid,'%e',nx);
        efit_data.nbndry=fscanf(fid,'%e',1);
        efit_data.nlim=fscanf(fid,'%e',1);
        bdry_temp=fscanf(fid,'%e',efit_data.nbndry*2);
        efit_data.xbndry=bdry_temp(1:2:efit_data.nbndry*2);
        efit_data.zbndry=bdry_temp(2:2:efit_data.nbndry*2);
        lim_temp=fscanf(fid,'%e',efit_data.nlim*2);
        efit_data.xlim=lim_temp(1:2:efit_data.nlim*2);
        efit_data.zlim=lim_temp(2:2:efit_data.nlim*2);
        if rotate
            efit_data.kvtor=fscanf(fid,'%d',1);
            efit_data.rvtor=fscanf(fid,'%e',1);
            efit_data.nmass=fscanf(fid,'%d',1);
            if efit_data.kvtor > 0
                efit_data.pressw=fscanf(fid,'%e',nx);
                efit_data.pwprim=fscanf(fid,'%e',nx);
            end
            if efit_data.nmass > 0
                efit_data.rho0=fscanf(fid,'%e',nx);
            end
        end
        fclose(fid);
        efit_data.out=read_namelist(filename,'OUT1');
        efit_data.basis=read_namelist(filename,'BASIS');
        efit_data.chiout=read_namelist(filename,'CHIOUT');
        efit_data.datatype='EFIT_G';
        efit_data.xgrid=efit_data.redge+efit_data.xdim.*(0:(nx-1))./(nx-1);
        efit_data.zgrid=efit_data.zdim.*(0:(nz-1))./(nz-1)-efit_data.zdim/2.;
end

end

