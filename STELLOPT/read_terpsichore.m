function data = read_terpsichore( filename )
%READ_TERPSICHORE(file) Reads the TERPSICHORE output files.
%   This function reads the files output by the TERPSICHORE code into the
%   Matlab workspace.  It can handle either files name 'fort.XX' or those
%   output by the STELLOPT interface 'terpsichore_XX'  the filenames
%   determine how the data is read so if filenames are changed they must
%   not alter the STELLOPT naming convention.
%
%   Example usage
%       data_16=read_terpsichore('fort.16');     % FORT.16 file
%       data_22=read_terpsichore('fort.22');     % FORT.22 file
%       data_23=read_terpsichore('terpsichore_22_test_file');     % 23 file
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0


% Figure out which file we're reading
if ~isempty(strfind(filename,'terpsichore_16')) || ~isempty(strfind(filename,'fort.16'))
    data = read_terp16(filename);
elseif ~isempty(strfind(filename,'terpsichore_17')) || ~isempty(strfind(filename,'fort.17'))
elseif ~isempty(strfind(filename,'terpsichore_19')) || ~isempty(strfind(filename,'fort.19'))
    data = read_terp19(filename);
elseif ~isempty(strfind(filename,'terpsichore_22')) || ~isempty(strfind(filename,'fort.22'))
    data = read_terp22(filename);
elseif ~isempty(strfind(filename,'terpsichore_23')) || ~isempty(strfind(filename,'fort.23'))
    data = read_terp23(filename);
else
    disp([' Unknown file type: ' filename]);
    data = -1;
end

return;

end

function data=read_terp16(filename)
fid = fopen(filename,'r');
line = fgetl(fid);
while ~feof(fid)
    if strfind(line,'NIM   IVAC   NJ    NK    MM  NMIN  NMAX')
        fgetl(fid);line=fgetl(fid);
        %[data.ni,data.ivac,data.nj,data.nk,data.m,data.nmin,data.nmax] = sscanf(line,'%i',7);
        temp = sscanf(line,'%i',7);
        data.ni = temp(1);
        data.ivac = temp(2);
        data.nj = temp(3);
        data.nk = temp(4);
        data.njk = data.nj*data.nk;
        data.nvi = data.ni+data.ivac;
        data.mm = temp(5);
        data.nmin = temp(6);
        data.nmax = temp(7);
    elseif strfind(line,'TABLE OF R/Z (PR/Z) COEFFICIENTS')
        fgetl(fid);fgetl(fid);fgetl(fid);
        for i=data.nmin:data.nmax
            line = fgetl(fid);
            temp = sscanf(line,'%i',38);
            data.rz_modes(:,i+abs(data.nmin)+1) = temp(1:37);
        end
        data.mb = 0:36;
        data.nb = data.nmin:data.nmax;
    elseif strfind(line,'IOTA           M           VP          P')
        for i = 1:data.ni+1
            line = fgetl(fid);
            temp = sscanf(line,'%e',4);
            data.iota(i) = temp(1);
            data.M(i) = temp(2);
            data.VP(i) = temp(3);
            data.P(i) = temp(4);
        end
    elseif strfind(line,' LMNV')
        temp = sscanf(line,'%*5c %i',3);
        data.lmnv = temp(1);
        data.lmnl = temp(2);
        data.lmnb = temp(3);
    elseif strfind(line,'  mwall  nwall      frwall      fzwall')
        line = fgetl(fid);
        i=1;
        while ~isempty(line)
            temp = sscanf(line,'%i %i %e %e',4);
            data.mwall(i) = temp(1);
            data.nwall(i) = temp(2);
            data.frwall(i) = temp(3);
            data.fzwall(i) = temp(4);
            line = fgetl(fid); i=i+1;
        end
    elseif strfind(line,'   I         VVP(I)     PVP(I)    PVPI(I)')
        for i =1:data.ni-1
            temp =fscanf(fid,'%*i %e %e %e %e %e %e %e',7);
            data.vvp(i) = temp(1);
            data.pvp(i) = temp(2);
            data.pvpi(i) = temp(3);
            data.parpvi(i) = temp(4);
            data.cjvp(i) = temp(5);
            data.civp(i) = temp(6);
            data.equi(i) = temp(7);
        end
    elseif strfind(line,'Modes in the Boozer Table in which')
        fgetl(fid);
        fgetl(fid);
        for i=data.nmin:data.nmax
            line = fgetl(fid);
            temp = sscanf(line,'%i',38);
            data.lrz_modes(:,i+abs(data.nmin)+1) = temp(1:37);
        end
    elseif strfind(line,'   NSTA     MMS   NSMIN   NSMAX   MODEL')
        fgetl(fid);
        temp = fscanf(fid,'%i',5);
        data.nsta=temp(1);
        data.mms=temp(2);
        data.nsmin=temp(3);
        data.nsmax=temp(4);
        data.model=temp(5);
        data.ms=0:55;
        data.ns=data.nsmin:data.nsmax;
    elseif strfind(line,'TABLE OF R AND Z COEFFICIENTS FOR STABILITY')
        fgetl(fid);fgetl(fid);fgetl(fid);
        for i=data.ns
            line = fgetl(fid);
            temp = sscanf(line,'%i',57);
            data.rzs_modes(:,i+abs(data.nsmin)+1) = temp(1:56);
        end
%      elseif strfind(line,'FOURIN')
%          n = sscanf(line,'%*s %i',1);
%          for i=1:n
%              fgetl(fid);
%              temp = fscanf(fid,'%e %e %e %e',4)
%              data.c1(1,1,i) = temp(1);
%              data.c1(1,2,i) = temp(2);
%              data.c1(2,1,i) = temp(3);
%              data.c1(2,2,i) = temp(4);
%              temp = fscanf(fid,'%e %e %e %e',4);
%              data.c2(1,1,i) = temp(1);
%              data.c2(1,2,i) = temp(2);
%              data.c2(2,1,i) = temp(3);
%              data.c2(2,2,i) = temp(4);
%              temp = fscanf(fid,'%e %e %e %e',4);
%              data.c3(1,1,i) = temp(1);
%              data.c3(1,2,i) = temp(2);
%              data.c3(2,1,i) = temp(3);
%              data.c3(2,2,i) = temp(4);
%              temp = fscanf(fid,'%e %e %e %e',4);
%              data.c4(1,1,i) = temp(1);
%              data.c4(1,2,i) = temp(2);
%              data.c4(2,1,i) = temp(3);
%              data.c4(2,2,i) = temp(4);
%              temp = fscanf(fid,'%e %e %e %e',4);
%              data.c5(1,1,i) = temp(1);
%              data.c5(1,2,i) = temp(2);
%              data.c5(2,1,i) = temp(3);
%              data.c5(2,2,i) = temp(4);
%              temp = fscanf(fid,'%e %e %e %e',4);
%              data.c6(1,1,i) = temp(1);
%              data.c6(1,2,i) = temp(2);
%              data.c6(2,1,i) = temp(3);
%              data.c6(2,2,i) = temp(4);
%              fgetl(fid); fgetl(fid);
%          end
%     elseif strfind(line,'LHSMAT:')
%         fgetl(fid);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.a(1,1,n)=temp(1);
%         data.a(1,2,n)=temp(2);
%         data.a(2,1,n)=temp(3);
%         data.a(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.c(1,1,n)=temp(1);
%         data.c(1,2,n)=temp(2);
%         data.c(2,1,n)=temp(3);
%         data.c(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.d(1,1,n)=temp(1);
%         data.d(1,2,n)=temp(2);
%         data.d(2,1,n)=temp(3);
%         data.d(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.e(1,1,n)=temp(1);
%         data.e(1,2,n)=temp(2);
%         data.e(2,1,n)=temp(3);
%         data.e(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.f(1,1,n)=temp(1);
%         data.f(1,2,n)=temp(2);
%         data.f(2,1,n)=temp(3);
%         data.f(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.g(1,1,n)=temp(1);
%         data.g(1,2,n)=temp(2);
%         data.g(2,1,n)=temp(3);
%         data.g(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.h(1,1,n)=temp(1);
%         data.h(1,2,n)=temp(2);
%         data.h(2,1,n)=temp(3);
%         data.h(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.bl(1,1,n)=temp(1);
%         data.bl(1,2,n)=temp(2);
%         data.bl(2,1,n)=temp(3);
%         data.bl(2,2,n)=temp(4);
%         temp = fscanf(fid,'%e %e %e %e %i',5);
%         n = temp(5);
%         data.bu(1,1,n)=temp(1);
%         data.bu(1,2,n)=temp(2);
%         data.bu(2,1,n)=temp(3);
%         data.bu(2,2,n)=temp(4);
    elseif strfind(line,'RHSMAT1:')
        fgetl(fid);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bb(1,1,n)=temp(1);
        data.bb(1,2,n)=temp(2);
        data.bb(2,1,n)=temp(3);
        data.bb(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.ba(1,1,n)=temp(1);
        data.ba(1,2,n)=temp(2);
        data.ba(2,1,n)=temp(3);
        data.ba(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bc(1,1,n)=temp(1);
        data.bc(1,2,n)=temp(2);
        data.bc(2,1,n)=temp(3);
        data.bc(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bd(1,1,n)=temp(1);
        data.bd(1,2,n)=temp(2);
        data.bd(2,1,n)=temp(3);
        data.bd(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.be(1,1,n)=temp(1);
        data.be(1,2,n)=temp(2);
        data.be(2,1,n)=temp(3);
        data.be(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bf(1,1,n)=temp(1);
        data.bf(1,2,n)=temp(2);
        data.bf(2,1,n)=temp(3);
        data.bf(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bg(1,1,n)=temp(1);
        data.bg(1,2,n)=temp(2);
        data.bg(2,1,n)=temp(3);
        data.bg(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bh(1,1,n)=temp(1);
        data.bh(1,2,n)=temp(2);
        data.bh(2,1,n)=temp(3);
        data.bh(2,2,n)=temp(4);   
    elseif strfind(line,'LHSMAT1:')
        fgetl(fid);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.ai(1,1,n)=temp(1);
        data.ai(1,2,n)=temp(2);
        data.ai(2,1,n)=temp(3);
        data.ai(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.cc(1,1,n)=temp(1);
        data.cc(1,2,n)=temp(2);
        data.cc(2,1,n)=temp(3);
        data.cc(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.di(1,1,n)=temp(1);
        data.di(1,2,n)=temp(2);
        data.di(2,1,n)=temp(3);
        data.di(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.ei(1,1,n)=temp(1);
        data.ei(1,2,n)=temp(2);
        data.ei(2,1,n)=temp(3);
        data.ei(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.fi(1,1,n)=temp(1);
        data.fi(1,2,n)=temp(2);
        data.fi(2,1,n)=temp(3);
        data.fi(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.gi(1,1,n)=temp(1);
        data.gi(1,2,n)=temp(2);
        data.gi(2,1,n)=temp(3);
        data.gi(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.hi(1,1,n)=temp(1);
        data.hi(1,2,n)=temp(2);
        data.hi(2,1,n)=temp(3);
        data.hi(2,2,n)=temp(4);
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bl1(1,1,n)=temp(1);
        data.bl1(1,2,n)=temp(2);
        data.bl1(2,1,n)=temp(3);
        data.bl1(2,2,n)=temp(4);   
        temp = fscanf(fid,'%e %e %e %e %i',5);
        n = temp(5);
        data.bu1(1,1,n)=temp(1);
        data.bu1(1,2,n)=temp(2);
        data.bu1(2,1,n)=temp(3);
        data.bu1(2,2,n)=temp(4);          
%     elseif strfind(line,'I        Q(SI)        DELTAW(SI)   MERCIER CRITERION')
%         fgetl(fid);fgetl(fid);fgetl(fid);
%         for i = 1:data.nvi
%             line = fgetl(fid);
%             temp = sscanf(line,'%*i %e %e %e %e %e %e',6);
%             data.Q(i) = temp(1);
%             data.DELTAW(i) = temp(2);
%             data.KO_STD(i) = temp(3);
%             data.KO_LAM(i) = temp(4);
%             data.NI_STD(i) = temp(5);
%             data.NI_LAM(i) = temp(6);
%         end
    elseif strfind(line,'EIGENVALUE SHIFT')
        data.eigen_shift = sscanf(line,'%*s %*s %*s %e',1);
    elseif strfind(line,'RAYLEIGH QUOTIENT')
        data.rayl_quot = sscanf(line,'%*s %*s %*s %e',1);
    elseif strfind(line,'EIGENVALUE FROM NORMALIZATION')
        data.eigen_norm = sscanf(line,'%*s %*s %*s %*s %e',1);
    elseif strfind(line,'ITERATIONS DONE')
        data.iter = sscanf(line,'%*s %*s %*s %i',1);
    elseif strfind(line,'NUMBER OF NEGATIVE EIGENVALUES')
        data.neg_eigen = sscanf(line,'%*s %*s %*s %*s %*s %i',1);
    elseif strfind(line,'EIGENVALUE FROM WP / WK')
        try
        temp =  sscanf(line,'%*6c %e %*6c %e %*27c %e',3);
        data.WP = temp(1);
        data.WK = temp(2);
        data.eigen = temp(3);
        catch
        end
    elseif strfind(line,' GROWTH RATE')
        data.growthrate = sscanf(line,'%*15c %e',1);
    elseif strfind(line,'PARVEC')
        line = fgetl(fid);
        temp = sscanf(line,'%i %i %i %e %e',5);
        lc = temp(1); 
        data.ms(lc)    = temp(2);
        data.ns(lc)    = temp(3);
        data.ximax(lc) = temp(4);
        data.etamax(lc) = temp(5);
        data.xi(lc,:) = fscanf(fid,'%e',data.nvi+1);
        data.eta(lc,:) = fscanf(fid,'%e',data.nvi);
    end
    line = fgetl(fid);
end
end

function data=read_terp17(filename)
end

function data=read_terp19(filename)
fid = fopen(filename,'r');
temp = fscanf(fid,'%i',3);
data.njk=temp(1);
data.nj=temp(2);
data.nk=temp(3);
temp = fscanf(fid,'%i %e %e %e %e',[5 data.njk]);
data.idex=temp(1,:);
data.rwall=reshape(temp(2,:),[data.nj data.nk]);
data.zwall=reshape(temp(3,:),[data.nj data.nk]);
data.rpvi=reshape(temp(4,:),[data.nj data.nk]);
data.zpvi=reshape(temp(5,:),[data.nj data.nk]);
fclose(fid);
end

function data=read_terp22(filename)
fid = fopen(filename,'r');
temp = fscanf(fid,'%i %i %e');
data.ni = temp(1);
data.lmns = temp(2);
data.omega2 = temp(3);
data.pth = fscanf(fid,'%e',data.ni);
data.aiota = fscanf(fid,'%e',data.ni);
data.wpsi = fscanf(fid,'%e',data.ni);
data.ms   = fscanf(fid,'%i',data.lmns);
data.ns   = fscanf(fid,'%i',data.lmns);
data.am = fscanf(fid,'%e',data.ni);
data.pvp = fscanf(fid,'%e',data.ni);
data.pvpi = fscanf(fid,'%e',data.ni);
fclose(fid);
end

function data=read_terp23(filename)
int_type = 'int16';
flt_type = 'single';
fid = fopen(filename);
temp = fread(fid,16,int_type);
data.ni = temp(3);
data.nj = temp(5);
data.nk = temp(7);
data.njk = data.nj.*data.nk;
data.nsta = temp(9);
data.nper = temp(11);
data.lmns = temp(13);
data.nvi = temp(15);
temp = fread(fid,data.nvi+3,flt_type);
data.s = temp(3:end); % 0:nvi
data.pth = fread(fid,data.ni,flt_type);
data.fpp = fread(fid,data.ni+1,flt_type);
data.ftp = fread(fid,data.ni+1,flt_type);
data.cj = fread(fid,data.ni,flt_type);
data.ci = fread(fid,data.ni,flt_type);
data.pp = fread(fid,data.ni,flt_type);
data.parity = fread(fid,1,flt_type);
temp = fread(fid,3,int_type);
temp = fread(fid,2.*data.lmns,int_type);
data.ms = temp(2:2:end);
temp = fread(fid,2.*data.lmns,int_type);
data.ns = temp(2:2:end);
data.ql = fread(fid,data.lmns,flt_type);
temp = fread(fid,1,int_type);
data.xi = fread(fid,data.lmns*(data.nvi+1),flt_type);
data.eta = fread(fid,data.lmns*(data.nvi),flt_type);
data.rmu = fread(fid,data.lmns*(data.nvi),flt_type);
temp = fread(fid,4,int_type);
data.r = permute(reshape(fread(fid,data.njk.*(data.ni+1),flt_type),[data.nj data.nk data.ni+1]),[3 1 2]);
data.z = permute(reshape(fread(fid,data.njk.*(data.ni+1),flt_type),[data.nj data.nk data.ni+1]),[3 1 2]);
data.phv = permute(reshape(fread(fid,data.njk.*(data.ni),flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.rs = permute(reshape(fread(fid,data.njk.*(data.ni),flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.zs = permute(reshape(fread(fid,data.njk.*(data.ni),flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.rt = permute(reshape(fread(fid,data.njk.*(data.ni),flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.zt = permute(reshape(fread(fid,data.njk.*(data.ni),flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.rp = permute(reshape(fread(fid,data.njk.*(data.ni),flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.zp = permute(reshape(fread(fid,data.njk.*(data.ni),flt_type),[data.nj data.nk data.ni]),[3 1 2]);
temp = fread(fid,4,int_type);
data.bjac = permute(reshape(fread(fid,data.njk.*(data.nvi+1),flt_type),[data.nj data.nk data.nvi+1]),[3 1 2]);
data.sigbs = permute(reshape(fread(fid,data.njk.*(data.ni+1),flt_type),[data.nj data.nk data.ni+1]),[3 1 2]);
data.gssl = permute(reshape(fread(fid,data.njk.*(data.nvi+1),flt_type),[data.nj data.nk data.nvi+1]),[3 1 2]);
data.gstl = permute(reshape(fread(fid,data.njk.*(data.nvi+1),flt_type),[data.nj data.nk data.nvi+1]),[3 1 2]);
data.gttl = permute(reshape(fread(fid,data.njk.*(data.nvi+1),flt_type),[data.nj data.nk data.nvi+1]),[3 1 2]);
temp = fread(fid,4,int_type);
data.lmnb = fread(fid,1,int_type);
temp = fread(fid,2.*data.lmnb,int_type);
data.mb = temp(2:2:end);
temp = fread(fid,2.*data.lmnb,int_type);
data.nb = temp(2:2:end);
temp = fread(fid,1,int_type);
data.fbjac = fread(fid,data.lmnb*(data.ni+1),flt_type);
data.bjacs = permute(reshape(fread(fid,data.njk.*(data.ni+1),flt_type),[data.nj data.nk data.ni+1]),[3 1 2]);
data.ftpp = fread(fid,data.ni+1,flt_type);
data.fppp = fread(fid,data.ni+1,flt_type);
data.cip = fread(fid,data.ni,flt_type);
data.cjp = fread(fid,data.ni,flt_type);
data.WP = fread(fid,1,flt_type);
data.WK = fread(fid,1,flt_type);
temp = fread(fid,4,int_type);
data.gparp = permute(reshape(fread(fid,data.njk.*data.ni,flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.gperp = permute(reshape(fread(fid,data.njk.*data.ni,flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.sigmab = permute(reshape(fread(fid,data.njk.*data.ni,flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.taub = permute(reshape(fread(fid,data.njk.*data.ni,flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.parkur = permute(reshape(fread(fid,data.njk.*data.ni,flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.parjp = permute(reshape(fread(fid,data.njk.*data.ni,flt_type),[data.nj data.nk data.ni]),[3 1 2]);
data.curfac = fread(fid,1,flt_type);
data.modelk = fread(fid,1,flt_type);

fclose(fid);
end