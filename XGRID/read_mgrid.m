function data=read_mgrid(filename,varargin)
%data=READ_MGRID(filename[,'quiet']) Reads data from mgrid file
% This function reads the mgrid file produced by the MAKEGRID Fortran
% routine.  To supress command window messages, run with 'quiet' option.
%   Options:
%       'quiet':    Supress all non-error messages.
%
%   Usage:
%       mgrid_data=read_mgrid('mgrid.test');
%
%   See also plot_mgrid.
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Date:       7/20/11
%   Version:    1.51

loutput=1;
lnewstyle=0;
% Handle varargin
if nargin > 1
    for i=1:nargin-1
        switch varargin{i}
            case 'quiet'
                loutput=0;
        end
    end
end
% Define Some File Constants
isize='int32';
fsize='float64';
csize='*char';
% Check the filetype
if contains(filename,'.nc')
    data=read_netcdf(filename,'flipdim');
    % First we rename some fields
    data.nr=data.ir;
    data.nz=data.jz;
    data.nphi=data.kp;
    data.curlabel=regexp(data.coil_group,'\s+(,)?','split');
    data=rmfield(data,'ir');
    data=rmfield(data,'jz');
    data=rmfield(data,'coil_group');
    data=rmfield(data,'kp');
    % In the netCDF file each current group has it's own r,phi,z field
    data.br=zeros(data.nr,data.nz,data.nphi,data.nextcur);
    data.bphi=zeros(data.nr,data.nz,data.nphi,data.nextcur);
    data.bz=zeros(data.nr,data.nz,data.nphi,data.nextcur);
    data.bx=zeros(data.nr,data.nz,data.nphi,data.nextcur);
    data.by=zeros(data.nr,data.nz,data.nphi,data.nextcur);
    for i=1:data.nextcur
        data.br(:,:,:,i)=data.(['br_' num2str(i,'%03d')]);
        data.bphi(:,:,:,i)=data.(['bp_' num2str(i,'%03d')]);
        data.bz(:,:,:,i)=data.(['bz_' num2str(i,'%03d')]);
        data=rmfield(data,['br_' num2str(i,'%03d')]);
        data=rmfield(data,['bp_' num2str(i,'%03d')]);
        data=rmfield(data,['bz_' num2str(i,'%03d')]);
    end
    data.phi=0:2*pi/double(data.nfp)/double(data.nphi):2*pi/double(data.nfp);
    data.phi=data.phi(1:data.nphi);
    data.raxis=data.rmin:(data.rmax-data.rmin)/double(data.nr-1):data.rmax;
    data.zaxis=data.zmin:(data.zmax-data.zmin)/double(data.nz-1):data.zmax;
    % Now make the bx and by arrays
    for k=1:data.nphi
        data.bx(:,:,k,:)=data.br(:,:,k,:).*cos(data.phi(k))...
            -data.bphi(:,:,k,:).*sin(data.phi(k));
        data.by(:,:,k,:)=data.br(:,:,k,:).*sin(data.phi(k))...
            +data.bphi(:,:,k,:).*cos(data.phi(k));
    end
    disp('NOTE:  structure elements may not correctly named!');
else
    % Open the File this will handle the new style file where nextcut<0
    fid=fopen(filename,'rb');
    % Read nr,nz,nphi,nfp,nextcur
    fread(fid,1,isize);
    temp=fread(fid,5,isize);
    nr=temp(1);
    nz=temp(2);
    nphi=temp(3);
    nfp=temp(4);
    nextcur=temp(5);
    if nextcur<0
        nextcur=-nextcur;
        lnewstyle=1;
        if loutput, disp(' - New Style binary file detected. (nextcur<0)'); end
    end
    if loutput, disp(strcat(' - Number of Current Systems:',num2str(nextcur))); end
    % Read in limits of domain
    fread(fid,1,fsize);
    temp=fread(fid,4,fsize);
    rmin=temp(1);
    zmin=temp(2);
    rmax=temp(3);
    zmax=temp(4);
    % Read in Current Labels 30 char strings
    curlabel=cell(1,nextcur);
    for i=1:nextcur
        curlabel{i}=strtrim(fread(fid,30,csize))';
    end
    % Now read in B-field data
    nt=nr*nz*nphi;
    if loutput, disp(strcat(' - Number of Gridpoints:',num2str(nt))); end
    % Now handle both file formats
    if lnewstyle
        %Lifted from mgrid_mod.f in the LIBSTELL package
        br=zeros(nr,nz,nphi,nextcur);
        bz=zeros(nr,nz,nphi,nextxur);
        bphi=zeros(nr,nz,nphi,nextcur);
        fread(fid,1,fsize);
        br=fread(fid,[nr nz nphi],fsize);
        bz=fread(fid,[nr nz nphi],fsize);
        bphi=fread(fid,[nr nz nphi],fsize);
        mgrid_mode=fread(fid,1,csize);
        extcur=fread(fid,nextcur,fsize);
    else
        br=zeros(nt,nextcur);
        bz=zeros(nt,nextcur);
        bphi=zeros(nt,nextcur);
        fread(fid,1,fsize);
        for i=1:nextcur
            fread(fid,1,fsize);
            if loutput, disp(strcat(' - Reading Current system:',num2str(i))); end
            temp=fread(fid,3*nt,fsize);
            % So now we need to reformulate the arrays.  Let us first extract
            % the x z and phi data
            br(:,i)=temp(1:3:nt*3);
            bz(:,i)=temp(2:3:nt*3);
            bphi(:,i)=temp(3:3:nt*3);
        end
        % Reformulate the Arrays
        tempx=zeros(nr,nz,nphi,nextcur);
        tempy=zeros(nr,nz,nphi,nextcur);
        tempz=zeros(nr,nz,nphi,nextcur);
        for i=1:nextcur
            tempx(:,:,:,i)=reshape(br(:,i),nr,nz,nphi);
            tempy(:,:,:,i)=reshape(bz(:,i),nr,nz,nphi);
            tempz(:,:,:,i)=reshape(bphi(:,i),nr,nz,nphi);
        end
        br=tempx;
        bz=tempy;
        bphi=tempz;
        mgrid_mode='N';
    end
    % Close the file
    fclose(fid);
    % Create Output Structure
    data.br=br;
    data.bz=bz;
    data.bphi=bphi;
    data.nr=nr;
    data.rmin=rmin;
    data.rmax=rmax;
    data.raxis=rmin:(rmax-rmin)/(nr-1):rmax;
    data.zaxis=zmin:(zmax-zmin)/(nz-1):zmax;
    data.nz=nz;
    data.zmin=zmin;
    data.zmax=zmax;
    data.phi=0:2*pi/(nfp)/(nphi-1):2*pi/nfp;
    data.nfp=nfp;
    data.nphi=nphi;
    data.nextcur=nextcur;
    data.curlabel=curlabel;
    % Now we need to make bx by arrays
    bx=zeros(nr,nz,nphi,nextcur);
    by=bx;
    for i=1:nr
        for j=1:nz
            for k=1:nphi
                bx(i,j,k,:)=br(i,j,k,:).*cos(data.phi(k))-bphi(i,j,k,:).*sin(data.phi(k));
                by(i,j,k,:)=br(i,j,k,:).*sin(data.phi(k))+bphi(i,j,k,:).*cos(data.phi(k));
            end
        end
    end
    data.bx=bx;
    data.by=by;
    data.curlabel=curlabel;
end
end

