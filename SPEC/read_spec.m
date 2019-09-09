function data = read_spec(filename)
%READ_SPEC Reads the HDF5 file produced by the SPEC code.
%   The READ_SPEC function reads the HDF5 file produced by the SPEC code
%   and returns it's contents as elements of a structure.
%
% Example usage
%      data=read_spec('spec_test.h5');  % Reads SPEC HDF5 file
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.0


% Read the HDF5 File
try
    data = read_hdf5(filename);
catch file_err
    disp(['Could not open file: ' filename]);
    data=-1;
    return
end

% Now adjust the variables
mu0=4*pi*1.0E-7;
data.mnmax    = data.mn;
data.phi_step = data.tflux;
data.iota_step = data.iota;
data.xm = double(data.im)';
data.xn = double(data.in)';
data.nu = 90;
data.mpol = max(data.xm);
data.nfp = round(min(data.xn(data.xn>0)));
if isempty(data.nfp), data.nfp = 1; end
if (max(data.xn) > 0)
    data.ntor = round(max(data.xn)/data.nfp);
else
    data.ntor = 0;
end
dex = max(strfind(filename,'.h5'))-1;
data.input_extension = filename(1:dex);
data.datatype='SPEC';

% Now create the grid
nsub=10.0;
%if nsub*data.Nvol < 100
%    nsub = round(100/data.Nvol);
%end
nvol = double(data.Nvol);
data.rmnc=data.Rbc(:,1);
data.zmns=data.Zbs(:,1);
data.rmns=data.Rbs(:,1);
data.zmnc=data.Zbc(:,1);
data.pressure=data.pressure./mu0;
data.presf=data.pressure(1);
data.muf=data.mu(1);
data.phif=0.0;
temp_flux = [0; data.tflux];
dx = 1./(nsub);
x=dx:dx:1;
nx=length(x);
data.rmncp(:,1) = 0.5*(data.Rbc(:,2)-data.Rbc(:,1));
data.zmnsp(:,1) = 0.5*(data.Zbs(:,2)-data.Zbs(:,1));
data.rmnsp(:,1) = 0.5*(data.Rbs(:,2)-data.Rbs(:,1));
data.zmncp(:,1) = 0.5*(data.Zbc(:,2)-data.Zbc(:,1));
for i=1:nvol
    if (i == 1)
        for mn=1:data.mnmax
            mval = 0.5*data.xm(mn);
            data.rmnc(mn,2:nx+1)  = data.Rbc(mn,2).*x.^mval + data.Rbc(mn,1).*(1-x.^mval);
            data.rmns(mn,2:nx+1)  = data.Rbs(mn,2).*x.^mval + data.Rbs(mn,1).*(1-x.^mval);
            data.zmnc(mn,2:nx+1)  = data.Zbc(mn,2).*x.^mval + data.Zbc(mn,1).*(1-x.^mval);
            data.zmns(mn,2:nx+1)  = data.Zbs(mn,2).*x.^mval + data.Zbs(mn,1).*(1-x.^mval);
            data.rmncp(mn,2:nx+1) = 0.5.*mval*x.^(mval-1).*data.Rbc(mn,2)-0.5.*mval*x.^(mval-1).*data.Rbc(mn,1);
            data.rmnsp(mn,2:nx+1) = 0.5.*mval*x.^(mval-1).*data.Rbs(mn,2)-0.5.*mval*x.^(mval-1).*data.Rbs(mn,1);
            data.zmncp(mn,2:nx+1) = 0.5.*mval*x.^(mval-1).*data.Zbc(mn,2)-0.5.*mval*x.^(mval-1).*data.Zbc(mn,1);
            data.zmnsp(mn,2:nx+1) = 0.5.*mval*x.^(mval-1).*data.Zbs(mn,2)-0.5.*mval*x.^(mval-1).*data.Zbs(mn,1);
        end
        data.rmncp(:,1)      = 2.0*data.rmncp(:,2)-data.rmncp(:,3);
        data.rmnsp(:,1)      = 2.0*data.rmnsp(:,2)-data.rmnsp(:,3);
        data.zmncp(:,1)      = 2.0*data.zmncp(:,2)-data.zmncp(:,3);
        data.zmnsp(:,1)      = 2.0*data.zmnsp(:,2)-data.zmnsp(:,3);
        
    else
        data.rmnc  = [data.rmnc  interp1([0 1],data.Rbc(:,i:i+1)',x,'linear')'];
        data.rmncp = [data.rmncp 0.5.*repmat((data.Rbc(:,i+1)-data.Rbc(:,i)),[1 nx])];
        data.zmns  = [data.zmns  interp1([0 1],data.Zbs(:,i:i+1)',x,'linear')'];
        data.zmnsp = [data.zmnsp 0.5.*repmat((data.Zbs(:,i+1)-data.Zbs(:,i)),[1 nx])];
        data.rmns  = [data.rmns  interp1([0 1],data.Rbs(:,i:i+1)',x,'linear')'];
        data.rmnsp = [data.rmnsp 0.5.*repmat((data.Rbs(:,i+1)-data.Rbs(:,i)),[1 nx])];
        data.zmnc  = [data.zmnc  interp1([0 1],data.Zbc(:,i:i+1)',x,'linear')'];
        data.zmncp = [data.zmncp 0.5.*repmat((data.Zbc(:,i+1)-data.Zbc(:,i)),[1 nx])];
    end
    data.presf= [data.presf zeros(1,length(dx:dx:1))+data.pressure(i)];
    data.muf= [data.muf zeros(1,length(dx:dx:1))+data.mu(i)];
    data.phif = [data.phif; interp1([0 1],temp_flux(i:i+1),x,'linear')'];
end
data.ideal_arr=nsub+1:nsub:size(data.rmnc,2);
data.phi = data.phif;
data.nsub=nsub;
data.ns  =size(data.rmnc,2);
data.pres_step=data.presf;
% Make the offset plots for stepped equilibria
offset = 1;
data.pres_plot= data.presf;
data.phi_plot = data.phif';
data.mu_plot  = data.muf;
for i=1:data.Nvol-1
    dex = i*data.nsub+offset;
    data.pres_plot = [data.pres_plot(1:dex) data.pres_plot(dex) data.pres_plot(dex+1:end)];
    data.mu_plot = [data.mu_plot(1:dex) data.mu_plot(dex) data.mu_plot(dex+1:end)];
    data.phi_plot = [data.phi_plot(1:dex) data.phi_plot(dex+1) data.phi_plot(dex+1:end)];
    offset = offset + 1;
end
%%% After this point read the 'other' files

% Read the poincare files
machine_format='a';
int_format='int32';
float_format='float64';
spacer_format='int32';
j=1;
for i=1:nvol
    offset = real(i-1)./real(nvol);
    try
        poincare_file=['.' filename(1:length(filename)-3) '.poincare.' num2str(i,'%4.4i')];
        fid = fopen(poincare_file,'r',machine_format);
        if (fid > 0)
            while ~feof(fid)
                % Interface Labels
                fread(fid,1,spacer_format);
                if (feof(fid)), break; end;
                temp=fread(fid,2,int_format);
                fread(fid,1,spacer_format);
                nzeta=temp(1); %NZ
                nppts=temp(2); %NPPTS
                data.npoinc=nppts;
                % blockdata
                fread(fid,1,spacer_format);
                temp=fread(fid,4*nzeta*nppts,float_format);
                fread(fid,1,spacer_format);
                if ~isempty(temp)
                    temp=reshape(temp,[4 nzeta nppts]);
                    data.th_lines(j,1:nzeta,1:nppts)=temp(1,1:nzeta,1:nppts);
                    rho = 0.5*(temp(2,1:nzeta,1:nppts)+1)./real(nvol);
                    data.rho_lines(j,1:nzeta,1:nppts)=rho+offset;
                    data.R_lines(j,1:nzeta,1:nppts)=temp(3,1:nzeta,1:nppts);
                    data.Z_lines(j,1:nzeta,1:nppts)=temp(4,1:nzeta,1:nppts);
                    j=j+1;
                end
            end
            fclose(fid);
        end
    catch
        disp(' - Could not read poincare file');
    end
end

% Read the magnetic field files
machine_format='a';
int_format='int32';
float_format='float64';
spacer_format='int32';
vecpot_file=['.' filename(1:length(filename)-3) '.AtAzmn'];
data.Asubumnc = zeros(size(data.rmnc));
data.Asubumns = zeros(size(data.rmnc));
data.Asubvmnc = zeros(size(data.rmnc));
data.Asubvmns = zeros(size(data.rmnc));
try
    fid = fopen(vecpot_file,'r',machine_format);
    fread(fid,1,spacer_format);
    temp=fread(fid,2,int_format);
    fread(fid,1,spacer_format);
    Nvol = temp(1);
    mnmax_b   = temp(2);
    fread(fid,1,spacer_format);
    im=fread(fid,mnmax_b,int_format);
    fread(fid,1,spacer_format);
    fread(fid,1,spacer_format);
    in=fread(fid,mnmax_b,int_format);
    fread(fid,1,spacer_format);
    for i = 1:Nvol
        fread(fid,1,spacer_format);
        lrad=fread(fid,1,int_format);
        fread(fid,1,spacer_format);
        for j=1:mnmax_b
            fread(fid,1,spacer_format);
            data.Ate(i,j,1:lrad+1)=fread(fid,lrad+1,float_format);
            fread(fid,1,spacer_format);
            fread(fid,1,spacer_format);
            data.Aze(i,j,1:lrad+1)=fread(fid,lrad+1,float_format);
            fread(fid,1,spacer_format);
            fread(fid,1,spacer_format);
            data.Ato(i,j,1:lrad+1)=fread(fid,lrad+1,float_format);
            fread(fid,1,spacer_format);
            fread(fid,1,spacer_format);
            data.Azo(i,j,1:lrad+1)=fread(fid,lrad+1,float_format);
            fread(fid,1,spacer_format);
        end
    end
    fclose(fid);
    % Now get into AsubXmnL where X is u/v and L is c/s
    data.Asubumncp = zeros(size(data.rmnc));
    data.Asubumnsp = zeros(size(data.rmnc));
    data.Asubvmncp = zeros(size(data.rmnc));
    data.Asubvmnsp = zeros(size(data.rmnc));
    xmreg=0.5*min(data.xm,4);
    for mn=1:mnmax_b
        for i=1:data.ns
            vol=sum(data.ideal_arr < i)+1;
            x=2.0*(i-data.ideal_arr(vol)-1)./(data.nsub+1)+1;
            if (vol == 1)
                x_bar = 0.5*(x+1.0);
                factor  = x_bar.^xmreg(mn);
                factorp = 0.5*xmreg(mn).*x_bar.^(xmreg(mn)-1);
                for j=1:size(data.Ate,3)
                    ch_poly              = ChebyshevPoly(j-1).*data.Ate(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubumnc(mn,i)  = data.Asubumnc(mn,i) + polyval(ch_poly,x).*factor;
                    data.Asubumncp(mn,i) = data.Asubumncp(mn,i) + polyval(chder_poly,x).*factor + polyval(ch_poly,x).*factorp;
                    ch_poly              = ChebyshevPoly(j-1).*data.Ato(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubumns(mn,i)  = data.Asubumns(mn,i) + polyval(ch_poly,x).*factor;
                    data.Asubumnsp(mn,i) = data.Asubumnsp(mn,i) + polyval(chder_poly,x).*factor + polyval(ch_poly,x).*factorp;
                    ch_poly              = ChebyshevPoly(j-1).*data.Aze(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubvmnc(mn,i)  = data.Asubvmnc(mn,i) + polyval(ch_poly,x).*factor;
                    data.Asubvmncp(mn,i) = data.Asubvmncp(mn,i) + polyval(chder_poly,x).*factor + polyval(ch_poly,x).*factorp;
                    ch_poly              = ChebyshevPoly(j-1).*data.Azo(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubvmns(mn,i)  = data.Asubvmns(mn,i) + polyval(ch_poly,x).*factor;
                    data.Asubvmnsp(mn,i) = data.Asubvmnsp(mn,i) + polyval(chder_poly,x).*factor + polyval(ch_poly,x).*factorp;
                    data.xvec(vol,i)=x;
                end
            else
                for j=1:size(data.Ate,3)
                    ch_poly              = ChebyshevPoly(j-1).*data.Ate(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubumnc(mn,i)  = data.Asubumnc(mn,i) + polyval(ch_poly,x);
                    data.Asubumncp(mn,i) = data.Asubumncp(mn,i) + polyval(chder_poly,x);
                    ch_poly              = ChebyshevPoly(j-1).*data.Ato(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubumns(mn,i)  = data.Asubumns(mn,i) + polyval(ch_poly,x);
                    data.Asubumnsp(mn,i) = data.Asubumnsp(mn,i) + polyval(chder_poly,x);
                    ch_poly              = ChebyshevPoly(j-1).*data.Aze(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubvmnc(mn,i)  = data.Asubvmnc(mn,i) + polyval(ch_poly,x);
                    data.Asubvmncp(mn,i) = data.Asubvmncp(mn,i) + polyval(chder_poly,x);
                    ch_poly              = ChebyshevPoly(j-1).*data.Azo(vol,mn,j);
                    chder_poly           = polyder(ch_poly);
                    data.Asubvmns(mn,i)  = data.Asubvmns(mn,i) + polyval(ch_poly,x);
                    data.Asubvmnsp(mn,i) = data.Asubvmnsp(mn,i) + polyval(chder_poly,x);
                    data.xvec(vol,i)=x;
                end
            end
        end
    end
    data.Asubumncp(:,1) = 2.0.*data.Asubumncp(:,2)-data.Asubumncp(:,3);
    data.Asubumnsp(:,1) = 2.0.*data.Asubumnsp(:,2)-data.Asubumnsp(:,3);
    data.Asubvmncp(:,1) = 2.0.*data.Asubvmncp(:,2)-data.Asubvmncp(:,3);
    data.Asubvmnsp(:,1) = 2.0.*data.Asubvmnsp(:,2)-data.Asubvmnsp(:,3);
catch
    disp(' - Could not read Vector Potential file');
end

% Read the transform files
machine_format='a';
int_format='int32';
float_format='float64';
spacer_format='int32';
data.iotaf=[];
dx = 1./(nsub);
x=dx:dx:1;
for i=1:nvol
    try
        transform_file=['.' filename(1:length(filename)-3) '.transform.' num2str(i,'%4.4i')];
        fid = fopen(transform_file,'r',machine_format);
        if (fid > 0)
            % Get subgrid
            fread(fid,1,spacer_format);
            npts=fread(fid,1,int_format);
            fread(fid,1,spacer_format);
            % Get Iota edges
            fread(fid,1,spacer_format);
            iotae=fread(fid,2,float_format);
            fread(fid,1,spacer_format);
            % Get Profile
            fread(fid,1,spacer_format);
            iota_temp=fread(fid,[2 npts],float_format);
            fread(fid,1,spacer_format);
            % Now make the grid
            if i == 1, data.iotaf = iota_temp(2,1); end
            data.iotaf = [data.iotaf pchip(0:1/(npts-1):1.,iota_temp(2,1:npts),x)];
        else
            data.iotaf = [iotaf zeros(1,nsub)];
        end
        fclose(fid);
    catch
        disp(' - Could not read transform file');
    end
end

    
return
end
% Note on Chebyshev polynomials
%   To(x) = 1
%   T1(x) = x
%   T(n+1) = 2xTn(x) - Tn-1(x)

