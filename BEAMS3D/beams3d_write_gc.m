function groupid = beams3d_write_gc(beam_data,dex,varargin)
%BEAMS3D_WRITE_GC Outputs a subset of particles as input
%   The BEAMS3D_WRITE_GC function outputs to a file a subset of particles
%   based on a beams3d run.  It takes a beams3d data strucuture as returned
%   by read_beams3d and an dex array as inputs.  The dex array should be
%   nparticles long and contain the index in terms of npoinc from which to
%   pull the particle (setting a value to <=0 will skip that particle).  If
%   a empty array is supplied to dex then the last datapoint will be used
%   for each particle.
%
%   If the optional argument 'ascot5' is provided then an HDF5 file with
%   the /marker/gc_ group will be produced.
%
%   If the optional argument (...,'ndiv',ndiv) is provided then ndiv will
%   be used to divide up each selected marker into submarkers.  These
%   submarkers have the same initial conditions as the orriginal marker but
%   a reduced weight.
%
%   Example
%       beam_data=read_beams3d('beams3d_test.h5');
%       beams3d_write_gc(beam_data,[]); % last points
%
%       dex = [1 1 0 1 1]; % skip 3rd particle
%       beams3d_write_gc(beam_data,dex);
%
%       beams3d_write_gc(beam_data,dex,'ascot5'); %ascot marker produced
%
%       beamd3d_write_gc(beam_data,dex,'ascot5','ndiv',25);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

% Defaults
lplot=0;
lbeams3d = 1;
lascot5  = 0;
ndiv = 1;
pertcent = 0;
pertdeg=0;
ec = 1.60217662E-19;
amu = 1.66053906660E-27;
aux_str='';
filename = 'beams3d_gc.txt';

% Handle varargin
if ~isempty(varargin)
    n=1;
    while n <= length(varargin)
        switch varargin{n}
            case 'ascot5'
                lbeams3d=0;
                lascot5=1;
            case 'pert'
                n=n+1;
                pertcent = varargin{n};
                aux_str=[aux_str ' pert: ',num2str(pertcent,'%5.4f') ';'];
            case 'pert_deg'
                n=n+1;
                pertdeg = varargin{n};
                aux_str=[aux_str ' pert: ',num2str(pertdeg,'%3.1f') ' ^o;'];
            case 'ndiv'
                n=n+1;
                ndiv = varargin{n};
                aux_str=[aux_str ' ndiv: ',num2str(pertcent,'%5i') ';'];
            case 'filename'
                n=n+1;
                filename = varargin{n};
            case {'plot','plots'}
                lplot=1;
        end
        n=n+1;
    end
end

% Handle empty dex
if isempty(dex)
    mask = ones(1,beam_data.nparticles)==1;
    for i=1:beam_data.nparticles
        dex_local(i) = find(beam_data.R_lines(:,i)>0,1,'last');
    end
else
    if length(dex) ~= beam_data.nparticles
        disp('ERROR: dex array must be beam_data.nparticles long.');
        disp(['     LENGTH(DEX): ' num2str(length(dex),'%i')]);
        disp(['     NPARTICLES:  ' num2str(beam_data.nparticles,'%i')]);
        return;
    end
    if any(dex > beam_data.npoinc+1)
        disp('ERROR: dex values must not be larger than beam_data.npoinc+1');
        disp(['     MAX(DEX): ' num2str(max(dex),'%i')]);
        disp(['     NPOINC+1:  ' num2str(beam_data.npoinc+1,'%i')]);
        return;
    end
    mask = dex>0;
    dex_local = dex;
end

U_lines=beams3d_fixUlines(beam_data);

% Now process subarrays
nnew=sum(dex_local>0);
nold   = nnew;
rho    = zeros(1,nnew); %temp
u      = zeros(1,nnew);
r      = zeros(1,nnew);
phi    = zeros(1,nnew);
z      = zeros(1,nnew);
mu     = zeros(1,nnew);
vll    = zeros(1,nnew);
b      = zeros(1,nnew);
mass   = beam_data.mass(mask)';
charge = beam_data.charge(mask)';
zatom  = beam_data.Zatom(mask)';
tend  = beam_data.t_end(mask)';
weight  = beam_data.Weight(mask)';
k=1;
for i=1:beam_data.nparticles
    if mask(i)
        j = dex_local(i);
        rho(k) = sqrt(beam_data.S_lines(j,i));
        u(k) = U_lines(j,i);
        r(k)   = beam_data.R_lines(j,i);
        phi(k) = beam_data.PHI_lines(j,i);
        z(k)   = beam_data.Z_lines(j,i);
        mu(k)  = beam_data.moment_lines(j,i);
        vll(k) = beam_data.vll_lines(j,i);
        b(k)   = beam_data.B_lines(j,i);
        k=k+1;
    end
end

% Now handle ndiv
if ndiv>1
    nnew   = ndiv*nnew;
    r      = reshape(repmat(r,[ndiv,1]),[1 nnew]);
    z      = reshape(repmat(z,[ndiv,1]),[1 nnew]);
    phi    = reshape(repmat(phi,[ndiv,1]),[1 nnew]);
    mu     = reshape(repmat(mu,[ndiv,1]),[1 nnew]);
    vll    = reshape(repmat(vll,[ndiv,1]),[1 nnew]);
    b      = reshape(repmat(b,[ndiv,1]),[1 nnew]);
    mass   = reshape(repmat(mass,[ndiv,1]),[1 nnew]);
    charge = reshape(repmat(charge,[ndiv,1]),[1 nnew]);
    zatom  = reshape(repmat(zatom,[ndiv,1]),[1 nnew]);
    tend   = reshape(repmat(tend,[ndiv,1]),[1 nnew]);
    weight = reshape(repmat(weight,[ndiv,1]),[1 nnew])./ndiv;
end

% Add perturbation
if or(pertcent > 0,pertdeg >0)
    vperp2 = 2.*mu.*b./mass;
    vtotal = sqrt(vll.*vll+vperp2);
    pitch  = vll./vtotal;
    if (max(vtotal) > 1E6)
        scale = 1E-6;
        units = ' x10^3 [km/s]';
    else
        scale = 1E-3;
        units = ' [km/s]';
    end
    fig = figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    plot(vll.*scale,sqrt(vperp2).*scale,'or');
    if pertcent >0
        temp   = 1 - pertcent + 2*pertcent.*rand(1,length(pitch));
        temp(1:ndiv:end) = 1; % Keep the original particle
        pitch  = pitch.*temp;
    elseif pertdeg >0
        temp   = cosd(acosd(pitch) - pertdeg + 2*pertdeg.*rand(1,length(pitch))); % pert in deg
        temp(1:ndiv:end) = pitch(1:ndiv:end); % Keep the original particle
        pitch  = temp;
    end
    vll    = pitch.*vtotal;
    vperp2 = abs(vtotal.*vtotal - vll.*vll);
    mu     = 0.5.*vperp2.*mass./b;
    hold on; plot(vll.*scale,sqrt(vperp2).*scale,'.b');
    set(gca,'FontSize',24);
    xlabel(['V_{||}' units]);
    ylabel(['V_\perp' units]);
    title('Particle Distribution');
    legend('Original (BEAMS3D)','New (ASCOT5)');
    xlim([-1.05 1.05].*max(abs(vtotal)).*scale);
    ylim([0 1.05].*max(abs(vtotal)).*scale);
    annotation('textbox',[0.15 0.15 0.1 0.1],...
        'String',['Original: ' num2str(nold,'%i') ' particles'],...
        'FontSize',24,'LineStyle','none');
    annotation('textbox',[0.15 0.10 0.1 0.1],...
        'String',['New: ' num2str(nnew,'%i') ' particles'],...
        'FontSize',24,'LineStyle','none');
    if pertcent >0
        annotation('textbox',[0.75 0.10 0.1 0.1],...
            'String',['Pert.: ' num2str(pertcent.*100,'%5.2f') ' %'],...
            'FontSize',24,'LineStyle','none');
    elseif pertdeg >=0
        annotation('textbox',[0.75 0.10 0.1 0.1],...
            'String',['Pert.: ' num2str(pertdeg,'%3.1f') ' ^o'],...
            'FontSize',24,'LineStyle','none');
    end
end

% Now output
if lbeams3d
    fid = fopen(filename,'w');
    fprintf(fid,'&BEAMS3D_INPUT\n');
    fprintf(fid,['! Created by beams3d_write_gc. ' aux_str]);
    fprintf(fid,['!   DATE: ' disp(datestr(now,'mm-dd-yyyy HH:MM:SS'))]);
    fprintf(fid,'  NR = %i\n',beam_data.nr);
    fprintf(fid,'  NZ = %i\n',beam_data.nz);
    fprintf(fid,'  NPHI = %i\n',beam_data.nphi);
    fprintf(fid,'  RMIN = %20.10E\n',beam_data.raxis(1));
    fprintf(fid,'  RMAX = %20.10E\n',beam_data.raxis(end));
    fprintf(fid,'  ZMIN = %20.10E\n',beam_data.zaxis(1));
    fprintf(fid,'  ZMAX = %20.10E\n',beam_data.zaxis(end));
    fprintf(fid,'  PHIMIN = %20.10E\n',beam_data.phiaxis(1));
    fprintf(fid,'  PHIMAX = %20.10E\n',beam_data.phiaxis(end));
    fprintf(fid,'  INT_TYPE = ''LSODE''\n');
    fprintf(fid,'  FOLLOW_TOL = 1.0E-8\n');
    fprintf(fid,'  VC_ADAPT_TOL = 1.0E-3\n');
    fprintf(fid,'  NPOINC = %i\n',beam_data.npoinc);
    if range(charge) == 0
        fprintf(fid,'  CHARGE_IN = %i*%-20.10E\n',nnew,mean(charge));
    else
        fprintf(fid,'  CHARGE_IN = ');
        fprintf(fid,' %20.10E ',charge);
        fprintf(fid,'\n');
    end
    if range(mass) == 0
        fprintf(fid,'  MASS_IN = %i*%-20.10E\n',nnew,mean(mass));
    else
        fprintf(fid,'  MASS_IN = ');
        fprintf(fid,' %20.10E ',mass);
        fprintf(fid,'\n');
    end
    if range(zatom) == 0
        fprintf(fid,'  ZATOM_IN = %i*%-20.10E\n',nnew,mean(zatom));
    else
        fprintf(fid,'  ZATOM_IN = ');
        fprintf(fid,' %20.10E ',zatom);
        fprintf(fid,'\n');
    end
    if range(tend) == 0
        fprintf(fid,'  T_END_IN = %i*%-20.10E\n',nnew,mean(tend));
    else
        fprintf(fid,'  T_END_IN = ');
        fprintf(fid,' %20.10E ',tend);
        fprintf(fid,'\n');
    end
    fprintf(fid,'  R_START_IN = ');
    fprintf(fid,' %20.10E ',r);
    fprintf(fid,'\n');
    fprintf(fid,'  Z_START_IN = ');
    fprintf(fid,' %20.10E ',z);
    fprintf(fid,'\n');
    fprintf(fid,'  PHI_START_IN = ');
    fprintf(fid,' %20.10E ',phi);
    fprintf(fid,'\n');
    fprintf(fid,'  VLL_START_IN = ');
    fprintf(fid,' %20.10E ',vll);
    fprintf(fid,'\n');
    fprintf(fid,'  MU_START_IN = ');
    fprintf(fid,' %20.10E ',mu);
    fprintf(fid,'\n');
    fprintf(fid,'/\n');
    fclose(fid);
elseif lascot5
    filename=strrep(filename,'beams3d','ascot5');
    filename=strrep(filename,'.txt','.h5');
    disp('WRITING ASCOT5 HDF GC File!');
    phi2   = rad2deg(phi);
    vperp2 = 2.*mu.*b./mass;
    vtotal = sqrt(vll.*vll+vperp2);
    energy = 0.5.*mass.*vtotal.*vtotal./ec;
    pitch  = vll./vtotal;
    %zeta   = zeros(1,nnew);
    zeta   = rand(1,nnew).*pi.*2;
    charge2 = round(charge./ec);
    anum    = round(mass./amu);
    znum    = round(charge./ec);
    mass2   = mass./amu;
    id      = 1:nnew;
    ascot5_writemarker_gc(filename,r,phi2,z,energy,pitch,...
        zeta,mass2,charge2,anum,znum,weight,tend,id);
    h5writeatt(filename,['/marker'],'notes',['Created in MATLAB via beams3d_write_gc. ' aux_str],'TextEncoding','system');
    if lplot
        figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
        subplot(2,2,1);
        plot3(r.*cosd(phi2),r.*sind(phi2),z,'.k');
        axis off;
        subplot(2,2,2);
        plot(rho.*cos(u),rho.*sin(u),'.k');
        hold on;
        th=0:360;
        plot(1.5.*cosd(th),1.25.*sind(th),'--k');
        plot(1.0.*cosd(th),1.0.*sind(th),'--k');
        plot(0.5.*cosd(th),0.5.*sind(th),'--k');
        set(gca,'FontSize',24);
        axis equal;
        title('\rho vs. \theta^*');
        subplot(2,2,4);
        plot(pitch,energy./1E3,'.k');
        set(gca,'FontSize',24);
        xlabel('V_{ll}/V_{total}');
        ylabel(' Energy [keV]');
    end

end

    

return;
end

