function beams3d_transport_plots(ext)
%BEAMS3D_TRASPORT_PLOTS Outputs transport data to a file.
%   This subroutine reads beams3d_ext_birth.h5 and beams3d_ext_slow.h5
%   files to create an overview plot and an overview text file. It assumes
%   the _birth file has a deposition BEAMS3D run and _slow is a resulting
%   slowing down run.
%
% Example usage
%      beams3d_transport_plots('test_run');  
%
% See also:
%   BEAMS3D_WRITE_TRANSPORT
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

file_birth = ['beams3d_' ext '_birth.h5'];
file_slow = ['beams3d_' ext '_slow.h5'];
out_file = ['overview_' ext '.txt'];

% Defaults
lplot=1;
nout = 64;

try
    input_file=['input.' ext '_birth'];
    nml=read_namelist(input_file,'BEAMS3D_INPUT');
    source_dex=nml.dex_beams;
catch
    disp('Warning Could not read inptu file for BEAM indicies');
    source_dex=[];
end

% Create output rho array
hs_out = 1./nout;
rho_out = 0:hs_out:1;
rho_out = 0.5.*(rho_out(1:end-1)+rho_out(2:end));
s_out = rho_out.^2;

% Make Birth Data
beam_data=read_beams3d(file_birth);
% Extract Profile Information
[s, ne, te, ti, zeff] = beams3d_profiles(beam_data);
[s_er, pot, dpotds] = beams3d_er(beam_data);
[s_vol, ~, dVds] = beams3d_volume(beam_data);
dVdr = 2.*sqrt(s_vol).*dVds;
Er = dpotds.*2.*sqrt(s_er); %missing Aminor

% Define number of Sources
epersource = round(length(unique(beam_data.Energy)));
nsources = round(double(beam_data.nbeams)./epersource);
if isempty(source_dex)
    source_dex=1:nsources;
else
    source_dex=source_dex(1:epersource:end);
end

% Birth Data
power = beams3d_powertable(beam_data,'quiet');
birth = beams3d_calc_depo(beam_data);
rho_birth = 0:1./(beam_data.ns_prof):1;
rho_birth = 0.5.*(rho_birth(1:end-1)+rho_birth(2:end));
power_out = zeros(8,nsources); birth_out=zeros(nout,nsources);
for i=1:nsources
    i1=(i-1)*epersource+1;
    i2=i1+epersource-1;
    power_out(:,i)=sum(power(:,i1:i2),2);
    birth_out(:,i)=pchip(rho_birth,sum(birth(i1:i2,:),1),rho_out);
end

% Output profiles
dataout = []; header_out=[];
dataout = [dataout; rho_out]; header_out=[header_out 'r/a; '];
dataout = [dataout; pchip(s,ne,s_out)]; header_out=[header_out 'Ne [m^-3]; '];
dataout = [dataout; pchip(s,te,s_out)]; header_out=[header_out 'Te [eV]; '];
dataout = [dataout; pchip(s,ti,s_out)]; header_out=[header_out 'Ti [eV]; '];
dataout = [dataout; pchip(s,zeff,s_out)]; header_out=[header_out ' Zeff ; '];
dataout = [dataout; pchip(s_er,pot,s_out)]; header_out=[header_out ' Ves [V]; '];
dataout = [dataout; pchip(s_er,Er,s_out)]; header_out=[header_out ' Er*A [V]; '];
dataout = [dataout; pchip(s_vol,dVdr,s_out)]; header_out=[header_out ' dV/dr [m^3]; '];
for i = 1:nsources
    dataout = [dataout; pchip(rho_birth,birth_out(:,i),rho_out)]; header_out=[header_out ' Birth (S' num2str(source_dex(i),'%1.1i') ') [m^-3/s]; '];
end

% Now work on slowing down
beam_data=read_beams3d(file_slow);

% Create thermalization array
therm_index=beams3d_finddex(beam_data,'therm_last');
therm=zeros(beam_data.nbeams,nout);
for i=1:beam_data.nparticles
    if therm_index(i)==0, continue; end
    rhot = sqrt(beam_data.S_lines(therm_index(i),i));
    i1 = max(round(rhot./hs_out),1);
    i2 = beam_data.Beam(i);
    therm(i2,i1)=therm(i2,i1)+beam_data.Weight(i);
end
power = beams3d_powertable(beam_data,'quiet');
dense_out=[]; qe_out=[]; qi_out=[]; therm_out=[]; j_out=[];
for i=1:nsources
    i1=(i-1)*epersource+1;
    i2=i1+epersource-1;
    %power_out(:,i)=sum(power(:,i1:i2),3);
    dense_out(:,i)=sum(beam_data.dense_prof(i1:i2,:),1);
    qe_out(:,i)=sum(beam_data.epower_prof(i1:i2,:),1);
    qi_out(:,i)=sum(beam_data.ipower_prof(i1:i2,:),1);
    therm_out(:,i)=sum(therm(i1:i2,:),1);
    j_out(:,i) = sum(beam_data.j_prof(i1:i2,:),1);
    power_out(4:8,i)=sum(power(4:8,i1:i2),2);
end
for i = 1:nsources
    dataout = [dataout; pchip(rho_birth,dense_out(:,i),rho_out)]; header_out=[header_out ' FI Density (S' num2str(source_dex(i),'%1.1i') ') [m^-3]; '];
end
for i = 1:nsources
    dataout = [dataout; pchip(rho_birth,qe_out(:,i),rho_out)]; header_out=[header_out ' Electron Heating (S' num2str(source_dex(i),'%1.1i') ') [W/m^3]; '];
end
for i = 1:nsources
    dataout = [dataout; pchip(rho_birth,qi_out(:,i),rho_out)]; header_out=[header_out ' Ion Heating (S' num2str(source_dex(i),'%1.1i') ') [W/m^3]; '];
end
for i = 1:nsources
    dataout = [dataout; pchip(rho_birth,therm_out(:,i),rho_out)]; header_out=[header_out ' Thermalized Particles (S' num2str(source_dex(i),'%1.1i') ') [m^-3/s]; '];
end
for i = 1:nsources
    dataout = [dataout; pchip(rho_birth,j_out(:,i).*beam_data.GFactor,rho_out)]; header_out=[header_out ' Current Drive (S' num2str(source_dex(i),'%1.1i') ') [A/m^-2]; '];
end

if lplot
    fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    subplot(2,3,1);
    plot(dataout(1,:),dataout(2,:).*1E-19,'b','LineWidth',3); hold on;
    plot(dataout(1,:),1E-19.*dataout(2,:)./dataout(5,:),'r','LineWidth',3);
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Density x10^{19} [m^{-3}]'); legend('Electrons','Ions');
    subplot(2,3,2);
    plot(dataout(1,:),dataout(3,:).*1E-3,'b','LineWidth',3); hold on;
    plot(dataout(1,:),1E-3.*dataout(4,:),'r','LineWidth',3);
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Temperature [keV]');
    subplot(2,3,3);
    plot(dataout(1,:),dataout(7,:).*1E-3,'k','LineWidth',3);
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Er*A_{minor} [kV]');
    subplot(2,3,4);
    j1=9;
    j2=j1+nsources-1;
    plot(dataout(1,:),sum(dataout(j1:j2,:),1).*1E-19,'k','LineWidth',3);
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Birth Rate x10^{19} [m^{-3}/s]');
    subplot(2,3,5);
    % Skip FI density
    j1=j2+1;
    j2=j1+nsources-1;
    j1=j2+1;
    j2=j1+nsources-1;
    plot(dataout(1,:),sum(dataout(j1:j2,:),1).*1E-3,'b','LineWidth',3); hold on;
    j1=j2+1;
    j2=j1+nsources-1;
    plot(dataout(1,:),sum(dataout(j1:j2,:),1).*1E-3,'r','LineWidth',3);
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Heating [kW/m^{-3}]');
    subplot(2,3,6);
    % Skip Thermalized Particles
    j1=j2+1;
    j2=j1+nsources-1;
    j1=j2+1;
    j2=j1+nsources-1;
    plot(dataout(1,:),sum(dataout(j1:j2,:),1).*1E-3,'k','LineWidth',3);
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Current Drive [kA/m^{-2}]');
    saveas(fig,['overview_' ext '.fig']);
    saveas(fig,['overview_' ext '.png']);
    close(fig);
end

fid=fopen(out_file,'w');
fprintf(fid,'BEAMS3D VERSION: %4.2f\n',beam_data.VERSION);
fprintf(fid,['Run Extension: ' ext '\n']);
fprintf(fid,'Plasma Mass: %20.10G\n',beam_data.plasma_mass);
fprintf(fid,'\n');
fprintf(fid,[header_out '\n']);
fmt=''; for i=1:size(dataout,1), fmt=[fmt ' %20.10E']; end
fmt = [fmt '\n'];
fprintf(fid,fmt,dataout);
fprintf(fid,'\n');
fprintf(fid,'    POWER TOTAL   PORTS   SHINE    IONS    ELEC    WALL   ORBIT   THERM\n');
for i=1:size(power_out,2)
    fprintf(fid,'BEAM %1.1i : %6i  %6i  %6i  %6i  %6i  %6i  %6i  %6i  kW\n',i,round(power_out(:,i)./1E3));
end
fclose(fid);

end

