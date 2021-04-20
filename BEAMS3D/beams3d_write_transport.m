function beams3d_write_transport(beam_data_in)
%BEAMS3D_WRITE_TRANSPORT Generate a transport file from BEAMS3D data
%   The BEAMS3D_WRITE_TRANSPORT routine is used for outputting quantities
%   for transport analysis.  The code takes either a BEAMS3D data structure
%   as returned by READ_BEAMS3D or a cell array of those data structures.
%   When a cell array is passed it is assumed you want to output the mean
%   values of those runs with error bar generated from the range of values
%   present in those simulations. If a single data structure is passed no
%   error bars are generated.  A text file with transport quantities is
%   written as is an overview plot.
%
% Example usage
%      beam_data{1}=read_beams3d('beams3d_run1.h5');
%      beam_data{2}=read_beams3d('beams3d_run1.h5');
%      beams3d_write_transport(beam_data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

beamdex=[];
lplot=1;
filename='beams3d_transport';

% Act like it's a cell array even if it's not
if ~iscell(beam_data_in)
    beam_data{1} = beam_data_in;
else
    beam_data = beam_data_in;
end

% Find number of runs
nruns = length(beam_data);

% Figure out rho size
nrho=[]; nbeams=[];
for i=1:nruns
    %rho = beam_data{i}.rho;
    nrho = [nrho beam_data{i}.ns_prof1];
    nbeams = [nbeams beam_data{i}.nbeams];
end
%rho_min = mean(rho_min);
%rho_max = mean(rho_max);
nrho = mean(nrho);
drho = 1.0/double(nrho);
rho_out = 0.0:drho:1.0;
rhol_out = rho_out(1:nrho-1);
rhor_out = rho_out(2:nrho);
rho_out = 0.5.*(rhol_out+rhor_out);
nrho = nrho - 1;

% Check nbeams
if (length(unique(nbeams))>1)
    disp('  One of the supplied runs has a different number or beams');
    return;
end

% Default to all beamlines
if isempty(beamdex)
    beamdex = 1:nbeams(1);
end

% Define output arrays
output=[];
for type={'ndot_prof','epower_prof','ipower_prof','j_prof'}
    output.(type{1})=zeros(nruns,nrho);
    for i=1:nruns
        rho = beam_data{i}.rho;
        val = beam_data{i}.(type{1});
        val = sum(val(beamdex,:),1); % sum over nbeams
        output.(type{1})(i,:) = pchip(rho,val,rho_out);
    end
end

% We need the deposition profile
output.depo=zeros(nruns,nrho);
output.vp=zeros(nruns,nrho);
for i=1:nruns
    d=1;
    if (beam_data{i}.ldepo), d=3; end
    dex = false(1,beam_data{i}.nparticles);
    for j=beamdex
        dex = or(dex,beam_data{i}.Beam'==j);
    end
    for j=1:nrho
        l = and(beam_data{i}.S_lines(d,dex)>=rhol_out(j),...
            beam_data{i}.S_lines(d,dex)<rhor_out(j));
        output.depo(i,j) = sum(beam_data{i}.Weight(l));
    end
    [s,~,Vp] = beams3d_volume(beam_data{i});
    rho = sqrt(s);
    Vp = Vp.*2.*rho;
    output.Vp(i,:) = pchip(rho,Vp,rho_out);
end
output.depo = output.depo./(output.Vp.*drho);

% Create output array
out_arr=[rho_out; mean(output.Vp,1); ...
    mean(output.depo,1); range(output.depo,1);...
    mean(output.ndot_prof,1); range(output.ndot_prof,1);...
    mean(output.epower_prof,1); range(output.epower_prof,1);...
    mean(output.ipower_prof,1); range(output.ipower_prof,1);...
    mean(output.j_prof,1); range(output.j_prof,1)];

% Output to a file
fid=fopen([filename '.txt'],'w');
fprintf(fid,['#           rho   dV/drho [m^-3]  ndote [part/m^3/s]'...
    '  err_ndote ndoti [part/m^3/s]    err_ndoti    P_e [W/m^3]'...
    '          err_P_e     P_i [W/m^3]         err_P_i'...
    '       j [A/m^2]           err_j\n']);
fprintf(fid,'%15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n',out_arr);
fclose(fid);


% Part where we do a plot
if lplot
    fig = figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    subplot(2,2,1);
    plot(out_arr(1,:),out_arr(3,:)./1E19,'b','LineWidth',2); hold on;
    plot(out_arr(1,:),out_arr(5,:)./1E19,'r','LineWidth',2);
    xlim([0 1]);
    if nruns>1
        fill([out_arr(1,:) fliplr(out_arr(1,:))],([-out_arr(4,:)./2+out_arr(3,:) fliplr(out_arr(4,:)./2+out_arr(3,:))])./1E19,'b','FaceAlpha',0.5);
        fill([out_arr(1,:) fliplr(out_arr(1,:))],([-out_arr(6,:)./2+out_arr(5,:) fliplr(out_arr(6,:)./2+out_arr(5,:))])./1E19,'r','FaceAlpha',0.5);
        text(0.1,max(ylim)*0.2,['Elec: ' num2str(sum(out_arr(3,:).*out_arr(2,:).*drho),'%9.2e') '\pm' num2str(sum(out_arr(4,:).*out_arr(2,:).*drho),'%9.2e') ' [part/s]'],'FontSize',14);
        text(0.1,max(ylim)*0.1,['Ion:  ' num2str(sum(out_arr(5,:).*out_arr(2,:).*drho),'%9.2e') '\pm' num2str(sum(out_arr(6,:).*out_arr(2,:).*drho),'%9.2e') ' [part/s]'],'FontSize',14);
    else
        text(0.1,max(ylim)*0.2,['Elec: ' num2str(sum(out_arr(3,:).*out_arr(2,:).*drho),'%9.2e')  ' [part/s]'],'FontSize',14);
        text(0.1,max(ylim)*0.1,['Ion:  ' num2str(sum(out_arr(5,:).*out_arr(2,:).*drho),'%9.2e')  ' [part/s]'],'FontSize',14);
    end
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('\Gamma_k [part/(m^3s)]'); legend('Electron (birth)','Ion (therm)');
    title('BEAMS3D Particle Source');
    subplot(2,2,2);
    plot(out_arr(1,:),out_arr(11,:)./1E3,'k','LineWidth',2); hold on;
    xlim([0 1]);
    if nruns>1
        fill([out_arr(1,:) fliplr(out_arr(1,:))],([-out_arr(12,:)./2+out_arr(11,:) fliplr(out_arr(12,:)./2+out_arr(11,:))])./1E3,'k','FaceAlpha',0.5);
        text(0.1,max(ylim)*0.1+min(ylim),['I:  ' num2str(sum(out_arr(11,:).*out_arr(2,:).*drho)./1E3,'%9.2e') '\pm' num2str(sum(out_arr(12,:).*out_arr(2,:).*drho)./1E3,'%9.2e') ' [kA]'],'FontSize',14);
    else
        text(0.1,max(ylim)*0.1+min(ylim),['I: ' num2str(sum(out_arr(11,:).*out_arr(2,:).*drho)./1E3,'%9.2e')  ' [kA]'],'FontSize',14);
    end
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('j_{fast} [kA/m^2]');
    title('Fast Ion Current');
    subplot(2,2,3);
    plot(out_arr(1,:),out_arr(7,:)./1E3,'b','LineWidth',2); hold on;
    xlim([0 1]);
    if nruns>1
        fill([out_arr(1,:) fliplr(out_arr(1,:))],([-out_arr(8,:)./2+out_arr(7,:) fliplr(out_arr(8,:)./2+out_arr(7,:))])./1E3,'b','FaceAlpha',0.5);
        text(0.1,max(ylim)*0.1+min(ylim),['Q_e:  ' num2str(sum(out_arr(7,:).*out_arr(2,:).*drho)./1E3,'%9.2e') '\pm' num2str(sum(out_arr(8,:).*out_arr(2,:).*drho)./1E3,'%9.2e') ' [kW]'],'FontSize',14);
    else
        text(0.1,max(ylim)*0.1+min(ylim),['Q_e: ' num2str(sum(out_arr(7,:).*out_arr(2,:).*drho)./1E3,'%9.2e')  ' [kW]'],'FontSize',14);
    end
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Heating [kW/m^3]');
    title('Electron Heating');
    subplot(2,2,4);
    plot(out_arr(1,:),out_arr(9,:)./1E3,'r','LineWidth',2); hold on;
    xlim([0 1]);
    if nruns>1
        fill([out_arr(1,:) fliplr(out_arr(1,:))],([-out_arr(10,:)./2+out_arr(9,:) fliplr(out_arr(10,:)./2+out_arr(9,:))])./1E3,'r','FaceAlpha',0.5);
        text(0.1,max(ylim)*0.1+min(ylim),['Q_i:  ' num2str(sum(out_arr(9,:).*out_arr(2,:).*drho)./1E3,'%9.2e') '\pm' num2str(sum(out_arr(10,:).*out_arr(2,:).*drho)./1E3,'%9.2e') ' [kW]'],'FontSize',14);
    else
        text(0.1,max(ylim)*0.1+min(ylim),['Q_i: ' num2str(sum(out_arr(9,:).*out_arr(2,:).*drho)./1E3,'%9.2e')  ' [kW]'],'FontSize',14);
    end
    set(gca,'FontSize',24);
    xlabel('r/a'); ylabel('Heating [kW/m^3]');
    title('Ion Heating');
    saveas(fig,[filename '.fig']);
    saveas(fig,[filename '.png']);
end

end

