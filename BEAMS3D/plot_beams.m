function [ output_args ] = plot_beams( beam_data,varargin)
%PLOT_BEAMS Makes plots of BEAMS3D data
%The PLOT_BEAMS function creates various canned plot of BEAMS3D data.  The
%particle trajectory data or diagnostic data can be read using the
%READ_BEAMS3D function and the resulting strucutre passed to PLOT_BEAMS.
%
% Example usage
%      beam_data = read_beams3d('beams3d_test.h5');
%      plot_beams(beam_data);
%      plot_beams(beam_data,'overview');
%      plot_beams(beam_data,'xyz');
%      plot_beams(beam_data,'xyz_therm');
%      plot_beams(beam_data,'xyz_lost');
%      plot_beams(beam_data,'xyz_initial');
%      plot_beams(beam_data,'xyz_therm_initial');
%      plot_beams(beam_data,'xyz_lost_initial');
%      plot_beams(beam_data,'xyz_birth');
%      plot_beams(beam_data,'flux_'); % note _opt as in xyz above
%      plot_beams(beam_data,'pitch_'); % note _opt as in xyz above
%      plot_beams(beam_data,'dist_'); % note _opt as in xyz above
%      plot_beams(beam_data,'rho_'); % note _opt as in xyz above
%      plot_beams(beam_data,'orbit_rz'); % Full particle Orbit in RZ
%      plot_beams(beam_data,'orbit_flux'); % Full particle Orbit in flux
%      plot_beams(beam_data,'orbit_vspace'); % Full particle Orbit in V
%      plot_beams(beam_data,'frac_loss'); % Wall Loss fraction
%      plot_beams(beam_data,'frac_therm'); % Thermalized fraction
%      plot_beams(beam_data,'distribution'); % 2D vll/vperp distribution
%      plot_beams(beam_data,'heating'); % 1D Heating profiles
%      plot_beams(beam_data,'current'); % 1D Current Density profiles
%      plot_beams(beam_data,'fueling'); % 1D EP fuelint profiles
%      plot_beams(beam_data,'injection'); % Make an NBI injection plot
%      plot_beams(beam_data,'birth_image'); % R/Z Image plot
%      plot_beams(beam_data,'wall_loss'); % # of particle lost
%      plot_beams(beam_data,'wall_shine'); % [W/m^2]
%      plot_beams(beam_data,'wall_heat'); % [W/m^2]
%      plot_beams(beam_data,'wall_heat_2d'); % [W/m^2]
%      plot_beams(beam_data,'bmir'); % Mirro magnetic field
%      plot_beams(beam_data,'benchmark'); % Old benchmark plots
%      beam_dex = [7:12 19:24];
%      plot_beams(beam_data,'wall_heat','beam',beam_dex); % plot subset.
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.70

ec=1.6021773300E-19; % Charge of an electron (leave alone)
camera=[];
% Handle varargin
%plot_type = 'overview';
plot_type = {};
beamdex = -1;
if nargin > 1
    i = 1;
    while i < nargin
        switch varargin{i}
            case {'overview',...
                    'xyz','xyz_therm','xyz_lost','xyz_birth','xyz_total',...
                    'xyz_initial','xyz_therm_initial','xyz_lost_initial','xyz_total_initial',...
                    'flux','flux_therm','flux_birth','flux_total',...
                    'flux_initial','flux_lost_initial','flux_therm_initial',...
                    'pitch','pitch_therm','pitch_lost', 'pitch_birth', ...
                    'pitch_initial','pitch_therm_initial','pitch_lost_initial',...
                    'dist','dist_therm','dist_lost','dist_birth', ...
                    'dist_initial','dist_therm_initial','dist_lost_initial',...
                    'orbit_rz','orbit_flux','orbit_vspace',...
                    'rho','rho_therm','rho_birth','rho_total',...
                    'rho_initial','rho_lost_initial','rho_therm_initial',...
                    'frac_loss', 'frac_therm',...
                    'distribution','heating','current','fueling',...
                    'injection','birth_image',...
                    'wall','wall_loss','wall_heat','wall_shine','benchmarks',...
                    'wall_loss_2d','wall_heat_2d','wall_shine_2d',...
                    'wall_loss_log10','wall_heat_log10','wall_shine_log10',...
                    'grid','grid_s',...
                    'camview','bmir'}
                plot_type{end+1}=varargin{i}; %Make multiple plots possible                    
            case 'beam'
                i=i+1;
                beamdex = varargin{i};
        end
        i=i+1;
    end
end
if isfield(beam_data,'vlldist')
    if beamdex == -1
        for beamdex=1:beam_data.nbeams
            figure('Color','white','Position',[1 -100 1024 768])
            subplot(2,1,1);
            nt     = size(beam_data.vlldist,3);
            t=1000.*(0:1/(nt-3):1).*beam_data.t_end(beamdex);
            colormap hot;
            temp_y = squeeze(beam_data.vllaxis);
            temp_z = squeeze(beam_data.vlldist(beamdex,:,3:end));
            maxy = max(temp_y);
            miny = min(temp_y);
            maxxy = max([ abs(maxy) abs(miny)]);
            ytext = 'V_{ll}';
            yticks = [miny 0 maxy];
            if maxxy > 1E6
                temp_y = temp_y./1E6;
                ytext = 'V_{ll} x10^6';
                yticks = [miny 0 maxy]./1E6;
            end
            pixplot(t(1:end-1),temp_y,temp_z(:,1:end-1)');
            cmap=caxis;
            caxis([0 max(cmap)*.5]);
            set(gca,'FontSize',24,'YTick',yticks,'YTickLabelMode','auto',...
                'XTickMode','auto','XTickLabelMode','auto');
            ylabel(ytext);
            title(['Distribution Function Evolution (BEAM: ' num2str(beamdex,'%2.2i') ')']);
            subplot(2,1,2);
            max_part = max(sum(temp_z));
            lost_perc = 100.*(max_part-sum(temp_z))./max_part;
            set(gca,'FontSize',24);
            if (lost_perc(1) >= 100) % Neutral Beam Injection
                plot(t(1:end-1),100.*sum(temp_z(1:end-1))./max_part)
                ylabel('Population (%)');
            else
                plot(t(1:end-1),lost_perc(1:end-1));
                ylabel('Lost (%)');
            end
            xlabel('Time (ms)');
        end
        if beam_data.nbeams == 1, return; end
        figure('Color','white','Position',[1 -100 1024 768])
        subplot(2,1,1);
        nt     = size(beam_data.vlldist,3);
        t=1000.*(0:1/(nt-3):1).*beam_data.t_end(beamdex);
        colormap hot;
        temp_y = squeeze(beam_data.vllaxis);
        temp_z = squeeze(sum(beam_data.vlldist(:,:,3:end)));
        maxy = max(temp_y);
        miny = min(temp_y);
        maxxy = max([ abs(maxy) abs(miny)]);
        ytext = 'V_{ll}';
        yticks = [miny 0 maxy];
        if maxxy > 1E6
            temp_y = temp_y./1E6;
            ytext = 'V_{ll} x10^6';
            yticks = [miny 0 maxy]./1E6;
        end
        pixplot(t(1:end-1),temp_y,temp_z(:,1:end-1)');
        cmap=caxis;
        caxis([0 max(cmap)*.5]);
        set(gca,'FontSize',24,'YTick',yticks,'YTickLabelMode','auto',...
            'XTickMode','auto','XTickLabelMode','auto');
        ylabel(ytext);
        title('Distribution Function Evolution (Total)');
        subplot(2,1,2);
        max_part = max(sum(temp_z));
        lost_perc = 100.*(max_part-sum(temp_z))./max_part;
        set(gca,'FontSize',24);
        if (lost_perc(1) >= 100) % Neutral Beam Injection
            plot(t(1:end-1),100.*sum(temp_z(1:end-1))./max_part)
            ylabel('Population (%)');
        else
            plot(t(1:end-1),lost_perc(1:end-1));
            ylabel('Lost (%)');
        end
        xlabel('Time (ms)');
    else
        figure('Color','white','Position',[1 -100 1024 768])
        subplot(2,1,1);
        nbeams = size(beam_data.vlldist,1);
        nt     = size(beam_data.vlldist,3);
        t=1000.*(0:1/(nt-1):1).*beam_data.t_end(beamdex);
        colormap hot;
        temp_y = squeeze(beam_data.vllaxis);
        temp_z = squeeze(beam_data.vlldist(beamdex,:,:));
        pixplot(t,temp_y,temp_z');
        cmap=caxis;
        caxis([0 max(cmap)*.5]);
        ylabel('V_{ll}');
        title('Distribution Function Evolution');
        subplot(2,1,2);
        plot(t,100.*(beam_data.npart(beamdex)-sum(temp_z))./beam_data.npart(beamdex));
        ylabel('Lost (%)');
        xlabel('Time (ms)');
    end
else
    % Setup indicies
    dex1 = 1;
    shine_dex = beam_data.end_state==3;
    lost_dex = beam_data.end_state==2;
    therm_dex = beam_data.end_state==1;
    orbit_dex = beam_data.end_state==0;
    born_dex = ones(beam_data.nparticles,1)==1;
    if any(shine_dex)
        dex1=3; 
        born_dex = ~shine_dex;
    end
    % Handle plotting only certain beams
    if beamdex == -1
        beamdex = 1:double(beam_data.nbeams);
    end
    beam_dex=zeros(beam_data.nparticles,1);
    for i=1:length(beamdex)
        beam_dex(beam_data.Beam==beamdex(i)) = 1;
    end
    beam_dex = beam_dex==1;
    shine_dex = and(shine_dex,beam_dex);
    lost_dex = and(lost_dex,beam_dex);
    therm_dex = and(therm_dex,beam_dex);
    orbit_dex = and(orbit_dex,beam_dex);
    born_dex = and(born_dex,beam_dex);
    % Calc Last Index for each particle
    last_dex=zeros(1,beam_data.nparticles);
    for i=1:beam_data.nparticles
        last_dex(i) = min([find(beam_data.R_lines(:,i)>0,1,'last') beam_data.npoinc+1]);
    end
    % Make plots
    for i = 1:size(plot_type,2)
    switch lower(plot_type{i})
        case 'xyz'
            x1  = beam_data.X_lines(:,orbit_dex);
            y1  = beam_data.Y_lines(:,orbit_dex);
            z1  = beam_data.Z_lines(:,orbit_dex);
            ld  = last_dex(orbit_dex);
            x=[]; y=[]; z=[];
            for i=1:length(ld)
                x = [x x1(ld(i),i)];
                y = [y y1(ld(i),i)];
                z = [z z1(ld(i),i)];
            end
            plot3(x,y,z,'.k');
            axis equal;
            axis off;
            title('Orbiting Particles');
        case 'xyz_therm'
            x1  = beam_data.X_lines(:,therm_dex);
            y1  = beam_data.Y_lines(:,therm_dex);
            z1  = beam_data.Z_lines(:,therm_dex);
            ld  = last_dex(therm_dex);
            x=[]; y=[]; z=[];
            for i=1:length(ld)
                x = [x x1(ld(i),i)];
                y = [y y1(ld(i),i)];
                z = [z z1(ld(i),i)];
            end
            plot3(x,y,z,'.k');
            axis equal;
            axis off;
            title('Thermalized Particles');
        case 'xyz_lost'
            x1  = beam_data.X_lines(:,lost_dex);
            y1  = beam_data.Y_lines(:,lost_dex);
            z1  = beam_data.Z_lines(:,lost_dex);
            ld  = last_dex(lost_dex);
            x=[]; y=[]; z=[];
            for i=1:length(ld)
                x = [x x1(ld(i),i)];
                y = [y y1(ld(i),i)];
                z = [z z1(ld(i),i)];
            end
            plot3(x,y,z,'.k');
            axis equal;
            axis off;
            title('Lost Particles');
        case 'xyz_birth'
            x  = beam_data.X_lines(dex1,born_dex);
            y  = beam_data.Y_lines(dex1,born_dex);
            z  = beam_data.Z_lines(dex1,born_dex);
            plot3(x,y,z,'.k');
            axis equal;
            axis off;
            title('Particle Birth');
        case {'overview','xyz_total'}
            leg_text={};
            x1  = beam_data.X_lines(:,therm_dex);
            y1  = beam_data.Y_lines(:,therm_dex);
            z1  = beam_data.Z_lines(:,therm_dex);
            ld  = last_dex(therm_dex);
            x=[]; y=[]; z=[];
            if ~isempty(ld)
                for i=1:length(ld)
                    x = [x x1(ld(i),i)];
                    y = [y y1(ld(i),i)];
                    z = [z z1(ld(i),i)];
                end
                plot3(x,y,z,'.b');
                hold on;
                leg_text=[leg_text; 'Thermalized'];
            end
            x1  = beam_data.X_lines(:,lost_dex);
            y1  = beam_data.Y_lines(:,lost_dex);
            z1  = beam_data.Z_lines(:,lost_dex);
            ld  = last_dex(lost_dex);
            x=[]; y=[]; z=[];
            if ~isempty(ld)
                for i=1:length(ld)
                    x = [x x1(ld(i),i)];
                    y = [y y1(ld(i),i)];
                    z = [z z1(ld(i),i)];
                end
                plot3(x,y,z,'.r');
                hold on;
                leg_text=[leg_text; 'Lost'];
            end
            x1  = beam_data.X_lines(:,orbit_dex);
            y1  = beam_data.Y_lines(:,orbit_dex);
            z1  = beam_data.Z_lines(:,orbit_dex);
            ld  = last_dex(orbit_dex);
            x=[]; y=[]; z=[];
            if ~isempty(ld)
                for i=1:length(ld)
                    x = [x x1(ld(i),i)];
                    y = [y y1(ld(i),i)];
                    z = [z z1(ld(i),i)];
                end
                plot3(x,y,z,'.k');
                hold on;
                leg_text=[leg_text; 'Orbiting'];
            end
            legend(leg_text); axis equal; title('Final State');
        case 'xyz_initial'
            x  = beam_data.X_lines(dex1,orbit_dex);
            y  = beam_data.Y_lines(dex1,orbit_dex);
            z  = beam_data.Z_lines(dex1,orbit_dex);
            plot3(x,y,z,'.k');
            axis equal;
            axis off;
            title('Initial Orbiting Particles');
        case 'xyz_therm_initial'
            x  = beam_data.X_lines(dex1,therm_dex);
            y  = beam_data.Y_lines(dex1,therm_dex);
            z  = beam_data.Z_lines(dex1,therm_dex);
            plot3(x,y,z,'.k');
            axis equal;
            axis off;
            title('Initial Thermalized Particles');
        case 'xyz_lost_initial'
            x  = beam_data.X_lines(dex1,lost_dex);
            y  = beam_data.Y_lines(dex1,lost_dex);
            z  = beam_data.Z_lines(dex1,lost_dex);
            plot3(x,y,z,'.k');
            axis equal;
            axis off;
            title('Initial Lost Particles');
        case {'xyz_total_initial'}
            leg_text={};
            x  = beam_data.X_lines(dex1,therm_dex);
            y  = beam_data.Y_lines(dex1,therm_dex);
            z  = beam_data.Z_lines(dex1,therm_dex);
            if ~isempty(x)
                plot3(x,y,z,'.b');
                hold on;
                leg_text=[leg_text; 'Thermalized'];
            end
            x  = beam_data.X_lines(dex1,lost_dex);
            y  = beam_data.Y_lines(dex1,lost_dex);
            z  = beam_data.Z_lines(dex1,lost_dex);
            if ~isempty(x)
                plot3(x,y,z,'.r');
                hold on;
                leg_text=[leg_text; 'Lost'];
            end
            x  = beam_data.X_lines(dex1,orbit_dex);
            y  = beam_data.Y_lines(dex1,orbit_dex);
            z  = beam_data.Z_lines(dex1,orbit_dex);
            if ~isempty(x)
                plot3(x,y,z,'.k');
                hold on;
                leg_text=[leg_text; 'Orbiting'];
            end
            legend(leg_text); axis equal; title('Initial State');
        case 'flux'
            u   = beams3d_fixUlines(beam_data);
            x1  = u(:,orbit_dex);
            y1  = beam_data.S_lines(:,orbit_dex);
            ld  = last_dex(orbit_dex);
            x=[]; y=[];
            for i=1:length(ld)
                x = [x x1(ld(i),i)];
                y = [y y1(ld(i),i)];
            end
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Orbiting Particles');
        case 'flux_therm'
            u   = beams3d_fixUlines(beam_data);
            x1  = u(:,therm_dex);
            y1  = beam_data.S_lines(:,therm_dex);
            ld  = last_dex(therm_dex);
            x=[]; y=[];
            for i=1:length(ld)
                x = [x x1(ld(i),i)];
                y = [y y1(ld(i),i)];
            end
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Thermalized Particles');
        case 'flux_birth'
            u   = beams3d_fixUlines(beam_data);
            x  = u(dex1,born_dex);
            y  = beam_data.S_lines(dex1,born_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Born Particles');
        case 'flux_initial'
            u   = beams3d_fixUlines(beam_data);
            x  = u(dex1,orbit_dex);
            y  = beam_data.S_lines(dex1,orbit_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Intial Orbiting Particles');
        case 'flux_lost_initial'
            u   = beams3d_fixUlines(beam_data);
            x  = u(dex1,lost_dex);
            y  = beam_data.S_lines(dex1,lost_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Intial Lost Particles');
        case 'flux_therm_initial'
            u   = beams3d_fixUlines(beam_data);
            x  = u(dex1,therm_dex);
            y  = beam_data.S_lines(dex1,therm_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Intial Thermalized Particles');
        case 'rho'
            cdex=pchip([0 1],[1 0 0; 0 1 0]',0:1./double(beam_data.npoinc):1)';
            colororder(cdex);
            for j=1:beam_data.npoinc+1
                x1 =  sqrt(beam_data.S_lines(j,orbit_dex));
                h  = 1.2/64.0;
                edges = 0:h:1.2;
                rho   = 0.5.*(edges(1:end-1)+edges(2:end));
                y=zeros(1,64);
                for i=1:length(edges)-1
                    sm = x1 > edges(i);
                    sp = x1 <= edges(i+1);
                    y(i) = sum(and(sm,sp));
                end
                hold on;
                plot(rho,y);
                hold off;
            end
            xlim([0 1.2]);
            xlabel('r/a');
            ylabel('# Particles');
            title('Orbiting Particles');
        case 'rho_therm'
            cdex=pchip([0 1],[1 0 0; 0 1 0]',0:1./double(beam_data.npoinc):1)';
            colororder(cdex);
            leg_text={};
            dt = double(beam_data.t_end(1))./double(beam_data.npoinc);
            for j=1:beam_data.npoinc+1
                leg_text{j} = ['t=' num2str(dt.*double(j-1).*1000,'%3.0f')];
                x1 =  sqrt(beam_data.S_lines(j,therm_dex));
                h  = 1.2/64.0;
                edges = 0:h:1.2;
                rho   = 0.5.*(edges(1:end-1)+edges(2:end));
                y=zeros(1,64);
                for i=1:length(edges)-1
                    sm = x1 > edges(i);
                    sp = x1 <= edges(i+1);
                    y(i) =sum(and(sm,sp));
                end
                plot(rho,y);
                hold on;
            end
            hold off;
            xlim([0 1.2]);
            xlabel('r/a');
            ylabel('# Particles');
            title('Thermalized Particles');
            legend(leg_text);
        case 'rho_birth'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            x1 =  sqrt(abs(beam_data.S_lines(dex1,born_dex)));
            w = double(beam_data.Weight(born_dex));
            nres = 64;
            y_size=[0 1.2];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, x1, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w,[NumBins_y 1]);
            plot(ybins,vals,'k','LineWidth',4);
            xlim([0 1.2]);
            set(gca,'FontSize',24);
            xlabel('r/a');
            ylabel('# Particles');
            title('Born Particles');
        case 'rho_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            x1 =  sqrt(abs(beam_data.S_lines(dex1,orbit_dex)));
            w = double(beam_data.Weight(orbit_dex));
            nres = 64;
            y_size=[0 1.2];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, x1, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w,[NumBins_y 1]);
            plot(ybins,vals,'k','LineWidth',4);
            xlim([0 1.2]);
            set(gca,'FontSize',24);
            xlabel('r/a');
            ylabel('# Particles');
            title('Intial Orbiting Particles');
        case 'rho_lost_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            x1 =  sqrt(abs(beam_data.S_lines(dex1,lost_dex)));
            w = double(beam_data.Weight(lost_dex));
            nres = 64;
            y_size=[0 1.2];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, x1, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w,[NumBins_y 1]);
            plot(ybins,vals,'k','LineWidth',4);
            xlim([0 1.2]);
            set(gca,'FontSize',24);
            xlabel('r/a');
            ylabel('# Particles');
            title('Intial Lost Particles');
        case 'rho_therm_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            x1 =  sqrt(abs(beam_data.S_lines(dex1,therm_dex)));
            w = double(beam_data.Weight(therm_dex));
            nres = 64;
            y_size=[0 1.2];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, x1, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w,[NumBins_y 1]);
            plot(ybins,vals,'k','LineWidth',4);
            xlim([0 1.2]);
            set(gca,'FontSize',24);
            xlabel('r/a');
            ylabel('# Particles');
            title('Intial Thermalized Particles');
        case 'pitch'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            w = double(beam_data.Weight(orbit_dex));
            vperp = vperp(:,orbit_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(orbit_dex);
            p2    = length(ld);
            for i=1:length(ld)
                p2(i) = pitch(ld(i),i);
            end
            nres = 128;
            y_size=[min(p2) max(p2)];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w);
            plot(ybins,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Orbiting Distribution');
        case 'pitch_therm'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            w = double(beam_data.Weight(therm_dex));
            vperp = vperp(:,therm_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(therm_dex);
            p2    = length(ld);
            for i=1:length(ld)
                p2(i) = pitch(ld(i),i);
            end
            nres = 128;
            y_size=[min(p2) max(p2)];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w);
            plot(ybins,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Thermalized Particle Pitch');
        case 'pitch_lost'
            disp('BROKEN')
            return;
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            w = double(beam_data.Weight(lost_dex));
            vperp = vperp(:,lost_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(lost_dex);
            p2    = length(ld);
            nres = 128;
            y_size=[min(p2) max(p2)];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w);
            plot(ybins,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Lost Particle Pitch');
        case 'pitch_birth'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(dex1,born_dex);
            w = double(beam_data.Weight(born_dex));
            vperp = vperp(dex1,born_dex);
            pitch = atan2d(vll,vperp);
            p2    = pitch(1,:);
            nres = 128;
            y_size=[min(p2) max(p2)];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w);
            plot(ybins,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Birth Distribution');
        case 'pitch_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            w = double(beam_data.Weight(orbit_dex));
            vperp = vperp(:,orbit_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(orbit_dex);
            p2    = pitch(1,:);
            nres = 128;
            y_size=[min(p2) max(p2)];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w);
            plot(ybins,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Initial Orbiting Distribution');
        case 'pitch_therm_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            w = double(beam_data.Weight(therm_dex));
            vperp = vperp(:,therm_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(therm_dex);
            p2    = pitch(1,:);
            nres = 128;
            y_size=[min(p2) max(p2)];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w);
            plot(ybins,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Initial Thermalized Distribution');
        case 'pitch_lost_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            w = double(beam_data.Weight(lost_dex));
            vperp = vperp(:,lost_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(lost_dex);
            p2    = pitch(1,:);
            nres = 128;
            y_size=[min(p2) max(p2)];
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_y = numel(ybins);
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray(Yi', w);
            plot(ybins,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Initial Lost Distribution');
        case 'dist'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            w = double(beam_data.Weight(orbit_dex));
            vperp = vperp(:,orbit_dex);
            ld    = last_dex(orbit_dex);
            p1    = zeros(1,length(ld));
            p2    = zeros(1,length(ld));
            for i=1:length(ld)
                p1(i) = vll(ld(i),i);
                p2(i) = vperp(ld(i),i);
            end
            nres = 256;
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            xbins = x_size(1):diff(x_size)./nres:x_size(2);
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_x = numel(xbins);
            NumBins_y = numel(ybins);
            Xi = round( interp1(xbins, 1:NumBins_x, p1, 'linear', 'extrap') );
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Xi = max( min(Xi,NumBins_x), 1);
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray([Xi; Yi]', w, [NumBins_x NumBins_y]);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'Mm/s';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'km/s';
            else
                norm = 1;
                units = 'm/s';
            end
            colormap jet;
            ha=pcolor(xbins./norm,ybins./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Oribiting Distribution');
        case 'dist_therm'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            w = double(beam_data.Weight(therm_dex));
            vperp = vperp(:,therm_dex);
            ld    = last_dex(therm_dex);
            p1    = zeros(1,length(ld));
            p2    = zeros(1,length(ld));
            for i=1:length(ld)
                p1(i) = vll(ld(i),i);
                p2(i) = vperp(ld(i),i);
            end
            nres = 256;
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            xbins = x_size(1):diff(x_size)./nres:x_size(2);
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_x = numel(xbins);
            NumBins_y = numel(ybins);
            Xi = round( interp1(xbins, 1:NumBins_x, p1, 'linear', 'extrap') );
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Xi = max( min(Xi,NumBins_x), 1);
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray([Xi; Yi]', w, [NumBins_x NumBins_y]);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'Mm/s';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'km/s';
            else
                norm = 1;
                units = 'm/s';
            end
            colormap jet;
            ha=pcolor(xbins./norm,ybins./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Thermalized Distribution');
        case 'dist_lost'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            w = double(beam_data.Weight(lost_dex));
            vperp = vperp(:,lost_dex);
            ld    = last_dex(lost_dex);
            p1    = zeros(1,length(ld));
            p2    = zeros(1,length(ld));
            for i=1:length(ld)
                p1(i) = vll(ld(i),i);
                p2(i) = vperp(ld(i),i);
            end
            nres = 256;
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            xbins = x_size(1):diff(x_size)./nres:x_size(2);
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_x = numel(xbins);
            NumBins_y = numel(ybins);
            Xi = round( interp1(xbins, 1:NumBins_x, p1, 'linear', 'extrap') );
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Xi = max( min(Xi,NumBins_x), 1);
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray([Xi; Yi]', w, [NumBins_x NumBins_y]);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'Mm/s';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'km/s';
            else
                norm = 1;
                units = 'm/s';
            end
            colormap jet;
            ha=pcolor(xbins./norm,ybins./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Lost Distribution');
        case 'dist_birth'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(dex1,orbit_dex);
            w = double(beam_data.Weight(orbit_dex));
            vperp = vperp(dex1,orbit_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            nres = 256;
            vmax = sqrt(max(p1.*p1+p2.*p2));
            x_size=[-1 1].*max(beam_data.partvmax,vmax);
            y_size=[ 0 1].*max(beam_data.partvmax,vmax);
            xbins = x_size(1):diff(x_size)./nres:x_size(2);
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_x = numel(xbins);
            NumBins_y = numel(ybins);
            Xi = round( interp1(xbins, 1:NumBins_x, p1, 'linear', 'extrap') );
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Xi = max( min(Xi,NumBins_x), 1);
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray([Xi; Yi]', w, [NumBins_x NumBins_y]);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'Mm/s';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'km/s';
            else
                norm = 1;
                units = 'm/s';
            end
            colormap jet;
            ha=pcolor(xbins./norm,ybins./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Birth Distribution');
        case 'dist_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            w = double(beam_data.Weight(orbit_dex));
            vperp = vperp(:,orbit_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            nres = 256;
            vmax = sqrt(max(p1.*p1+p2.*p2));
            x_size=[-1 1].*max(beam_data.partvmax,vmax);
            y_size=[ 0 1].*max(beam_data.partvmax,vmax);
            xbins = x_size(1):diff(x_size)./nres:x_size(2);
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_x = numel(xbins);
            NumBins_y = numel(ybins);
            Xi = round( interp1(xbins, 1:NumBins_x, p1, 'linear', 'extrap') );
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Xi = max( min(Xi,NumBins_x), 1);
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray([Xi; Yi]', w, [NumBins_x NumBins_y]);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'Mm/s';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'km/s';
            else
                norm = 1;
                units = 'm/s';
            end
            colormap jet;
            ha=pcolor(xbins./norm,ybins./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Initial Oribiting Distribution');
        case 'dist_therm_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            w = double(beam_data.Weight(therm_dex));
            vperp = vperp(:,therm_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            nres = 256;
            vmax = sqrt(max(p1.*p1+p2.*p2));
            x_size=[-1 1].*max(beam_data.partvmax,vmax);
            y_size=[ 0 1].*max(beam_data.partvmax,vmax);
            xbins = x_size(1):diff(x_size)./nres:x_size(2);
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_x = numel(xbins);
            NumBins_y = numel(ybins);
            Xi = round( interp1(xbins, 1:NumBins_x, p1, 'linear', 'extrap') );
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Xi = max( min(Xi,NumBins_x), 1);
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray([Xi; Yi]', w, [NumBins_x NumBins_y]);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'Mm/s';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'km/s';
            else
                norm = 1;
                units = 'm/s';
            end
            colormap jet;
            ha=pcolor(xbins./norm,ybins./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Initial Thermalized Distribution');
        case 'dist_lost_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            w = double(beam_data.Weight(lost_dex));
            vperp = vperp(:,lost_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            nres = 256;
            vmax = sqrt(max(p1.*p1+p2.*p2));
            x_size=[-1 1].*max(beam_data.partvmax,vmax);
            y_size=[ 0 1].*max(beam_data.partvmax,vmax);
            xbins = x_size(1):diff(x_size)./nres:x_size(2);
            ybins = y_size(1):diff(y_size)./nres:y_size(2);
            NumBins_x = numel(xbins);
            NumBins_y = numel(ybins);
            Xi = round( interp1(xbins, 1:NumBins_x, p1, 'linear', 'extrap') );
            Yi = round( interp1(ybins, 1:NumBins_y, p2, 'linear', 'extrap') );
            Xi = max( min(Xi,NumBins_x), 1);
            Yi = max( min(Yi,NumBins_y), 1);
            vals = accumarray([Xi; Yi]', w, [NumBins_x NumBins_y]);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'Mm/s';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'km/s';
            else
                norm = 1;
                units = 'm/s';
            end
            colormap jet;
            ha=pcolor(xbins./norm,ybins./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Initial Lost Distribution');
        case 'orbit_rz'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            R = beam_data.R_lines(:,beam_dex);
            Z = beam_data.Z_lines(:,beam_dex);
            plot(R,Z);
            axis equal;
            set(gca,'FontSize',24);
            xlabel('R [m]');
            ylabel('Z [m]');
            title('Particle Orbits');
        case 'orbit_rz_lost'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            R = beam_data.R_lines(:,lost_dex);
            Z = beam_data.Z_lines(:,lost_dex);
            plot(R,Z);
            axis equal;
            set(gca,'FontSize',24);
            xlabel('R [m]');
            ylabel('Z [m]');
            title('Lost Particle Orbits');
        case 'orbit_flux'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            S = beam_data.S_lines(:,beam_dex);
            U = beams3d_fixUlines(beam_data);
            U = U(:,beam_dex);
            polarplot(U,sqrt(S),'.');
            set(gca,'FontSize',24);
            title('Particle Orbits (rho,u)');
        case 'orbit_flux_lost'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            S = beam_data.S_lines(:,lost_dex);
            U = beams3d_fixUlines(beam_data);
            U = U(:,lost_dex);
            polarplot(U,sqrt(S),'.');
            set(gca,'FontSize',24);
            title('Lost Particle Orbits (rho,u)');
        case 'orbit_vspace'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vll = beam_data.vll_lines(:,beam_dex);
            vperp = beams3d_calc_vperp(beam_data);
            vperp = vperp(:,beam_dex);
            if beam_data.partvmax>1E6
                units='x1000 [km/s]';
                factor = 1E-6;
            elseif beam_data.partvmax>1E3
                units='[km/s]';
                factor = 1E-3;
            else
                units='[m/s]';
                factor=1;
            end
            plot(vll.*factor,vperp.*factor,'.');
            xlim([-1 1].*beam_data.partvmax.*factor);
            ylim([0 1].*beam_data.partvmax.*factor);
            set(gca,'FontSize',24);
            title('Particle Orbits Velocity Space');
            xlabel(['Parallel Velocity (v_{||}) ' units ]);
            ylabel(['Perpendicuarl Velocity (v_\perp) ' units ]);
        case 'distribution'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            n1 = double(beam_data.ns_prof4-1);
            n2 = double(beam_data.ns_prof5-1);
            h1 = 2*beam_data.partvmax./n1;
            h2 = beam_data.partvmax./n2;
            x1 = -beam_data.partvmax:h1:beam_data.partvmax;
            x2 = 0:h2:beam_data.partvmax;
            if beam_data.partvmax>1E6
                units='x1000 [km/s]';
                factor = 1E-6;
            elseif beam_data.partvmax>1E3
                units='[km/s]';
                factor = 1E-3;
            else
                units='[m/s]';
                factor=1;
            end
            pixplot(x1.*factor,x2.*factor,squeeze(sum(beam_data.dist2d_prof(beamdex,:,:),1)));
            vmax = round(beam_data.partvmax.*factor);
            xtick = -vmax:1:vmax;
            ytick = 0:1:vmax;
            set(gca,'FontSize',24,'XTick',xtick,'XTickLabelMode','auto',...
                'YTick',ytick,'YTickLabelMode','auto');
            title('BEAMS3D Distribution Function');
            xlabel(['Parallel Velocity (v_{||}) ' units ]);
            ylabel(['Perpendicuarl Velocity (v_\perp) ' units ]);
            colormap jet;
        case 'heating'
            h = 1./double(beam_data.ns_prof1-1);
            rho = 0:h:1;
            imax = 2.*max(sum(beam_data.ipower_prof));
            if imax>1E6
                factor=1e-6;
                units = '[MW/m^3]';
            elseif imax>1E3
                factor=1e-3;
                units = '[kW/m^3]';
            else
                factor = 1;
                units = '[W/m^3]';
            end
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            plot(rho,sum(beam_data.epower_prof(beamdex,:),1).*factor,'b','LineWidth',4); hold on;
            plot(rho,sum(beam_data.ipower_prof(beamdex,:),1).*factor,'--r','LineWidth',4);
            ylim([0 imax.*factor]);
            set(gca,'FontSize',24);
            legend('Electrons','Ions');
            xlabel('Effective Minor Radius (r/a) [norm]');
            ylabel(['Power Deposition ' units]);
            title('BEAMS3D Total Power Deposition');
        case 'current'
            h = 1./double(beam_data.ns_prof1-1);
            rho = 0:h:1;
            jmax = 2.*max(sum(beam_data.j_prof));
            if jmax>1E6
                factor=1e-6;
                units = '[MA/m^2]';
            elseif jmax>1E3
                factor=1e-3;
                units = '[kA/m^2]';
            else
                factor = 1;
                units = '[A/m^2]';
            end
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            plot(rho,sum(beam_data.j_prof(beamdex,:),1).*factor,'k','LineWidth',4);
            ylim([0 jmax.*factor]);
            set(gca,'FontSize',24);
            xlabel('Effective Minor Radius (r/a) [norm]');
            ylabel(['Current Density ' units]);
            title('BEAMS3D Total Current Density (uncorr.)');
        case 'fueling'
            h = 1./double(beam_data.ns_prof1-1);
            rho = 0:h:1;
            jmax = 2.*max(sum(beam_data.ndot_prof));
            if jmax>1E15
                factor=1e-15;
                units = 'x10^{15} [part/(s*m^3)]';
            elseif jmax>1E12
                factor=1e-12;
                units = 'x10^{12} [part/(s*m^3)]';
            elseif jmax>1E9
                factor=1e-9;
                units = 'x10^9 [part/(s*m^3)]';
            elseif jmax>1E6
                factor=1e-6;
                units = 'x10^6 [part/(s*m^3)]';
            elseif jmax>1E3
                factor=1e-3;
                units = 'x10^3 [part/(s*m^3)]';
            else
                factor = 1;
                units = '[part/(s*m^3)]';
            end
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            plot(rho,sum(beam_data.ndot_prof(beamdex,:),1).*factor,'k','LineWidth',4);
            ylim([0 jmax.*factor]);
            set(gca,'FontSize',24);
            xlabel('Effective Minor Radius (r/a) [norm]');
            ylabel(['Fueling ' units]);
            title('BEAMS3D Fast-Ion Fueling');
        case 'injection'
            nmax=512;
            x_shine = beam_data.X_lines([1 2],shine_dex);
            y_shine = beam_data.Y_lines([1 2],shine_dex);
            z_shine = beam_data.Z_lines([1 2],shine_dex);
            n = size(x_shine,2);
            dn = max([round(n/nmax) 1]);
            dl = 1:dn:n;
            x_shine = x_shine(:,dl);
            y_shine = y_shine(:,dl);
            z_shine = z_shine(:,dl);
            temp_dex = and(~shine_dex,beam_dex);
            x = beam_data.X_lines([1 2],temp_dex);
            y = beam_data.Y_lines([1 2],temp_dex);
            z = beam_data.Z_lines([1 2],temp_dex);
            n = size(x,2);
            dn = max([round(n/nmax) 1]);
            dl = 1:dn:n;
            x = x(:,dl);
            y = y(:,dl);
            z = z(:,dl);
            leg_text={};
            hold on;
            if ~isempty(x_shine), plot3(x_shine,y_shine,z_shine,'r'); leg_text=[leg_text; 'Shinethrough']; end
            if ~isempty(x), plot3(x,y,z,'g'); leg_text=[leg_text; 'Injected'];end
            legend(leg_text);
        case 'birth_image'
            r_nb = [];
            z_nb = [];
            temp_dex = and(~shine_dex,beam_dex);
            r_nb = beam_data.R_lines(3,temp_dex);
            z_nb = beam_data.Z_lines(3,temp_dex);
            nres=0.01;
            x_size=[min(beam_data.raxis) max(beam_data.raxis)];
            y_size=[min(beam_data.zaxis) max(beam_data.zaxis)];
            edges={x_size(1):nres:x_size(2) ...
                y_size(1):nres:y_size(2)};
            vals = hist3([r_nb; z_nb]',edges);
            ha=pcolor(edges{1},edges{2},vals');
            set(ha,'LineStyle','none');
            colormap hot;
            caxis([0 mean(max(vals))])
            set(gcf,'Color','white','Position',[1 -100 1024 768]);
            set(gca,'FontSize',24);
            xlabel('R [m]');
            ylabel('Z [m]');
         case 'wall'
             output_args{1}=patch('Vertices',beam_data.wall_vertex,'Faces',beam_data.wall_faces,'FaceVertexCData',zeros(beam_data.nvertex,3));
             set(output_args{1},'EdgeColor','black','FaceColor','blue','FaceAlpha',0.33);
        case {'wall_loss','wall_shine','wall_heat'}
            switch lower(plot_type{i})
                case 'wall_loss'
                    val=beam_data.wall_strikes;
                case 'wall_shine'
                    val=sum(beam_data.wall_shine(beamdex,:),1)';
                case 'wall_heat'
                    val=sum(beam_data.wall_load(beamdex,:),1)';
            end
            output_args{1}=patch('Vertices',beam_data.wall_vertex,'Faces',beam_data.wall_faces,'FaceVertexCData',val,'LineStyle','none','CDataMapping','scaled','FaceColor','flat');
            cmap = colormap('hot');
            cmap(1,:) = [0.2 0.2 0.2]; % grey
            colormap(cmap);
        case {'wall_loss_log10','wall_shine_log10','wall_heat_log10'}
            switch lower(plot_type)
                case 'wall_loss_log10'
                    val=beam_data.wall_strikes;
                case 'wall_shine_log10'
                    val=sum(beam_data.wall_shine(beamdex,:),1)';
                case 'wall_heat_log10'
                    val=sum(beam_data.wall_load(beamdex,:),1)';
            end
            val = log10(val);
            output_args{1}=patch('Vertices',beam_data.wall_vertex,'Faces',beam_data.wall_faces,'FaceVertexCData',val,'LineStyle','none','CDataMapping','scaled','FaceColor','flat');
            cmap = colormap('hot');
            cmap(1,:) = [0.2 0.2 0.2]; % grey
            colormap(cmap);
        case {'wall_heat_2d','wall_shine_2d','wall_loss_2d'}
            verts = beam_data.wall_vertex;
            faces = beam_data.wall_faces;
            d1 = faces(:,1); d2=faces(:,2); d3=faces(:,3);
            x=verts(:,1); y=verts(:,2); z=verts(:,3); 
            phi = atan2(y,x); phi(phi<0) = phi(phi<0) + 2*pi;
            r = sqrt(x.*x+y.*y);
            [raxis, zaxis] = beams3d_magaxis(beam_data);
            rax = pchip(beam_data.phiaxis,raxis,mod(phi,max(beam_data.phiaxis)));
            zax = pchip(beam_data.phiaxis,zaxis,mod(phi,max(beam_data.phiaxis)));
            a = r - rax;
            b = z - zax;
            th = atan2(b,a);
            x = (phi(d1)+phi(d2)+phi(d3))./3;
            y = (th(d1)+th(d2)+th(d3))./3;
            factor=1;
            clabel='hits';
            switch lower(plot_type{i})
                case 'wall_heat_2d'
                    val = sum(beam_data.wall_load(beamdex,:))';
                    if max(val) > 1E6
                        factor = 1E6;
                        clabel = 'MW/m^2';
                    elseif max(val) > 1E3
                        factor = 1E3;
                        clabel = 'kW/m^2';
                    else
                        factor = 1E3;
                        clabel = 'W/m^2';
                    end
                case 'wall_shine_2d'
                    val = sum(beam_data.wall_shine(beamdex,:))';
                    if max(val) > 1E6
                        factor = 1E6;
                        clabel = 'MW/m^2';
                    elseif max(val) > 1E3
                        factor = 1E3;
                        clabel = 'kW/m^2';
                    else
                        factor = 1E3;
                        clabel = 'W/m^2';
                    end
                case 'wall_loss_2d'
                    val = beam_data.wall_strikes;
            end
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            F = scatteredInterpolant(x,y,val);
            xp = 0:2*pi./256:2*pi;
            yp = -pi:2*pi./256:pi;
            Fp = F({xp,yp});
            pixplot(xp,yp,Fp./factor);
            colormap hot;
            caxis([0 max(caxis)]);
            h = max(beam_data.phiaxis);
            %xtick=0:h:2*pi;
            x0 = h/2;
            xf = 2*pi;
            xtick=x0:h:xf;
            xticklabel={};
            for i=1:length(xtick)
                xticklabel{i} = num2str(i,'%i');
            end
            for i=1:length(xtick)-1
                hold on;
                plot([i i].*h,ylim,'w');
            end
            set(gca,'YTick',[-pi/2 0 pi/2],'YTickLabel',{'Lower','Outboard','Upper'},...
                'XTick',xtick,'XTickLabel',xticklabel,'FontSize',24);
            xlabel('Field Period');
            ylabel('Poloidal Extent');
            title('BEAMS3D 2D Wall Losses');
            ha = colorbar;
            ylabel(ha,clabel);
        case 'grid'
            x=[]; y=[]; z=[];
            raxis = beam_data.raxis;
            paxis = beam_data.phiaxis;
            zaxis = beam_data.zaxis;
            for k = 1:beam_data.nphi
                for j = 1:beam_data.nz
                    x = [x; raxis(1)*cos(paxis(k)) raxis(end)*cos(paxis(k))];
                    y = [y; raxis(1)*sin(paxis(k)) raxis(end)*sin(paxis(k))];
                    z = [z; zaxis(j) zaxis(j)];
                end
                for j = 1:beam_data.nr
                    x = [x; raxis(j)*cos(paxis(k)) raxis(j)*cos(paxis(k))];
                    y = [y; raxis(j)*sin(paxis(k)) raxis(j)*sin(paxis(k))];
                    z = [z; zaxis(1) zaxis(end)];
                end
            end
            plot3(x',y',z','k');
        case 'grid_s'
            x=[]; y=[]; z=[]; s=[];
            raxis = beam_data.raxis;
            paxis = beam_data.phiaxis;
            zaxis = beam_data.zaxis;
            
            for k = 1:beam_data.nphi
                for j = 1:beam_data.nz
                    x = [x; raxis.*cos(paxis(k))];
                    y = [y; raxis.*sin(paxis(k))];
                    z = [z; raxis*0.0+zaxis(j)];
                    s = [s; beam_data.S_ARR(:,k,j)];
                end
            end
            scatter3(x,y,z,s.*0.0+1,s,'o');
        case 'frac_loss'
            ftotal = sum(beam_data.Weight(beam_dex));
            w   = repmat(beam_data.Weight',[beam_data.npoinc+1,1]);
            endcond = repmat(beam_data.end_state',[beam_data.npoinc+1,1]);
            dex = beam_data.S_lines;
            dex(dex<=1) = 0;
            dex(dex>1)  = 1;
            dex(endcond~=2) = 0; %Only consider wall hits
            dex = dex.*w;
            f = sum(dex(:,beam_dex)');
            tend = max(beam_data.t_end(beam_dex));
            t = 0:tend./double(beam_data.npoinc):tend;
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            units = '[ms]'; factor =1E3;
            if (tend < 1E-6)
                units = '[ns]'; factor = 1E9;
            elseif (tend<1E-3)
                units = '[µs]'; factor = 1E6;
            end
            plot(t.*factor,100.*f./ftotal,'k','LineWidth',4);
            xlabel(['Time ' units]);
            ylabel('Percentage Lost [%]');
            title('Loss Fraction Evolution');
            ylim([0 min(1.2*max(ylim),100)]);
            set(gca,'FontSize',24);
        case 'frac_therm'
            ftotal = sum(beam_data.Weight(beam_dex));
            w   = repmat(beam_data.Weight',[beam_data.npoinc+1,1]);
            endcond = repmat(beam_data.end_state',[beam_data.npoinc+1,1]);
            dex = beam_data.S_lines;
            dex(dex<=1) = 0;
            dex(dex>1)  = 1;
            dex(endcond~=1) = 0; %Only consider thermalized markers
            dex = dex.*w;
            f = sum(dex(:,beam_dex)');
            tend = max(beam_data.t_end(beam_dex));
            t = 0:tend./double(beam_data.npoinc):tend;
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            units = '[ms]'; factor =1E3;
            if (tend < 1E-6)
                units = '[ns]'; factor = 1E9;
            elseif (tend<1E-3)
                units = '[µs]'; factor = 1E6;
            end
            plot(t.*factor,100.*f./ftotal,'k','LineWidth',4);
            xlabel(['Time ' units]);
            ylabel('Percentage Thermalized [%]');
            title('Thermalized Fraction Evolution');
            ylim([0 min(1.2*max(ylim),100)]);
            set(gca,'FontSize',24);            
        case 'bmir'
            vperp = beams3d_calc_vperp(beam_data);
            B_mir = beam_data.B_lines.*(vperp.^2+beam_data.vll_lines.^2)./(vperp.^2);
            B    = sqrt(beam_data.B_R.^2+beam_data.B_PHI.^2+beam_data.B_Z.^2);
            mask = beam_data.S_ARR(:,1,:)<1;
            temp = B(:,1,:);
            B0_min = min(temp(mask),[],'all');
            B0_max = max(temp(mask),[],'all');
            n=round(beam_data.nphi/2);
            mask = beam_data.S_ARR(:,n,:)<1;
            temp = B(:,n,:);
            B1_min = min(temp(mask),[],'all');
            B1_max = max(temp(mask),[],'all');
            edges=min([B0_min B1_min]):0.01:max(B_mir(dex1,beam_dex));
            [N,edges] = histcounts(B_mir(dex1,born_dex),edges);
            fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardcopy','off');
            centers = 0.5.*(edges(1:end-1)+edges(2:end));
            plot(centers,N,'k','LineWidth',4);
            axis tight;
            hold on;
            y=ylim;
            fill([B0_min B0_max B0_max B0_min],[y(1) y(1) y(2) y(2)],'red','LineStyle','none','FaceAlpha','0.3')
            fill([B1_min B1_max B1_max B1_min],[y(1) y(1) y(2) y(2)],'blue','LineStyle','none','FaceAlpha','0.3')
            set(gca,'FontSize',24);
            xlabel('B [T]');
            ylabel('Marker Count');
            title('Mirror Magnetic Field')
            legend('B_{mirror}=E/\mu','|B| \phi=0','|B| \phi=\phi_{max}/2');
        case 'camview'
            if isempty(camera), camera=[1024 768]; end
            x_cam = campos;
            a_cam = camva;
            u_cam = camup;
            t_cam = camtarget;
            n_cam = t_cam-x_cam;
            n_cam = n_cam./sqrt(sum(n_cam.*n_cam));
            syn = zeros(camera);
            X   = beam_data.X_lines(:,2);
            Y   = beam_data.Y_lines(:,2);
            Z   = beam_data.Z_lines(:,2);
            [x_temp] = points_to_camera(X(2:end),Y(2:end),Z(2:end),...
                'camera',camera,...
                'fov',a_cam,'camera_pos',x_cam,'camera_normal',n_cam,...
                'camera_up',u_cam);
            x_im = x_temp(:,1);
            y_im = x_temp(:,2);
            dex   = x_im <= camera(1);
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = y_im <= camera(2);
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = x_im >= 1;
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            dex   = y_im >= 1;
            x_im  = x_im(dex);
            y_im  = y_im(dex);
            x_max = max(x_im);
            y_max = max(y_im);
            x_min = min(x_im);
            y_min = min(y_im);
            syn_temp=hist3([x_im y_im],'nbins',[round(x_max-x_min) round(y_max-y_min)]);
            xb=linspace(x_min,x_max,size(syn_temp,1));
            yb=linspace(y_min,y_max,size(syn_temp,2));
            syn(round(xb),round(yb))=syn_temp./double(i)+syn(round(xb),round(yb));
            
            % Smooth
            sigma = 1; % set sigma to the value you need
            sz = 2*ceil(2.6 * sigma) + 1; % See note below
            mask = fspecial('gauss', sz, sigma);
            syn2 = conv2(syn, mask, 'same');
            syn=syn2; syn2=[];
            % New Stuff
            pixplot(syn)
            caxis([0 max(mean(syn))]);
            set(gcf,'Units','pixels','Position',[1 1 camera]);
            set(gca,'Units','pixels','Color','black','Position',[1 1 camera]);
            xlim([1 camera(1)]);
            ylim([1 camera(2)]);
            colormap hot;
        case 'benchmarks'
            % First make sorting arrays
            for i = 1:beam_data.nbeams
                beam_dex(:,i) = beam_data.Beam==i;
            end
            depo_dex = beam_data.neut_lines(2,:) == 0; % Particles which neutralize
            % R deposition
            figure('Position',[1 -100 1024 768],'Color','white');
            cmap = colormap('lines');
            for i=1:beam_data.nbeams
                dex = and(depo_dex',beam_dex(:,i));
                [counts,bins]=hist(beam_data.R_lines(2,dex),100);
                hold on;
                plot(bins,counts,'color',cmap(i,:));
                hold off;
            end
            set(gca,'FontSize',24);
            xlabel('R [m]');
            ylabel('# of Particles');
            title('Radial Birth Distribtuion');
            % Z deposition
            figure('Position',[1 -100 1024 768],'Color','white');
            cmap = colormap('lines');
            for i=1:beam_data.nbeams
                dex = and(depo_dex',beam_dex(:,i));
                [counts,bins]=hist(beam_data.Z_lines(2,dex),100);
                hold on;
                plot(bins,counts,'color',cmap(i,:));
                hold off;
            end
            set(gca,'FontSize',24);
            xlabel('Z [m]');
            ylabel('# of Particles');
            title('Vertical Birth Distribtuion');
            % PHI deposition
            figure('Position',[1 -100 1024 768],'Color','white');
            cmap = colormap('lines');
            for i=1:beam_data.nbeams
                dex = and(depo_dex',beam_dex(:,i));
                [counts,bins]=hist(beam_data.PHI_lines(2,dex),100);
                hold on;
                plot(bins,counts,'color',cmap(i,:));
                hold off;
            end
            set(gca,'FontSize',24);
            xlabel('Phi [rad]');
            ylabel('# of Particles');
            title('Toroidal Birth Distribtuion');
            % RHO deposition
            figure('Position',[1 -100 1024 768],'Color','white');
            cmap = colormap('lines');
            for i=1:beam_data.nbeams
                dex = and(depo_dex',beam_dex(:,i));
                [counts,bins]=hist(sqrt(beam_data.S_lines(2,dex)),0:0.01:1.25);
                hold on;
                plot(bins,counts,'color',cmap(i,:));
                hold off;
            end
            set(gca,'FontSize',24);
            xlabel('Rho [norm]');
            ylabel('# of Particles');
            title('Radial Birth Distribtuion');
            % vll deposition
            figure('Position',[1 -100 1024 768],'Color','white');
            cmap = colormap('lines');
            for i=1:beam_data.nbeams
                dex = and(depo_dex',beam_dex(:,i));
                [counts,bins]=hist(beam_data.vll_lines(3,dex)./1E6,100);
                hold on;
                plot(bins,counts,'color',cmap(i,:));
                hold off;
            end
            set(gca,'FontSize',24);
            xlabel('Parallel Velocity (x10^6) [m/s]]');
            ylabel('# of Particles');
            title('Parallel Velocity Birth Distribution');
            % mu deposition
            figure('Position',[1 -100 1024 768],'Color','white');
            cmap = colormap('lines');
            for i=1:beam_data.nbeams
                dex = and(depo_dex',beam_dex(:,i));
                [counts,bins]=hist(beam_data.moment_lines(3,dex).*1E16,100);
                hold on;
                plot(bins,counts,'color',cmap(i,:));
                hold off;
            end
            set(gca,'FontSize',24);
            xlabel('Magnetic Moment (x10^{16}) [J/T]');
            ylabel('# of Particles');
            title('Magnetic Moment Birth Distribution');
            % E deposition
            figure('Position',[1 -100 1024 768],'Color','white');
            cmap = colormap('lines');
            for i=1:beam_data.nbeams
                dex = and(depo_dex',beam_dex(:,i));
                E_birth=0.5.*(beam_data.mass(dex)'.*beam_data.vll_lines(3,dex).^2+...
                    abs(beam_data.moment_lines(3,dex).*beam_data.B_lines(3,dex)));
                [counts,bins]=hist(1E-3.*E_birth./beam_data.charge(dex)',100);
                hold on;
                plot(bins,counts,'color',cmap(i,:));
                hold off;
            end
            set(gca,'FontSize',24);
            xlabel('Energy [keV]');
            ylabel('# of Particles');
            title('Birth Energy Distribution');
            % Birth Dist Func
            colormap hot;
            for i=1:beam_data.nbeams
                figure('Position',[1 -100 1024 768],'Color','white');
                dex = and(depo_dex',beam_dex(:,i));
                v_perp = sqrt(2.*beam_data.moment_lines(3,dex).*beam_data.B_lines(3,dex)./beam_data.mass(dex)');
                v_ll   = beam_data.vll_lines(3,dex);
                vlmax = max(abs(v_ll));
                vpmax = max(v_perp).*2;
                vmax  = max([vlmax vpmax])./1E6;
                [vals, C]=hist3([v_ll; v_perp]'./1E6,{-vmax:2*vmax/1023:vmax, 0:vmax./1023:vmax});
                hold on;
                ha=pcolor(C{1},C{2},vals');
                set(ha,'EdgeColor','none');
                set(gca,'FontSize',24);
                xlabel('V_{ll} [Mm/s]');
                ylabel('V_{perp} [Mm/s]');
                title(['Birth Distribution BEAM:' num2str(i,'%3.3i')])
                axis tight;
                hold off;
            end
            figure('Position',[1 -100 1024 768],'Color','white');
            dex = depo_dex';
            v_perp = sqrt(2.*beam_data.moment_lines(3,dex).*beam_data.B_lines(3,dex)./beam_data.mass(dex)');
            v_ll   = beam_data.vll_lines(3,dex);
            vlmax = max(abs(v_ll));
            vpmax = max(v_perp).*2;
            vmax  = max([vlmax vpmax])./1E6;
            [vals, C]=hist3([v_ll; v_perp]'./1E6,{-vmax:2*vmax/1023:vmax, 0:vmax./1023:vmax});
            hold on;
            ha=pcolor(C{1},C{2},vals');
            set(ha,'EdgeColor','none');
            set(gca,'FontSize',24);
            xlabel('V_{ll} [Mm/s]');
            ylabel('V_{perp} [Mm/s]');
            title('Birth Distribution Total')
            axis tight;
            hold off;
    end
    hold on;
    end
end
return
end



