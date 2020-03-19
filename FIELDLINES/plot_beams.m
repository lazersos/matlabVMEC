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
%      plot_beams(beam_data,'flux'); % note _opt as in xyz above
%      plot_beams(beam_data,'pitch'); % note _opt as in xyz above
%      plot_beams(beam_data,'dist'); % note _opt as in xyz above
%      plot_beams(beam_data,'orbit_rz'); % Full particle Orbit in RZ
%      plot_beams(beam_data,'distribution'); % 2D vll/vperp distribution
%      plot_beams(beam_data,'heating'); % 1D Heating profiles
%      plot_beams(beam_data,'current'); % 1D Current Density profiles
%      plot_beams(beam_data,'fueling'); % 1D EP fuelint profiles
%      plot_beams(beam_data,'injection'); % Make an NBI injection plot
%      plot_beams(beam_data,'birth_image'); % R/Z Image plot
%      plot_beams(beam_data,'wall_loss'); % # of particle lost
%      plot_beams(beam_data,'wall_shine'); % [W/m^3]
%      plot_beams(beam_data,'wall_heat'); % [W/m^3]
%      plot_beams(beam_data,'benchmark'); % Old benchmark plots
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.70

ec=1.6021773300E-19; % Charge of an electron (leave alone)
camera=[];
% Handle varargin
plot_type = 'overview';
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
                    'orbit_rz','orbit_flux',...
                    'distribution','heating','current','fueling',...
                    'injection','birth_image',...
                    'wall_loss','wall_heat','wall_shine','benchmarks',...
                    'grid','grid_s',...
                    'camview'}
                plot_type=varargin{i};
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
        last_dex(i) =find(beam_data.R_lines(:,i)>0,1,'last');
    end
    % Make plots
    switch lower(plot_type)
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
        case 'xyz_thermalized'
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
            x1  = beam_data.U_lines(:,orbit_dex);
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
            x1  = beam_data.U_lines(:,therm_dex);
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
            x  = beam_data.U_lines(dex1,born_dex);
            y  = beam_data.S_lines(dex1,born_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Born Particles');
        case 'flux_initial'
            x  = beam_data.U_lines(dex1,orbit_dex);
            y  = beam_data.S_lines(dex1,orbit_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Intial Orbiting Particles');
        case 'flux_lost_initial'
            x  = beam_data.U_lines(dex1,lost_dex);
            y  = beam_data.S_lines(dex1,lost_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Intial Lost Particles');
        case 'flux_therm_initial'
            x  = beam_data.U_lines(dex1,therm_dex);
            y  = beam_data.S_lines(dex1,therm_dex);
            polarplot(x,y,'.k');
            rlim([0 1.5]);
            title('Intial Thermalized Particles');
        case 'pitch'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            vperp = vperp(:,orbit_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(orbit_dex);
            p2    = length(ld);
            for i=1:length(ld)
                p2(i) = pitch(ld(i),i);
            end
            y_size=[min(p2) max(p2)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(p2,edges);
            plot(edges,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Orbiting Distribution');
        case 'pitch_therm'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            vperp = vperp(:,therm_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(therm_dex);
            p2    = length(ld);
            for i=1:length(ld)
                p2(i) = pitch(ld(i),i);
            end
            y_size=[min(p2) max(p2)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(p2,edges);
            plot(edges,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Thermalized Particle Pitch');
        case 'pitch_lost'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            vperp = vperp(:,lost_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(lost_dex);
            p2    = length(ld);
            for i=1:length(ld)
                p2(i) = pitch(ld(i),i);
            end
            y_size=[min(p2) max(p2)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(p2,edges);
            plot(edges,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Lost Particle Pitch');
        case 'pitch_birth'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(dex1,born_dex);
            vperp = vperp(dex1,born_dex);
            pitch = atan2d(vll,vperp);
            p2    = pitch(1,:);
            y_size=[min(p2) max(p2)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(p2,edges);
            plot(edges,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Birth Distribution');
        case 'pitch_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            vperp = vperp(:,orbit_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(orbit_dex);
            p2    = pitch(1,:);
            y_size=[min(p2) max(p2)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(p2,edges);
            plot(edges,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Initial Orbiting Distribution');
        case 'pitch_therm_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            vperp = vperp(:,therm_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(therm_dex);
            p2    = pitch(1,:);
            y_size=[min(p2) max(p2)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(p2,edges);
            plot(edges,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Initial Thermalized Distribution');
        case 'pitch_lost_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            vperp = vperp(:,lost_dex);
            pitch = atan2d(vll,vperp);
            ld    = last_dex(lost_dex);
            p2    = pitch(1,:);
            y_size=[min(p2) max(p2)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(p2,edges);
            plot(edges,vals,'k','LineWidth',4);
            set(gca,'FontSize',24);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Initial Lost Distribution');
        case 'dist'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            vperp = vperp(:,orbit_dex);
            ld    = last_dex(orbit_dex);
            p1    = length(ld);
            p2    = length(ld);
            for i=1:length(ld)
                p1(i) = vll(ld(i),i);
                p2(i) = vperp(ld(i),i);
            end
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([p1; p2]',edges);
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
            ha=pcolor(edges{1}./norm,edges{2}./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Oribiting Distribution');
        case 'dist_therm'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            vperp = vperp(:,therm_dex);
            ld    = last_dex(therm_dex);
            p1    = length(ld);
            p2    = length(ld);
            for i=1:length(ld)
                p1(i) = vll(ld(i),i);
                p2(i) = vperp(ld(i),i);
            end
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([p1; p2]',edges);
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
            ha=pcolor(edges{1}./norm,edges{2}./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Thermalized Distribution');
        case 'dist_lost'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            vperp = vperp(:,lost_dex);
            ld    = last_dex(lost_dex);
            p1    = length(ld);
            p2    = length(ld);
            for i=1:length(ld)
                p1(i) = vll(ld(i),i);
                p2(i) = vperp(ld(i),i);
            end
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([p1; p2]',edges);
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
            ha=pcolor(edges{1}./norm,edges{2}./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Lost Distribution');
        case 'dist_birth'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(dex1,orbit_dex);
            vperp = vperp(dex1,orbit_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([p1; p2]',edges);
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
            ha=pcolor(edges{1}./norm,edges{2}./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Birth Distribution');
        case 'dist_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,orbit_dex);
            vperp = vperp(:,orbit_dex);
            ld    = last_dex(orbit_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([p1; p2]',edges);
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
            ha=pcolor(edges{1}./norm,edges{2}./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Initial Oribiting Distribution');
        case 'dist_therm_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,therm_dex);
            vperp = vperp(:,therm_dex);
            ld    = last_dex(therm_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([p1; p2]',edges);
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
            ha=pcolor(edges{1}./norm,edges{2}./norm,vals');
            set(ha,'LineStyle','none');
            set(gca,'FontSize',24);
            xlabel(['Parallel Velocity (v_{||}) [' units ']']);
            ylabel(['Perpendicuarl Velocity (v_\perp) [' units ']']);
            title('Initial Thermalized Distribution');
        case 'dist_lost_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vperp = beams3d_calc_vperp(beam_data);
            vll = beam_data.vll_lines(:,lost_dex);
            vperp = vperp(:,lost_dex);
            ld    = last_dex(lost_dex);
            p1    = vll(1,:);
            p2    = vperp(1,:);
            x_size=[-1 1].*beam_data.partvmax;
            y_size=[ 0 1].*beam_data.partvmax;
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([p1; p2]',edges);
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
            ha=pcolor(edges{1}./norm,edges{2}./norm,vals');
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
            U = beam_data.U_lines(:,beam_dex);
            polarplot(U,sqrt(S),'.');
            set(gca,'FontSize',24);
            title('Particle Orbits (rho,u)');
        case 'orbit_flux_lost'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            S = beam_data.S_lines(:,lost_dex);
            U = beam_data.U_lines(:,lost_dex);
            polarplot(U,sqrt(S),'.');
            set(gca,'FontSize',24);
            title('Lost Particle Orbits (rho,u)');
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
        case 'wall_loss'
            output_args{1}=patch('Vertices',beam_data.wall_vertex,'Faces',beam_data.wall_faces,'FaceVertexCData',beam_data.wall_strikes,'LineStyle','none','CDataMapping','scaled','FaceColor','flat');
        case 'wall_shine'
            output_args{1}=patch('Vertices',beam_data.wall_vertex,'Faces',beam_data.wall_faces,'FaceVertexCData',sum(beam_data.wall_shine(beamdex,:))','LineStyle','none','CDataMapping','scaled','FaceColor','flat');
            colormap hot;
        case 'wall_heat'
            output_args{1}=patch('Vertices',beam_data.wall_vertex,'Faces',beam_data.wall_faces,'FaceVertexCData',sum(beam_data.wall_load(beamdex,:))','LineStyle','none','CDataMapping','scaled','FaceColor','flat');
            colormap hot;
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
end
return
end



