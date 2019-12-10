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
%      plot_beams(beam_data,'lost_len');
%      plot_beams(beam_data,'lost_initial');
%      plot_beams(beam_data,'E-alpha');
%      plot_beams(beam_data,'dist');
%      plot_beams(beam_data,'dist_initial');
%      plot_beams(beam_data,'injection');
%      plot_beams(beam_data,'birth_image');
%      plot_beams(beam_data,'birth_xyz');
%      plot_beams(beam_data,'birth_xyz_b');
%      plot_beams(beam_data,'birth_xyz_s');
%      plot_beams(beam_data,'wall_loss');
%      plot_beams(beam_data,'benchmark');
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.50

ec=1.6021773300E-19; % Charge of an electron (leave alone)
camera=[];
% Handle varargin
plot_type = 'overview';
beamdex = -1;
if nargin > 1
    for i=1:nargin-1
        switch varargin{i}
            case {'overview','lost_len','lost_inital','lost_flux',...
                    'deposition','injection','birth_image','birth_xyz',...
                    'birth_xyz_s','birth_xyz_b','wall_loss','benchmarks',...
                    'grid','grid_s','dist','dist_initial','e-alpha','pitch',...
                    'camview'}
                plot_type=varargin{i};
            case 'beam'
                beamdex = varargin{i+1};
        end
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
    % First determine type of run
    shine_dex = zeros(1,beam_data.nparticles);
    last_dex = zeros(1,beam_data.nparticles)+double(beam_data.npoinc);
    lost_dex = zeros(1,beam_data.nparticles);
    therm_dex = zeros(1,beam_data.nparticles);
    % Set last dex
    % Calculate the last index of each particle
    % Note that the indicies run from 1 to npoinc+1
    % Note that the last point recorded is the wall point
    for i = 1:beam_data.nparticles
        dex = find(beam_data.R_lines(:,i)>0,1,'last');
        last_dex(i) = dex;
    end
    % Handle Injection vs. non-Injection
    if any(beam_data.neut_lines(1,:) > 0) %then NBI run
        % Mark Shinethrough Particles
        % Calculate which particles just shine through
        % 1: Launch Point
        % 2: Last point or wall point
        % 3: Only recorded if particle hits the wall
        shine_dex(beam_data.neut_lines(2,:) > 0) = 2;
        for i=1:beam_data.nparticles
            s = beam_data.S_lines(last_dex(i),i);
            if (s>1.1)
                lost_dex(i) = last_dex(i);
            elseif last_dex(i) < beam_data.npoinc
                therm_dex(i) = last_dex(i);
            end
        end
        dex = shine_dex==2;
        lost_dex(dex) = 0;
        % Handle a depo run
        if beam_data.ldepo
            lost_dex = zeros(1,beam_data.nparticles);
            therm_dex = zeros(1,beam_data.nparticles);
        end
    else % non-nbi run
        for i=1:beam_data.nparticles
            s = beam_data.S_lines(last_dex(i),i);
            if (s>1.1)
                lost_dex(i) = last_dex(i);
            elseif last_dex(i) < beam_data.npoinc
                therm_dex(i) = last_dex(i);
            end
        end
    end
    
% Now handle plots
    switch lower(plot_type)
        case 'overview' % Simple endpoint plot
            x=[]; y=[]; z=[];
            x_lost=[]; y_lost=[]; z_lost=[];
            x_therm=[]; y_therm=[]; z_therm=[];
            x_shine=[]; y_shine=[]; z_shine=[];
            for i=1:beam_data.nparticles
                if last_dex(i) < beam_data.npoinc
                    if shine_dex(i)
                        x_shine = [x_shine beam_data.X_lines(shine_dex(i),i)];
                        y_shine = [y_shine beam_data.Y_lines(shine_dex(i),i)];
                        z_shine = [z_shine beam_data.Z_lines(shine_dex(i),i)];
                    elseif lost_dex(i)
                        x_lost = [x_lost beam_data.X_lines(lost_dex(i),i)];
                        y_lost = [y_lost beam_data.Y_lines(lost_dex(i),i)];
                        z_lost = [z_lost beam_data.Z_lines(lost_dex(i),i)];
                    elseif therm_dex(i)
                        x_therm = [x_therm beam_data.X_lines(last_dex(i),i)];
                        y_therm = [y_therm beam_data.Y_lines(last_dex(i),i)];
                        z_therm = [z_therm beam_data.Z_lines(last_dex(i),i)];
                    else
                        x = [x beam_data.X_lines(last_dex(i),i)];
                        y = [y beam_data.Y_lines(last_dex(i),i)];
                        z = [z beam_data.Z_lines(last_dex(i),i)];
                    end
                end
            end
            hold on;
            leg_text={};
            if ~isempty(x_therm), plot3(x_therm,y_therm,z_therm,'.b'); leg_text=[leg_text; 'Thermalized'];end
            if ~isempty(x_lost), plot3(x_lost,y_lost,z_lost,'.r'); leg_text=[leg_text; 'Lost']; end
            if ~isempty(x_shine), plot3(x_shine,y_shine,z_shine,'.g'); leg_text=[leg_text; 'Shinethrough']; end
            if ~isempty(x), plot3(x,y,z,'.k'); leg_text=[leg_text; 'Other']; end
            legend(leg_text); axis equal; title('Final State');
        case 'lost_len' % Lost Length
            line_len=[]; x=[]; y=[]; z=[];
            dex1 = 1;
            if any(shine_dex > 0), dex1=3; end
            for i=1:beam_data.nparticles
                if last_dex(i) < beam_data.npoinc
                    if lost_dex(i)
                        x = [x beam_data.X_lines(last_dex(i),i)];
                        y = [y beam_data.Y_lines(last_dex(i),i)];
                        z = [z beam_data.Z_lines(last_dex(i),i)];
                        dx = beam_data.X_lines(dex1+1:last_dex(i),i)-beam_data.X_lines(dex1:last_dex(i)-1,i);
                        dy = beam_data.Y_lines(dex1+1:last_dex(i),i)-beam_data.Y_lines(dex1:last_dex(i)-1,i);
                        dz = beam_data.Z_lines(dex1+1:last_dex(i),i)-beam_data.Z_lines(dex1:last_dex(i)-1,i);
                        dl = sqrt(dx.*dx+dy.*dy+dz.*dz);
                        line_len = [line_len sum(dl)];
                    end
                end
            end
            scatter3(x,y,z,1,line_len./1000);
            colormap hot;
            axis equal;
            axis off;
            title('Lost Particles');
            ha = colorbar;
            ylabel(ha,'Connection Length [km]');
        case 'lost_initial'
            line_len=[]; x=[]; y=[]; z=[];
            dex1 = 1;
            if any(shine_dex > 0), dex1=3; end
            for i=1:beam_data.nparticles
                if last_dex(i) < beam_data.npoinc
                    if lost_dex(i)
                        dex = dex - offset;
                        dx = beam_data.X_lines(dex1+1:last_dex(i),i)-beam_data.X_lines(dex1:last_dex(i)-1,i);
                        dy = beam_data.Y_lines(dex1+1:last_dex(i),i)-beam_data.Y_lines(dex1:last_dex(i)-1,i);
                        dz = beam_data.Z_lines(dex1+1:last_dex(i),i)-beam_data.Z_lines(dex1:last_dex(i)-1,i);
                        dl = sqrt(dx.*dx+dy.*dy+dz.*dz);
                        line_len = [line_len sum(dl)];
                        x = [x beam_data.X_lines(dex1,i)];
                        y = [y beam_data.Y_lines(dex1,i)];
                        z = [z beam_data.Z_lines(dex1,i)];
                    end
                end
            end
            scatter3(x,y,z,1,line_len./1000);
            colormap hot;
            axis equal;
            axis off;
            title('Initial Position');
            ha = colorbar;
            ylabel(ha,'Connection Length [km]');
            title('Connection Length');
        case 'pitch'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vll=[]; vperp=[]; mass=[];
            for i=1:beam_data.nparticles
                if and(lost_dex(i) == 0,shine_dex(i)==0)
                    vll = [vll beam_data.vll_lines(last_dex(i),i)]; 
                    mu = beam_data.moment_lines(last_dex(i),i);
                    b  = beam_data.B_lines(last_dex(i),i);
                    vperp = [vperp sqrt(2.*mu.*b./beam_data.mass(i))];
                    mass = [mass beam_data.mass(i)];
                end
            end
            pitch = atan2d(vll,vperp);
            %E     = 0.5.*mass.*(vll.^2+vperp.^2)./ec;
            %x_size=[0 max(E).*1.25];
            y_size=[min(pitch) max(pitch)];
            nres = 128;
            edges= y_size(1):diff(y_size)./nres:y_size(2);
            vals = hist(pitch,edges);
            %colormap jet;
            plot(edges,vals);
            %ha=pcolor(edges{2},edges{1}./norm,vals);
            %set(ha,'LineStyle','none');
            set(gca,'FontSize',18);
            xlabel('Pitch Angle [^o]');
            ylabel('Counts');
            title('Final Distribution');
        case 'e-alpha'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vll=[]; vperp=[]; mass=[];
            for i=1:beam_data.nparticles
                if and(lost_dex(i) == 0,shine_dex(i)==0)
                    vll = [vll beam_data.vll_lines(last_dex(i),i)]; 
                    mu = beam_data.moment_lines(last_dex(i),i);
                    b  = beam_data.B_lines(last_dex(i),i);
                    vperp = [vperp sqrt(2.*mu.*b./beam_data.mass(i))];
                    mass = [mass beam_data.mass(i)];
                end
            end
            pitch = atan2d(vll,vperp);
            E     = 0.5.*mass.*(vll.^2+vperp.^2)./ec;
            x_size=[0 max(E).*1.25];
            y_size=[-180 180];
            nres = 256;
            edges={x_size(1):diff(x_size)./nres:x_size(2) ...
                y_size(1):diff(y_size)./nres:y_size(2)};
            vals = hist3([E; pitch]',edges);
            if max(abs(x_size)) > 1E6
                norm = 1E6;
                units = 'MeV';
            elseif max(abs(x_size)) > 1E3
                norm = 1E3;
                units = 'keV';
            else
                norm = 1;
                units = 'eV';
            end
            colormap jet;
            ha=pcolor(edges{2},edges{1}./norm,vals);
            set(ha,'LineStyle','none');
            set(gca,'FontSize',18);
            xlabel('Pitch Angle [^o]');
            ylabel(['Energy [' units ']']);
            title('Final Distribution');
        case 'dist'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            vll_lost =[]; vperp_lost=[];
            vll_therm =[]; vperp_therm=[];
            vll_other =[]; vperp_other=[];
            for i=1:beam_data.nparticles
                if last_dex(i) <= beam_data.npoinc
                    if lost_dex(i) > 0
                        vll_lost = [vll_lost beam_data.vll_lines(lost_dex(i),i)];
                        vperp_lost  = [vperp_lost  beam_data.moment_lines(lost_dex(i),i).*2.*beam_data.B_lines(lost_dex(i),i)./beam_data.mass(i)];
                    elseif therm_dex(i) > 0
                        vll_therm = [vll_therm beam_data.vll_lines(therm_dex(i),i)];
                        vperp_therm  = [vperp_therm  beam_data.moment_lines(therm_dex(i),i).*2.*beam_data.B_lines(therm_dex(i),i)./beam_data.mass(i)];
                    elseif (shine_dex(i) ==0)
                        vll_other = [vll_other beam_data.vll_lines(last_dex(i),i)];
                        vperp_other  = [vperp_other  beam_data.moment_lines(last_dex(i),i).*2.*beam_data.B_lines(last_dex(i),i)./beam_data.mass(i)];
                    end
                end
            end
            hold on;
            leg_text={};
            if ~isempty(vll_therm), plot(vll_therm./1E3,sqrt(vperp_therm)./1E3,'.b'); leg_text=[leg_text; 'Thermalized'];end
            if ~isempty(vll_lost), plot(vll_lost./1E3,sqrt(vperp_lost)./1E3,'.r'); leg_text=[leg_text; 'Lost'];end
            if ~isempty(vll_other), plot(vll_other./1E3,sqrt(vperp_other)./1E3,'ok'); leg_text=[leg_text; 'Circulating'];end
            vmax = max(ylim);
            vmax = max([vmax max(abs(xlim))]);
            ylim([0 vmax]);
            xlim([-1 1].*max(ylim));
            axis equal;
            axis square;
            legend(leg_text);
            set(gca,'FontSize',18);
            xlabel('Parallel Velocity [km/s]');
            ylabel('Perp. Velocity [km/s]');
            title('Final Particle Distribution');
        case 'dist_initial'
            figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
            dex1 = 1;
            if any(shine_dex > 0), dex1=3; end
            vll = beam_data.vll_lines(dex1,:);
            mu  = beam_data.moment_lines(dex1,:);
            B   = beam_data.B_lines(dex1,:);
            vperp = sqrt(2.*mu.*B./beam_data.mass');
            dex =mu > 0;
            vll = vll(dex);
            vperp = vperp(dex);
            if isfield(beam_data,'Beam')
                beam  = beam_data.Beam(dex);
                leg_text=[];
                for i = 1:max(beam)
                    dex = beam==i;
                    hold on;
                    plot(vll(dex)./1E3,vperp(dex)./1E3,'.');
                    hold off;
                    leg_text=[leg_text; ['Beam #' num2str(i,'%i')]];
                end
            else
                hold on;
                plot(vll(:)./1E3,vperp(:)./1E3,'.');
                hold off;
                leg_text='Particles';
            end
            vmax = max(ylim);
            vmax = max([vmax max(abs(xlim))]);
            ylim([0 vmax]);
            xlim([-1 1].*max(ylim));
            %axis equal;
            %axis square;
            legend(leg_text);
            set(gca,'FontSize',18);
            xlabel('Parallel Velocity [km/s]');
            ylabel('Perp. Velocity [km/s]');
            title('Initial Particle Distribution');
        case 'lost_flux'
            u=[]; v=[]; line_len=[];
            dex1 = 1;
            if any(shine_dex > 0), dex1=3; end
            for i=1:beam_data.nparticles
                if last_dex(i) < beam_data.npoinc
                    if lost_dex(i) > 0
                        u = [u beam_data.U_lines(dex1,i)];
                        v = [v mod(beam_data.PHI_lines(dex1,i),2*pi)];
                        dx = beam_data.X_lines(dex1+1:lost_dex(i),i)-beam_data.X_lines(dex1:lost_dex(i)-1,i);
                        dy = beam_data.Y_lines(dex1+1:lost_dex(i),i)-beam_data.Y_lines(dex1:lost_dex(i)-1,i);
                        dz = beam_data.Z_lines(dex1+1:lost_dex(i),i)-beam_data.Z_lines(dex1:lost_dex(i)-1,i);
                        dl = sqrt(dx.*dx+dy.*dy+dz.*dz);
                        line_len = [line_len sum(dl)];
                    end
                end
            end
            scatter(v,u,1,line_len./1000,'o');
            xlim([0 2*pi]);
            ylim([0 2*pi]);
            axis equal;
            title('Lost Particle Flux Trajectories');
            xlabel('VMEC Poloidal Angle (u) [rad]');
            ylabel('Toroidal Angle (\phi) [rad]');
            ha = colorbar;
            ylabel(ha,'Connection Length [km]');
        case 'deposition'
            slen = size(beam_data.ndot_prof,2);
            s    = 0:1.0/(slen-1):1;
            plot(sqrt(s),beam_data.ndot_prof);
            xlabel('Effective Radius [\rho/a]');
            ylabel(beam_data.ndot_prof_description);
            title('Particle Source ');
            xlim([0,1]);
        case 'injection'
            x=[]; y=[]; z=[];
            x_shine=[]; y_shine=[]; z_shine=[];
            step = round(beam_data.nparticles/1000);
            for i=1:step:beam_data.nparticles
                if shine_dex(i) > 0
                    x_shine = [x_shine; beam_data.X_lines(1,i), beam_data.X_lines(2,i)];
                    y_shine = [y_shine; beam_data.Y_lines(1,i), beam_data.Y_lines(2,i)];
                    z_shine = [z_shine; beam_data.Z_lines(1,i), beam_data.Z_lines(2,i)];
                else
                    x = [x; beam_data.X_lines(1,i), beam_data.X_lines(2,i)];
                    y = [y; beam_data.Y_lines(1,i), beam_data.Y_lines(2,i)];
                    z = [z; beam_data.Z_lines(1,i), beam_data.Z_lines(2,i)];
                end
            end
            leg_text={};
            hold on;
            if ~isempty(x_shine), plot3(x_shine',y_shine',z_shine','r'); leg_text=[leg_text; 'Shinethrough']; end
            if ~isempty(x), plot3(x',y',z','g'); leg_text=[leg_text; 'Injected'];end
        case 'birth_image'
            r_nb = [];
            z_nb = [];
            for i = 1:beam_data.nparticles
                if shine_dex(i) == 0
                    r_nb = [r_nb beam_data.R_lines(3,i)];
                    z_nb = [z_nb beam_data.Z_lines(3,i)];
                end
            end
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
        case 'birth_xyz'
            x = []; y=[]; z=[];
            for i = 1:beam_data.nparticles
                if shine_dex(i) == 0
                    x = [x beam_data.X_lines(3,i)];
                    y = [y beam_data.Y_lines(3,i)];
                    z = [z beam_data.Z_lines(3,i)];
                end
            end
            plot3(x,y,z,'.b');
        case 'birth_xyz_s'
            x = []; y=[]; z=[]; s=[];
            for i = 1:beam_data.nparticles
                if shine_dex(i) == 0
                    x = [x beam_data.X_lines(3,i)];
                    y = [y beam_data.Y_lines(3,i)];
                    z = [z beam_data.Z_lines(3,i)];
                    s = [s beam_data.S_lines(3,i)];
                end
            end
            scatter3(x,y,z,s.*0.0+1,s,'.')
        case 'birth_xyz_b'
            x = []; y=[]; z=[]; b=[];
            for i = 1:beam_data.nparticles
                if shine_dex(i) == 0
                    x = [x beam_data.X_lines(3,i)];
                    y = [y beam_data.Y_lines(3,i)];
                    z = [z beam_data.Z_lines(3,i)];
                    b = [b beam_data.B_lines(3,i)];
                end
            end
            scatter3(x,y,z,b.*0.0+1,b,'.')
        case 'wall_loss'
            output_args{1}=patch('Vertices',beam_data.wall_vertex,'Faces',beam_data.wall_faces,'FaceVertexCData',beam_data.wall_strikes,'LineStyle','none','CDataMapping','scaled','FaceColor','flat');
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



