% For analyzing the wall hit location after a FIELDLINES run
% For a single dataset
% Based on the w7x_divertor_op12b_fullres mesh (and accelerated versions)
% Can check for each divertor wall hits and split them up in radial or
% toroidal direction

% Run by setting the filename at the top (aim to h5 file with results)

% Make sure you have the file in the MATLAB path, this is a requirement of
% the read_fieldlines script

%   Created by: D.J. Engels (d.j.engels@student.tue.nl)
%   Version:    1.0
%   Date:       May 2021

%% initialize
close all;
main_dir = split(pwd, '/');
main_dir = strjoin(main_dir(1:end-1), '/');
addpath(genpath(main_dir));

l_plot_div = 0;  % whether or not you want extra plots
l_plot_box = 1;

% quick analysis
l_quick_save = 1;

% error field without trim coils
%filename = '/home/dion/Internship_local/Data/05-18/14_20/fieldlines_w7x_eim+250_new.h5';
% no error field
%filename = '/home/dion/Internship_local/Data/05-17/12_00/fieldlines_w7x_eim+250_new.h5';

if l_quick_save
    
    l_plot_div = 0;  % whether or not you want extra plots
    l_plot_box = 0;
    
    phase = 150;
    amp = 116.75;
    
    % any error field with trim coils
    filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-19/16_50_%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase);
    
    % set save information
    dir_save = '~/Dropbox/__Internship/Data/05-19 Divertor loads/';
    save_name = sprintf('I%d_phase%d', amp, phase);
end

data = read_fieldlines(filename);
%% plot
plot_fieldlines(data, 'wall_strike');
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)'); title('Hits per face');
pause(0.001);

figure;
plot3(data.X_lines(:,2), data.Y_lines(:,2), data.Z_lines(:,2), 'o');
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)'); title('Hit locations');

%% Check incomplete lines
if sum(data.X_lines(:,3) == 0) > 0
    fprintf("Warning %d incomplete particle traces found\n", sum(data.X_lines(:,3) == 0));
else
    fprintf("Success, no incomplete particle traces!\n");
end

%% Run divertor heat load analysis
splits = 5;

% For each wall hit, get angle and distance
angle = atan2(data.Y_lines(:,2), data.X_lines(:,2));
angle(angle<-pi/splits) = angle(angle<-pi/splits)+2*pi;
angle(angle>2*pi - pi/splits) = angle(angle>2*pi - pi/splits)-2*pi;
dist = sqrt(data.Y_lines(:,2).^2 + data.X_lines(:,2).^2);

result = struct();
result.particle_load(splits) = struct();


if l_plot_div
    figure;
    hold on;
end

bool = false(length(data.Z_lines(:,2)), 1);
for i=1:splits
    lb = -pi/splits + (i-1)*2*pi/splits;
    ub = -pi/splits + i*2*pi/splits;
    upper_div = angle < ub & angle >= lb & data.Z_lines(:,2) >= 0;
    lower_div = angle < ub & angle >= lb & data.Z_lines(:,2) < 0;
    
    bool = bool | upper_div | lower_div;
    
    result.particle_load(i).lower_div = sum(lower_div);
    result.particle_load(i).upper_div = sum(upper_div);
    
    if l_plot_div
        title('Wall hits per divertor');
        plot3(data.X_lines(upper_div,2), data.Y_lines(upper_div,2), data.Z_lines(upper_div,2), 'o');
        plot3(data.X_lines(lower_div,2), data.Y_lines(lower_div,2), data.Z_lines(lower_div,2), 'o');
        xlabel('x (m)');ylabel('y (m)');zlabel('z (m)'); title('Hit locations');
    end
    
end

clear bool i lb splits ub lower_div upper_div
%% Calculate divertor result
tmp = struct();
tmp.upper_mean = mean([result.particle_load.upper_div]);
tmp.lower_mean = mean([result.particle_load.lower_div]);
tmp.upper_std = std([result.particle_load.upper_div]);
tmp.lower_std = std([result.particle_load.lower_div]);
tmp.upper_std_rel = tmp.upper_std / tmp.upper_mean;
tmp.lower_std_rel = tmp.lower_std / tmp.lower_mean;

result.particle_load_summary = tmp;

clear tmp
%% Plot divertor result
hold off;
figure;
bar([result.particle_load.upper_div]);
names={'Div 1'; 'Div 2'; 'Div 3'; 'Div 4'; 'Div 5';};
set(gca,'xticklabel',names)
hold on;
bar(-[result.particle_load.lower_div]);
yline(result.particle_load_summary.upper_mean, '-', 'Mean (Upper)', 'LineWidth',3);
yline(-result.particle_load_summary.lower_mean, '-', 'Mean (Lower)', 'LineWidth',3);
ym = yline(mean([result.particle_load.upper_div, result.particle_load.lower_div]), '-', 'Mean', 'LineWidth',3);
ym2 = yline(-mean([result.particle_load.upper_div, result.particle_load.lower_div]), '-', 'Mean', 'LineWidth',3);
ym.LabelHorizontalAlignment = 'left';
ym.Color = [1 0 0];
ym2.LabelHorizontalAlignment = 'left';
ym2.Color = [1 0 0];
title('Particle hits per divertor'); xlabel('Divertor Nr.'); ylabel({'Hits (-)';'Positive: upper. Negative: lower divertor'});
pause(0.001);

clear ym ym2 names
%% Analysis distance from pumping gap
% Fill in variables for upper divertor around 0 degrees, rest will be
% mirrored automatically
phi_start = 13;   % in degrees
phi_size = 2;  % in degrees
r_min = 4.95;
r_max = 5.4;
z_min = 0.9;
z_max = 1.05;

res = 100;

%% Run analysis distance from pumping gap
phi_start = phi_start/180*pi;   % in rad
phi_size = phi_size/180*pi;  % in rad

if l_plot_box
    %     figure;
    %     plot_fieldlines(data, 'wall_strike');
    %     xlabel('x (m)');ylabel('y (m)');zlabel('z (m)'); title('Hits per face');
    figure;plot3(data.X_lines(:,2), data.Y_lines(:,2), data.Z_lines(:,2), 'o'); xlabel('x (m)');ylabel('y (m)');zlabel('z (m)')
    hold on;
end

% For each divertor
for i=1:5
    % First make boolean of hits in upper divertor
    phi_min = phi_start + (i-1)*pi/2.5;
    
    % plot box if wanted
    if l_plot_box
        phi = linspace(phi_min, phi_min+phi_size, res);
        [x_inner,y_inner] = pol2cart(phi,r_min.*ones(1,res));
        [x_outer,y_outer] = pol2cart(phi,r_max.*ones(1,res));
        
        [x_left,y_left] = pol2cart(phi_min.*ones(1,res),linspace(r_min, r_max, res));
        [x_right,y_right] = pol2cart((phi_min+phi_size).*ones(1,res),linspace(r_min, r_max, res));
        z_move = linspace(z_min, z_max, res);
        
        x_bounds = [x_inner(1).*ones(1,res), x_inner, x_right(1).*ones(1,res), x_right, fliplr(x_outer), fliplr(x_left), x_inner, x_right, x_right(end).*ones(1,res), fliplr(x_outer), x_outer(1).*ones(1,res), fliplr(x_left), x_inner, x_right];
        y_bounds = [y_inner(1).*ones(1,res), y_inner, y_right(1).*ones(1,res), y_right, fliplr(y_outer), fliplr(y_left), y_inner, y_right, y_right(end).*ones(1,res), fliplr(y_outer), y_outer(1).*ones(1,res), fliplr(y_left), y_inner, y_right];
        z_bounds = [z_move, z_max.*ones(1,res*1), z_move, z_max.*ones(1,res*5), fliplr(z_move), z_min.*ones(1,res), fliplr(z_move), z_min.*ones(1,res*3)];
        plot3(x_bounds, y_bounds, z_bounds);
    end
    
    upper_div = angle < phi_min + phi_size & angle >= phi_min ...
        & data.Z_lines(:,2) >= z_min & data.Z_lines(:,2) < z_max ...
        & dist >= r_min & dist < r_max;
    
    % Then boolean of hits in lower divertor
    phi_min = (i-1)*pi/2.5 - phi_start;
    lower_div = angle >= phi_min - phi_size & angle < phi_min ...
        & data.Z_lines(:,2) < -z_min & data.Z_lines(:,2) >= -z_max ...
        & dist >= r_min & dist < r_max;
    
    % Plot if wanted
    if l_plot_box
        phi = linspace(phi_min, phi_min-phi_size, res);
        [x_inner,y_inner] = pol2cart(phi,r_min.*ones(1,res));
        [x_outer,y_outer] = pol2cart(phi,r_max.*ones(1,res));
        
        [x_left,y_left] = pol2cart(phi_min.*ones(1,res),linspace(r_min, r_max, res));
        [x_right,y_right] = pol2cart((phi_min-phi_size).*ones(1,res),linspace(r_min, r_max, res));
        z_move = linspace(-z_min, -z_max, res);
        
        x_bounds = [x_inner(1).*ones(1,res), x_inner, x_right(1).*ones(1,res), x_right, fliplr(x_outer), fliplr(x_left), x_inner, x_right, x_right(end).*ones(1,res), fliplr(x_outer), x_outer(1).*ones(1,res), fliplr(x_left), x_inner, x_right];
        y_bounds = [y_inner(1).*ones(1,res), y_inner, y_right(1).*ones(1,res), y_right, fliplr(y_outer), fliplr(y_left), y_inner, y_right, y_right(end).*ones(1,res), fliplr(y_outer), y_outer(1).*ones(1,res), fliplr(y_left), y_inner, y_right];
        z_bounds = [z_move, -z_max.*ones(1,res*1), z_move, -z_max.*ones(1,res*5), fliplr(z_move), -z_min.*ones(1,res), fliplr(z_move), -z_min.*ones(1,res*3)];
        plot3(x_bounds, y_bounds, z_bounds);
    end
    
    % Find hits in box
    upper_hits = [data.X_lines(upper_div, 2), data.Y_lines(upper_div, 2), data.Z_lines(upper_div, 2)];
    lower_hits = [data.X_lines(lower_div, 2), data.Y_lines(lower_div, 2), data.Z_lines(lower_div, 2)];
    
    result.hit_locations(i).upper_div = upper_hits;
    result.hit_locations(i).lower_div = lower_hits;
    
    % Save distances from pumping gap for each divertor
    result.hit_from_gap(i).upper_div = dist(upper_div) - r_min;
    result.hit_from_gap(i).lower_div = dist(lower_div) - r_min;
end

clear r_min r_max z_min z_max phi_start phi_size phi_min res phi
clear x_bounds x_inner x_left x_outer x_right
clear y_bounds y_inner y_left y_outer y_right
clear z_bounds z_move
clear upper_div lower_div i upper_hits lower_hits
%% Create summary hit from gap
tmp = struct();
tmp.upper_mean = mean(cat(1,result.hit_from_gap.upper_div));
tmp.lower_mean = mean(cat(1,result.hit_from_gap.lower_div));
tmp.upper_std = std(cat(1,result.hit_from_gap.upper_div));
tmp.lower_std = std(cat(1,result.hit_from_gap.lower_div));
tmp.upper_std_rel = tmp.upper_std / tmp.upper_mean;
tmp.lower_std_rel = tmp.lower_std / tmp.lower_mean;

result.hit_from_gap_summary = tmp;

clear tmp

%% save result
if l_quick_save
    save(strcat(dir_save, save_name, '.mat'), 'result');
end
