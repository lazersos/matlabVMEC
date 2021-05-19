% Main for analyzing the wall hit location after a FIELDLINES run
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
clear all; close all; hold off;
dir = split(pwd, '/');
dir = strjoin(dir(1:end-1), '/');
addpath(genpath(dir));

l_plot_div = 0;  % whether or not you want extra plots
l_plot_box = 1;

% error field without trim coils
%filename = '/home/dion/Internship_local/Data/05-18/14_20/fieldlines_w7x_eim+250_new.h5';
% no error field
%filename = '/home/dion/Internship_local/Data/05-17/12_00/fieldlines_w7x_eim+250_new.h5';
% any error field with trim coils
filename = '/home/dion/Dropbox/__Internship/Internship_local/Data/05-18/17_15_100A/144/fieldlines_w7x_eim+250_new.h5';


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
result.divertor_heat_load(splits) = struct();


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
    
    result.divertor_heat_load(i).lower_div = sum(lower_div);
    result.divertor_heat_load(i).upper_div = sum(upper_div);
    
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
tmp.upper_mean = mean([result.divertor_heat_load.upper_div]);
tmp.lower_mean = mean([result.divertor_heat_load.lower_div]);
tmp.upper_std = std([result.divertor_heat_load.upper_div]);
tmp.lower_std = std([result.divertor_heat_load.lower_div]);
tmp.upper_std_rel = tmp.upper_std / tmp.upper_mean;
tmp.lower_std_rel = tmp.lower_std / tmp.lower_mean;

result.divertor_heat_result = tmp;

clear tmp
%% Plot divertor result
hold off;
figure;
bar([result.divertor_heat_load.upper_div]);
names={'Div 1'; 'Div 2'; 'Div 3'; 'Div 4'; 'Div 5';};
set(gca,'xticklabel',names)
hold on;
bar(-[result.divertor_heat_load.lower_div]);
yu = yline(result.divertor_heat_result.upper_mean, '-', 'Mean (Upper)', 'LineWidth',3);
yl = yline(-result.divertor_heat_result.lower_mean, '-', 'Mean (Lower)', 'LineWidth',3);
ym = yline(mean([result.divertor_heat_load.upper_div, result.divertor_heat_load.lower_div]), '-', 'Mean', 'LineWidth',3);
ym2 = yline(-mean([result.divertor_heat_load.upper_div, result.divertor_heat_load.lower_div]), '-', 'Mean', 'LineWidth',3);
ym.LabelHorizontalAlignment = 'left';
ym.Color = [1 0 0];
ym2.LabelHorizontalAlignment = 'left';
ym2.Color = [1 0 0];
title('Particle hits per divertor'); xlabel('Divertor Nr.'); ylabel({'Hits (-)';'Positive: upper. Negative: lower divertor'});
pause(0.001);

clear yl ym ym2 yu names
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
    % First do upper divertor
    phi_min = phi_start + (i-1)*pi/2.5;
    
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
    
    % Then lower
    phi_min = (i-1)*pi/2.5 - phi_start;
    lower_div = angle >= phi_min - phi_size & angle < phi_min ...
        & data.Z_lines(:,2) < -z_min & data.Z_lines(:,2) >= -z_max ...
        & dist >= r_min & dist < r_max;
    
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
    
    % TO DO: find hits in boxes, proper placing of boxes
    
end

clear r_min r_max z_min z_max phi_start phi_size phi_min res phi
clear x_bounds x_inner x_left x_outer x_right
clear y_bounds y_inner y_left y_outer y_right
clear z_bounds z_move
clear upper_div lower_div i
%% Run head load analysis for horizontal divertor
% dex = zeros(3,1);
% cutoff = 0.2;
% height = 0.88;
%
% % find faces that are flat in z
% bool = false(length(data.wall_faces), 1);
% for i=1:length(data.wall_faces)
%     dex = data.wall_faces(i,:);
%
%     zdif = abs(max(data.wall_vertex(dex,3)) - min(data.wall_vertex(dex,3)));
%     ydif = abs(max(data.wall_vertex(dex,2)) - min(data.wall_vertex(dex,2)));
%     xdif = abs(max(data.wall_vertex(dex,1)) - min(data.wall_vertex(dex,1)));
%
%     rel = zdif / sqrt(ydif^2 + xdif^2);
%     bool(i) = rel < cutoff & (max(data.wall_vertex(dex,3)) < -height | min(data.wall_vertex(dex,3)) > height);
% end
% data2 = data;
% data2.wall_faces = data.wall_faces(bool, :);
% data2.wall_strikes = data.wall_strikes(bool);
%
% figure;
% plot_fieldlines(data2, 'wall_strike');
% hold on;
% xlabel('x (m)');ylabel('y (m)');zlabel('z (m)'); title('Hits per face');
%
%
% M = [data2.wall_faces(:,[1,2]) ; data2.wall_faces(:,[2,3]); data2.wall_faces(:,[3,1])];   % Find all edges in mesh, note internal edges are repeated
% [u,~,n] = unique(sort(M,2),'rows'); % determine uniqueness of edges
% counts = accumarray(n(:), 1);   % determine counts for each unique edge
% O = u(counts==1,:); % extract edges that only occurred once
% I = O(:,1:2)';
% [x0,y0,z0] = deal(data2.wall_vertex(I,1),data2.wall_vertex(I,2),data2.wall_vertex(I,3));
%
% % Plot Results
% % plot3(x0,y0,z0, 'o')
% %%
%
% r = sqrt(x0.^2+y0.^2);
%
% r_min = 5.15;
% r_max = 5.65;
% mask =r > r_min & r <r_max;
% plot3(x0(mask), y0(mask), z0(mask), 'o')
%
% figure;hist(r,50);
%
% angle2 = atan2(y0, x0);
% angle2(angle2<-pi/5) = angle2(angle2<-pi/5)+2*pi;
% angle2(angle2>2*pi - pi/5) = angle2(angle2>2*pi - pi/5)-2*pi;
%
% figure;histogram(angle2(z0>0),300); title('Upper divertors');
% figure;histogram(angle2(z0<0),300); title('Lower divertors');

% %%
% start_low = -0.34;
% start_high = -0.12;
% width = 0.3;
% splits=5;
%
% angle2 = atan2(data.wall_vertex(:,2), data.wall_vertex(:,1));
% angle2(angle2<-pi/5) = angle2(angle2<-pi/5)+2*pi;
% angle2(angle2>2*pi - pi/5) = angle2(angle2>2*pi - pi/5)-2*pi;
%
% bool = false(length(data.wall_vertex(:,2)), 1);
% for i=1:splits
%     lb = start_low + (i-1)*2*pi/splits;
%     ub = start_high + (i-1)*2*pi/splits;
%     upper_div = angle2 < ub+width & angle2 >= ub & data.wall_vertex(:,3) >= 0;
%     lower_div = angle2 < lb+width & angle2 >= lb & data.wall_vertex(:,3) < 0;
%
%     bool = bool | upper_div | lower_div;
% end
%
% data3 = data;
% data3.wall_faces = data.wall_faces(bool, :);
% data3.wall_strikes = data.wall_strikes(bool);
%
% figure;
% plot_fieldlines(data3, 'wall_strike');
% hold on;
% xlabel('x (m)');ylabel('y (m)');zlabel('z (m)'); title('Hits per face');