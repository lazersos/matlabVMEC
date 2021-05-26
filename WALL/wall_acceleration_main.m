% Main for making an accelerated wall
% Run by setting the filename at the top and giving it a block size
% Basically a helper script for when using wall_acceleration.m

%   Created by: D.J. Engels (d.j.engels@student.tue.nl)
%   Version:    1.0
%   Date:       May 2021

%% initialize
close all; addpath('c');

filename = '/home/dion/Dropbox/__Internship/Internship_local/W7X_FILES/w7x_divertor_op12b_fullres_0_025.dat';
%filename = '/home/dion/WALL_TEST/a10_tokamak_wall_0_25.dat';
%filename = '/home/dion/Dropbox/__Internship/STELLOPT/BENCHMARKS/FIELDLINES_TEST/NCSX_wall_trimesh_0_01.dat';

block_size = 0.1;

wall = read_wall(filename);
%% plot
plot_wall(wall);
pause(0.001);

%% run acceleration
wall_accelerated = wall_acceleration(wall, block_size);

%% plot result
plot_wall_accelerated(wall_accelerated);
pause(0.0001);

%% check output
while 1
    happy = input("Are you happy with this mesh. Y/N: ", 's');
    if happy == "Y"
        fprintf("Nice. Continuing to saving of mesh\n");
        break
    elseif happy == "N"
        fprintf("You can change the settings and try again\n")
        return
    else
        fprintf("Invalid input given\n");
    end 
end

%% Save new mesh
filename_new = strcat(filename(1:end-4), '_acc', num2str(block_size), 'v4', filename(end-3:end));
write_wall_accelerated(wall_accelerated, filename_new);
