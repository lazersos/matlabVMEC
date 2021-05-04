%% initialize
clear all; close all; clc;


filename = '/home/dion/WALL_TEST/a10_tokamak_wall.dat';
block_size = 5;



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
filename_new = strcat(filename(1:end-4), '_acc', str(block_size), filename(end-3:end));
write_wall_accelerated(wall_accelerated, filename_new);
