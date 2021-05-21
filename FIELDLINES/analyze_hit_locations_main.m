% Main script used to call the single and set functions
% Used by Dion Engels for quick data analysis
% For a single dataset
% WARNING
% Assumed that you used the w7x_divertor_op12b_fullres mesh (or accelerated versions)

%   Created by: D.J. Engels (d.j.engels@student.tue.nl)
%   Version:    1.0
%   Date:       May 2021

%% initialize
close all;
main_dir = split(pwd, '/');
main_dir = strjoin(main_dir(1:end-1), '/');
addpath(genpath(main_dir));

clear main_dir
%% Settings

l_more_plot = 0;  % Want extra plots or not
l_run_set = 1; % Run set analysis as well or not

% 0: no error field
% 1: self-applied error field without trim coils
% 2: self-applied error field with trim coils
% 3: self-applied error field "optimum" run
% 4: as-built without trim coils
% 5: as-built with trim coils
data_run = 3;

phase = 0;
amp = 200;

%% Run individual
for i=1:length(phase)
    switch data_run
        case 0
            % no error field
            filename = '/home/dion/Internship_local/Data/05-17/12_00/fieldlines_w7x_eim+250_new.h5';
            dir_save = '';
            save_name = '';
            l_save = 0;  % save results or not
        case 1
            % self-applied error field without trim coils
            filename = '/home/dion/Internship_local/Data/05-18/14_20/fieldlines_w7x_eim+250_new.h5';
            dir_save = '';
            save_name = '';
            l_save = 0;  % save results or not
        case 2
            % Any self-applied error field with trim coils
            if amp == 100
                filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-18/17_15_%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            else
                filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-19/10_00_%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            end
            dir_save = '~/Dropbox/__Internship/Data/05-19 Divertor loads self-made error field/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 3
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-19/16_50_116.75A/150/fieldlines_w7x_eim+250_new.h5');
            dir_save = '~/Dropbox/__Internship/Data/05-19 Divertor loads self-made error field/';
            save_name = sprintf('I116.75_phase150');
        case 4
            dir_save = '';
            save_name = '';
            l_save = 0;  % save results or not
        case 5
            l_save = 1;  % save results or not
        case default
            error("Unknown data_run set")
    end
    
    result = analyze_hit_locations_single(filename, l_more_plot, l_save, dir_save, save_name);
    
    if ~l_save  % break out if not saved run, then the phases don't matter
        break
    end
end
clear Icoil i I_more_plot I_save filename data_run amp save_name
%% Analyze set
if l_run_set
    result_set = analyze_hit_locations_set(dir_save);
end