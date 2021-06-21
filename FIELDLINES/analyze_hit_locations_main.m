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

% 000: no error field
% 001: self-applied n=1 error field without trim coils
% 002: Sweep n=1: self-applied n=1 error field with trim coils
% 003: self-applied n=1 error field "optimum" runs
% 011: Sweep n=1: self-applied n=1 & n=2 error field
% 051: Sweep n=2: self-applied n=1 & n=2 error field without any trim coils
% 061: Sweep n=2: self-applied n=1 & n=2 error field with optimum trim coils
% 062: Sweep n=2: self-applied n=1 & n=2 error field with optimum trim coils (optima)
% 071: Sweep n=2: self-applied n=2 error field only
% 101: as-built without trim coils
% 102: as-built with trim coils
% 103: as-built with 35 A 36 deg trim coils + sweep coils sweep
% 104: as-built with 35 A 36 deg trim coils -> 2 optima for sweep coils
data_run = 001;

switch data_run
    case {000, 001, 101}
        phase = 0;  % arbitrary, won't be used
        amp = 0; 
    case 104
        phase= [56, 174];
        amp = [118, 88];
    case 003
        phase = [212, 211, 216, 210, 214];
        amp = [111, 113, 113, 115, 115];
    case 062
        phase = [270, 250, 270];
        amp = [71, 91, 91];
    otherwise
        phase = 0:36:350;
        amp = 0;
end

%% Run individual
for i=1:length(phase)
    switch data_run
        case 000
            filename = '/home/dion/Internship_local/Data/05-17/12_00/fieldlines_w7x_eim+250_new.h5';
            dir_save = '';
            save_name = '';
            l_save = 0;  % save results or not
        case 001
            filename = '/home/dion/Internship_local/Data/05-18/14_20/fieldlines_w7x_eim+250_new.h5';
            dir_save = '~/Dropbox/__Internship/Data/05-29 Divertor loads self-applied n=1 v2/';
            save_name = 'I0_phase0';
            l_save = 1;  % save results or not
        case 002
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-28_SelfApplied/n1/%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            dir_save = '~/Dropbox/__Internship/Data/05-29 Divertor loads self-applied n=1 v2/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 003
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-29_SelfApplied/n1_optima/%dA%d/fieldlines_w7x_eim+250_new.h5', amp(i), phase(i));
            dir_save = '~/Dropbox/__Internship/Data/05-29 Divertor loads self-applied n=1 v2/';
            save_name = sprintf('I%d_phase%d', amp(i), phase(i));
            l_save = 1;  % save results or not
        case 011
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-29_SelfApplied/n1_sweep_with_n2/%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            dir_save = '~/Dropbox/__Internship/Data/05-31 Divertor loads self-applied sweep for n=1 with n=2/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 051
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-29_SelfApplied/n2_no_trims/%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            dir_save = '~/Dropbox/__Internship/Data/05-29 Divertor loads self-applied n=1 & n=2 no trims/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 061
            if amp == 75
                filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-31_SelfApplied/%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            else
                filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-29_SelfApplied/n2/%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            end
            dir_save = '~/Dropbox/__Internship/Data/05-29 Divertor loads self-applied n=1 & n=2 optimum trims v1/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 062
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-31_SelfApplied/%dA%d/fieldlines_w7x_eim+250_new.h5', amp(i), phase(i));
            dir_save = '~/Dropbox/__Internship/Data/05-29 Divertor loads self-applied n=1 & n=2 optimum trims v1/';
            save_name = sprintf('I%d_phase%d', amp(i), phase(i));
            l_save = 1;  % save results or not
        case 071
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/06-01_SelfApplied/n2_without_n1/%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            dir_save = '~/Dropbox/__Internship/Data/06-01 Divertor loads self-applied sweep for n=2 with only n=2/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 101
            filename = '/home/dion/Dropbox/__Internship/Internship_local/Data/05-20/14_00_00A/fieldlines_w7x_eim+250_new.h5';
            dir_save = '~/Dropbox/__Internship/Data/05-21 Divertor loads as-built/';
            save_name = 'I0_phase0';
            l_save = 1;  % save results or not
        case 102
            if amp == 30 || amp == 35 || amp == 38
                filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-21/17_00_as_built/%dA_%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            else
                filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-20/Scans/14_05_%dA/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            end
            dir_save = '~/Dropbox/__Internship/Data/05-21 Divertor loads as-built/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 103
            if amp == 0
                filename = '/home/dion/Dropbox/__Internship/Internship_local/Data/05-21/17_00_as_built/35A_36/fieldlines_w7x_eim+250_new.h5';
            else
                filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-26_n2/35A_36_%d/%d/fieldlines_w7x_eim+250_new.h5', amp, phase(i));
            end
            dir_save = '~/Dropbox/__Internship/Data/05-27 Divertor loads as-built sweep coils/';
            save_name = sprintf('I%d_phase%d', amp, phase(i));
            l_save = 1;  % save results or not
        case 104
            filename = sprintf('/home/dion/Dropbox/__Internship/Internship_local/Data/05-27_n2_as_built_opt/%dA%d/fieldlines_w7x_eim+250_new.h5', amp(i), phase(i));
            dir_save = '~/Dropbox/__Internship/Data/05-27 Divertor loads as-built sweep coils/';
            save_name = sprintf('I%d_phase%d', amp(i), phase(i));
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
while 1
    happy = input("Do you want to also run set analysis. Y/N: ", 's');
    if happy == "Y"
        fprintf("Continuing to set analysis.\n");
        break
    elseif happy == "N"
        fprintf("Stopping.\n")
        return
    else
        fprintf("Invalid input given\n");
    end 
end
result_set = analyze_hit_locations_set(dir_save);
