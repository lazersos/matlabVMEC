% For analyzing the wall hit location after a FIELDLINES run
% Can be run after doing the "single" variant for many dataset
% This one takes those results and makes compass plots and such

% Run by aiming the script at a folder with your results from
% analyze_hit_locations_single

%   Created by: D.J. Engels (d.j.engels@student.tue.nl)
%   Version:    1.0
%   Date:       May 2021

%% initialize
clear dir; close all;
main_dir = split(pwd, '/');
main_dir = strjoin(main_dir(1:end-1), '/');
addpath(genpath(main_dir));

% aim this at your folder with results
data_dir = '~/Dropbox/__Internship/Data/05-19 Divertor loads/';

%% load in all datasets
files = dir(data_dir); files = files(~[files.isdir]);
files = {files.name};

results = struct();
for i=1:length(files)
    % find phase and amp from name
    file_parts = split(files{i}, '_');
    amp = str2double(file_parts{1}(2:end));
    phase = split(file_parts{2}, '.'); phase = strcat(phase(1:end-1));
    phase = str2double(phase{1}(6:end));
    filename = strcat(data_dir, files{i}, '/');
    res = load(filename(1:end-1));
    results(i).name = files{i};
    results(i).results = res.result;
    results(i).phase = phase;
    results(i).I = amp;
end

clear res filename phase amp file_parts files i

%% analyze datasets 
r1 = analyze_variable(results, 'particle_load_summary', 'particle load');
pause(0.001)
r2 = analyze_variable(results, 'hit_from_gap_summary', 'distance from pumping gap');
% merge results into one struct
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
results_new.lower_div = mergestructs(r1.lower_div,r2.lower_div);
results_new.upper_div = mergestructs(r1.upper_div,r2.upper_div);

clear r1 r2 mergestructs

function r = analyze_variable(results, var_name, res_name)
name_I_opt = sprintf("%s_I_opt", var_name);
name_I = sprintf("%s_I", var_name);
name_full = sprintf("%s_full", var_name);
%% analyze datasets
options = optimoptions('lsqcurvefit', 'MaxIterations', 1000);
lb = [-2000 -2000 0 0 0];
ub = [2000 2000 1 1 1];
phase = [results.phase];
I = [results.I];
I_complex = complex(cosd(phase).*I, sind(phase).*I);
std_upper = [results.results]; std_upper = [std_upper.(var_name)];
std_lower = [std_upper.lower_std_rel]; std_upper = [std_upper.upper_std_rel];


meshI = 0:1:400; meshPhase = 0:1:360;
[meshI2, meshPhase2] = meshgrid(meshI,meshPhase);
meshIcomplex = complex(cosd(meshPhase2).*meshI2, sind(meshPhase2).*meshI2);

surfit = @(B, I) B(3) * (real(I) - B(1)).^2 + B(4) * (imag(I) - B(2)).^2 + B(5);
res = lsqcurvefit(surfit, [-80 60 0.01 0.01 0], I_complex', std_lower', lb, ub, options);
res(2) = -res(2);
r.lower_div.(name_full) = res;
r.lower_div.(name_I) = complex(res(1), res(2));
r.lower_div.(name_I_opt) = [abs(r.lower_div.(name_I)), rad2deg(angle(r.lower_div.(name_I)))];

mesh_res_lower = surfit(res, meshIcomplex);

res = lsqcurvefit(surfit, [-80 60 1 1 0], I_complex', std_upper', lb, ub, options);
res(2) = -res(2);
r.upper_div.(name_full) = res;
r.upper_div.(name_I) = complex(res(1), res(2));
r.upper_div.(name_I_opt) = [abs(r.upper_div.(name_I)), rad2deg(angle(r.upper_div.(name_I)))];

mesh_res_upper = surfit(res, meshIcomplex);
clear grid
%% Plot results
n_rticks = 4;

% polar figures
figure;
subplot(2,2,1)
polarscatter(deg2rad([results.phase]), [results.I], 30, std_upper);
title(sprintf('Upper divertor %s asymmetry', res_name));
rticks(max(I)/n_rticks:max(I)/n_rticks:max(I)+max(I)/n_rticks);
rticklabels({sprintf('%d A', max(I)/4*1), sprintf('%d A', max(I)/4*2), sprintf('%d A', max(I)/4*3), sprintf('%d A', max(I))})
h = colorbar; ylabel(h, 'Relative asymmetry');
hold on
polarscatter(deg2rad(r.upper_div.(name_I_opt)(2)), r.upper_div.(name_I_opt)(1), 200, '+', 'k');
hold off

subplot(2,2,2)
polarscatter(deg2rad([results.phase]), [results.I], 30, std_lower);
title(sprintf('Lower divertor %s asymmetry', res_name));
rticks(max(I)/n_rticks:max(I)/n_rticks:max(I)+max(I)/n_rticks);
rticklabels({sprintf('%d A', max(I)/4*1), sprintf('%d A', max(I)/4*2), sprintf('%d A', max(I)/4*3), sprintf('%d A', max(I))})
h = colorbar; ylabel(h, 'Relative asymmetry');
hold on
polarscatter(deg2rad(r.lower_div.(name_I_opt)(2)), r.lower_div.(name_I_opt)(1), 200, '+', 'k');
hold off

%Re/Im figures
subplot(2,2,3)
s = pcolor(real(meshIcomplex), imag(meshIcomplex), mesh_res_upper);
s.EdgeColor = 'none';
title(sprintf('Upper divertor %s asymmetry', res_name));
xlabel('Re(I)');ylabel('Im(I)'); h = colorbar; ylabel(h, 'Relative asymmetry');
xlim([min(real(I_complex))*1.2, max(real(I_complex))*1.2]); ylim([min(imag(I_complex))*1.2, max(imag(I_complex))*1.2])
caxis([0, 1.5*max(std_upper)])
hold on
scatter(real(I_complex), imag(I_complex), 30, std_upper, 'filled', 'MarkerEdgeColor', 'k')
scatter(real(r.upper_div.(name_I)), imag(r.upper_div.(name_I)), 200, '+', 'k');
hold off

subplot(2,2,4)
s = pcolor(real(meshIcomplex), imag(meshIcomplex), mesh_res_lower);
s.EdgeColor = 'none';
title(sprintf('Lower divertor %s asymmetry', res_name));
xlabel('Re(I)');ylabel('Im(I)'); h = colorbar; ylabel(h, 'Relative asymmetry');
xlim([min(real(I_complex))*1.2, max(real(I_complex))*1.2]); ylim([min(imag(I_complex))*1.2, max(imag(I_complex))*1.2])
caxis([0, 1.5*max(std_lower)])
hold on
scatter(real(I_complex), imag(I_complex), 30, std_lower, 'filled', 'MarkerEdgeColor', 'k')
scatter(real(r.lower_div.(name_I)), imag(r.lower_div.(name_I)), 200, '+', 'k');
hold off
end
