function [results_new] = analyze_hit_locations_set(data_dir)
% For analyzing the wall hit location after a FIELDLINES run
% Can be run after doing the "single" variant for many dataset
% This one takes those results and makes compass plots and such

% Run by aiming the script at a folder with your results from
% analyze_hit_locations_single

% Inputs:
% - data_dir: directory with all the results of hit_locations_single in it

%   Created by: D.J. Engels (d.j.engels@student.tue.nl)
%   Version:    1.0
%   Date:       May 2021

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
end

function r = analyze_variable(results, var_name, res_name)
name_I_opt = sprintf("%s_I_opt", var_name);
name_I = sprintf("%s_I", var_name);
name_full = sprintf("%s_full", var_name);
name_input = sprintf("%s_input", var_name);
%% analyze datasets
lb = [-2000 -2000 -200 -200 0];
ub = [2000 2000 200 200 0.5];
options = fitoptions('Method', 'NonlinearLeastSquares');
options = fitoptions(options, 'MaxIter',10000,'MaxFunEvals',1E6,...
    'DiffMaxChange',1E-3,'DiffMinChange',1E-14, 'Display', 'iter', 'TolFun', 1E-14, 'TolX', 1E-14, ...
    'StartPoint',[500 500 0 0 0.02], 'Lower',lb,'Upper',ub);
phase = [results.phase];
I = [results.I];
I_complex = complex(cosd(phase).*I, sind(phase).*I);
std_upper = [results.results]; std_upper = [std_upper.(var_name)];
std_lower = [std_upper.lower_std_rel]; std_upper = [std_upper.upper_std_rel];

meshI = 0:1:max(I)*1.8; meshPhase = 0:1:360;
[meshI2, meshPhase2] = meshgrid(meshI,meshPhase);
meshIcomplex = complex(cosd(meshPhase2).*meshI2, sind(meshPhase2).*meshI2);

%Fit
f= @(a,b,c,d,e,x,y) ((x-c)./a).^2 + ((y-d)./b).^2 + e;
res = fit([real(I_complex)' imag(I_complex)'],std_lower',f, options);

r.lower_div.(sprintf("%s_I", name_input)) = I;
r.lower_div.(sprintf("%s_phase", name_input)) = phase;
r.lower_div.(sprintf("%s_std", name_input)) = std_lower;
r.lower_div.(name_full) = res;
r.lower_div.(name_I) = complex(res.c, res.d);
r.lower_div.(name_I_opt) = [abs(r.lower_div.(name_I)), rad2deg(angle(r.lower_div.(name_I)))];

mesh_res_lower = f(res.a, res.b, res.c, res.d, res.e, real(meshIcomplex), imag(meshIcomplex));

res = fit([real(I_complex)' imag(I_complex)'],std_upper',f, options);

r.upper_div.(sprintf("%s_I", name_input)) = I;
r.upper_div.(sprintf("%s_phase", name_input)) = phase;
r.upper_div.(sprintf("%s_std", name_input)) = std_upper;
r.upper_div.(name_full) = res;
r.upper_div.(name_I) = complex(res.c, res.d);
r.upper_div.(name_I_opt) = [abs(r.upper_div.(name_I)), rad2deg(angle(r.upper_div.(name_I)))];

mesh_res_upper = f(res.a, res.b, res.c, res.d, res.e, real(meshIcomplex), imag(meshIcomplex));
clear grid
%% Plot results
n_rticks = 4;

% polar figures
figure;
% subplot(2,2,1)
% % subplot(1,2,1)
% polarscatter(deg2rad([results.phase]), [results.I], 30, std_upper);
% title(sprintf('Upper divertor %s asymmetry', res_name));
% rticks(max(I)/n_rticks:max(I)/n_rticks:max(I)+max(I)/n_rticks);
% rticklabels({sprintf('%d A', max(I)/4*1), sprintf('%d A', max(I)/4*2), sprintf('%d A', max(I)/4*3), sprintf('%d A', max(I))})
% h = colorbar; ylabel(h, 'Relative asymmetry');
% hold on
% polarscatter(deg2rad(r.upper_div.(name_I_opt)(2)), r.upper_div.(name_I_opt)(1), 200, '+', 'k');
% hold off
% 
% subplot(2,2,2)
% % subplot(1,2,2)
% polarscatter(deg2rad([results.phase]), [results.I], 30, std_lower);
% title(sprintf('Lower divertor %s asymmetry', res_name));
% rticks(max(I)/n_rticks:max(I)/n_rticks:max(I)+max(I)/n_rticks);
% rticklabels({sprintf('%d A', max(I)/4*1), sprintf('%d A', max(I)/4*2), sprintf('%d A', max(I)/4*3), sprintf('%d A', max(I))})
% h = colorbar; ylabel(h, 'Relative asymmetry');
% hold on
% polarscatter(deg2rad(r.lower_div.(name_I_opt)(2)), r.lower_div.(name_I_opt)(1), 200, '+', 'k');
% hold off

%Re/Im figures
% subplot(2,2,3)
subplot(1,2,1)
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

% subplot(2,2,4)
subplot(1,2,2)
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