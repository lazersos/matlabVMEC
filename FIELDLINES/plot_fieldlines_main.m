%% initialize
main_dir = split(pwd, '/');
main_dir = strjoin(main_dir(1:end-1), '/');
addpath(genpath(main_dir));

clear main_dir
load('/home/dion/Dropbox/__Internship/Internship_local/Data/06-02_Poincare/AllLoaded.mat')
%% plot
close all;
% settings
n_points = 1;
zoom = 0; % 1=zoom in on 5/5, 0: full plot

for i=1:n_points
    % init
    phase = (i-1)*2*pi/n_points;
    phase_deg = rad2deg(phase);
    
    % base
%     figure;
%     plot_fieldlines(base, 'phi', phase);
%     title({'Ideal coils', sprintf('%d degrees', phase_deg)});   
    
    % n=1
    figure;
    ax = tight_subplot(1,3, [0.05, 0.10], [0.10, 0.10], [0.10, 0.05]);
    axes(ax(1))
    plot_fieldlines(base, 'phi', phase);
    title({'Ideal', sprintf('%d degrees', phase_deg)});   
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(2))
    plot_fieldlines(n1, 'phi', phase);
    title({'n=1 error', sprintf('%d degrees', phase_deg)});     
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(3))
    plot_fieldlines(n1_trim, 'phi', phase);
    title({'n=1 (compensated)', sprintf('%d degrees', phase_deg)});  
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    linkaxes(ax, 'xy');
    
    % n=2
    figure;
    ax = tight_subplot(1,3, [0.05, 0.10], [0.10, 0.10], [0.10, 0.05]);
    axes(ax(1))
    plot_fieldlines(base, 'phi', phase);
    title({'Ideal', sprintf('%d degrees', phase_deg)});   
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(2))
    plot_fieldlines(n2, 'phi', phase);
    title({'n=2 error', sprintf('%d degrees', phase_deg)});     
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(3))
    plot_fieldlines(n2control, 'phi', phase);
    title({'n=2 (compensated)', sprintf('%d degrees', phase_deg)});  
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    linkaxes(ax, 'xy');
    
    % n=1 & n=2
    figure;
    ax = tight_subplot(1,4, [0.05, 0.10], [0.10, 0.15], [0.10, 0.05]);
    axes(ax(1))
    plot_fieldlines(base, 'phi', phase);
    title({'Ideal', sprintf('%d degrees', phase_deg)});   
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(2))
    plot_fieldlines(n1n2, 'phi', phase);
    title({'n=1 & n=2 error', sprintf('%d degrees', phase_deg)});     
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(3))
    plot_fieldlines(n1n2trim, 'phi', phase);
    title({'n=1 & n=2', 'n=1 compensated', sprintf('%d degrees', phase_deg)});  
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(4))
    plot_fieldlines(n1n2trim_control, 'phi', phase);
    title({'n=1 & n=2', 'both compensated', sprintf('%d degrees', phase_deg)});  
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    linkaxes(ax, 'xy');
    
    % as-built
    figure;
    ax = tight_subplot(1,5, [0 0.10], [0.10 0.20], [0.10, 0.05]);
    axes(ax(1))
    plot_fieldlines(base, 'phi', phase);
    title({'Ideal', sprintf('%d degrees', phase_deg)});   
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(2))
    plot_fieldlines(as_built, 'phi', phase);
    title({'As-built', sprintf('%d degrees', phase_deg)});   
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(3))
    plot_fieldlines(as_built_trim, 'phi', phase);
    title({'As-built', 'trim coils', sprintf('%d degrees', phase_deg)});     
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(4))
    plot_fieldlines(as_built_trim_control1, 'phi', phase);
    title({'As-built', 'trim & control coils', '(upper optimum)', sprintf('%d degrees', phase_deg)});  
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    axes(ax(5))
    plot_fieldlines(as_built_trim_control2, 'phi', phase);
    title({'As-built coils', 'trim & control coils', '(lower optimum)', sprintf('%d degrees', phase_deg)});  
    if zoom
        ylabel('z (m)'); xlabel('R (m)'); ylim([-0.2 0.2]); xlim([6.1 6.4]);
    else
        ylabel('z (m)'); xlabel('R (m)'); ylim([0 1.3]); xlim([5 6.5]);
    end
    
    linkaxes(ax, 'xy');
end
