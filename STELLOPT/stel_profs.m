function [ output_args ] = stel_profs(data)
%STEL_PROFS(stel_data) Outputs profile spline knots
%   Detailed explanation goes here

phi = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];
phi = [0 0.1 0.25 0.5 0.75 0.9 0.95 1.0];
npoly = 7;
ltarget=1;
s_edge=1.0;
disp('!-----------------------------------------------------------------------');
disp('!          PROFILE PARAMETERS');
disp('!-----------------------------------------------------------------------');
if isfield(data,'TE_S')
    % Get DATA
    j  = size(data.TE_equil,1);
    te = data.TE_equil(j,:);
    if ltarget; te = data.TE_target(j,:);end
    te_s = data.TE_S(j,:);
    subplot(3,1,1);
    errorbar(te_s,data.TE_target(j,:),data.TE_sigma(j,:),'ok');
    hold on;
    plot(te_s,te,'or');
    xlim([0 1.2]);
    temp = ylim; ylim([0 temp(2)]);
    % Filter Data
    te = te(te_s <= 1.0);
    te_s = te_s(te_s <= 1.0);
    % Sort Data
    [te_s, dex] = sort(te_s);
    te = te(dex);
    te_s = te_s./s_edge;
    % Adjust endpoints
    if min(te_s) > 0
        te_s = [0.0 te_s];
        te   = [te(1) te];
    end
    if max(te_s) < 1
        te_s = [te_s 1.0];
        te   = [te min(te_s)];
    end
    if s_edge < 1.0
        te(te_s > 1.0) = 0.0;
        te_s(te_s > 1.0) = 0.0;
    end
    te_spl = pchip(te_s,te);
    te_pp  = polyfit(te_s,te,npoly);
    te_pp(npoly+1) = max(te);
    te_pp(1)  = -sum(te_pp(2:npoly+1));
    disp('  TE_TYPE = ''akima_spline''');
    disp(['  TE_OPT   = ' sprintf('  %20.10E',fliplr(te_pp))]);
    disp(['  TE_AUX_S = ' sprintf('  %20.10E',phi)]);
    disp(['  TE_AUX_F = ' sprintf('  %20.10E',ppval(te_spl,phi))]);
    plot(phi,ppval(te_spl,phi),'b');
end
if isfield(data,'NE_S')
    % Get DATA
    j  = size(data.NE_equil,1);
    ne = data.NE_equil(j,:);
    if ltarget; ne = data.NE_target(j,:); end
    ne_s = data.NE_S(j,:);
    subplot(3,1,2);
    errorbar(ne_s,data.NE_target(j,:),data.NE_sigma(j,:),'ok');
    hold on;
    plot(ne_s,ne,'or');
    xlim([0 1.2]);
    temp = ylim; ylim([0 temp(2)]);
    % Filter Data
    ne = ne(ne_s <= 1.0);
    ne_s = ne_s(ne_s <= 1.0);
    % Sort Data
    [ne_s, dex] = sort(ne_s);
    ne = ne(dex);
    % Adjust endpoints
    if min(ne_s) > 0
        ne_s = [0.0 ne_s];
        ne   = [ne(1) ne];
    end
    if max(ne_s) < 1
        ne_s = [ne_s 1.0];
        ne   = [ne min(ne_s)];
    end
    % Handle Normalization
    factor = 1.0;
    if (ltarget == 1)
        if (max(ne) > 10.0)
            ne = ne./1.0E18;
        end
    elseif (max(ne) > 1.0)
        ne = ne./1E18;
        factor = 1.0E18;
    else
        if isfield(data,'NELINE_target')
            neline = max(data.NELINE_target(j,:));
            dex    = find(data.NELINE_target(j,:) == neline);
            r0     = data.NELINE_R0(j,dex);
            phi0     = data.NELINE_PHI0(j,dex);
            z0     = data.NELINE_Z0(j,dex);
            r1     = data.NELINE_R1(j,dex);
            phi1     = data.NELINE_PHI1(j,dex);
            z1     = data.NELINE_Z1(j,dex);
            x0     = r0*cos(phi0);   x1     = r1*cos(phi1);
            y0     = r0*sin(phi0);   y1     = r1*sin(phi1);
            dl     = sqrt((x0-x1)^2+(y0-y1)^2+(z0-z1)^2);
            ne_norm = 0.9*neline./dl;
        end
        ne = ne_norm*ne./(max(ne)*1E18);
    end
    ne_spl = pchip(ne_s,ne);
    ne_pp  = polyfit(ne_s,ne,npoly);
    ne_pp(npoly+1) = max(ne);
    ne_pp(1)  = -sum(ne_pp(2:npoly+1));
    disp('  NE_TYPE = ''akima_spline''');
    disp(['  NE_OPT   = ' sprintf('  %20.10E',fliplr(ne_pp))]);
    disp(['  NE_AUX_S = ' sprintf('  %20.10E',phi)]);
    disp(['  NE_AUX_F = ' sprintf('  %20.10E',ppval(ne_spl,phi))]);
    plot(phi,ppval(ne_spl,phi)*factor,'b');
    plot(phi,ppval(ne_spl,phi)./ppval(ne_spl,0.0),'b');
end
if isfield(data,'TI_S')
    % Get DATA
    j  = size(data.TI_equil,1);
    ti = data.TI_equil(j,:);
    if ltarget; ti = data.TI_target(j,:); end;
    ti_s = data.TI_S(j,:);
    subplot(3,1,3);
    errorbar(ti_s,data.TI_target(j,:),data.TI_sigma(j,:),'ok');
    hold on;
    plot(ti_s,ti,'or');
    xlim([0 1.2]);
    temp = ylim; ylim([0 temp(2)]);
    % Filter Data
    ti = ti(ti_s <= 1.0);
    ti_s = ti_s(ti_s <= 1.0);
    % Sort Data
    [ti_s, dex] = sort(ti_s);
    ti= ti(dex);
    % Adjust endpoints
    if min(ti_s) > 0
        ti_s = [0.0 ti_s];
        ti   = [ti(1) ti];
    end
    if max(ti_s) < 1
        ti_s = [ti_s 1.0];
        ti   = [ti min(ti_s)];
    end
    ti_spl = pchip(ti_s,ti);
    ti_pp  = polyfit(ti_s,ti,npoly);
    ti_pp(npoly+1) = max(ti);
    ti_pp(1)  = -sum(ti_pp(2:npoly+1));
    disp('  TI_TYPE = ''akima_spline''');
    disp(['  TI_OPT   = ' sprintf('  %20.10E',fliplr(ti_pp))]);
    disp(['  TI_AUX_S = ' sprintf('  %20.10E',phi)]);
    disp(['  TI_AUX_F = ' sprintf('  %20.10E',ppval(ti_spl,phi))]);
    plot(phi,ppval(ti_spl,phi),'b');
end



end

