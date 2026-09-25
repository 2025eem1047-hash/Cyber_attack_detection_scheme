
clear;
close all; 
clc;
rng(0);
busdata = bus1();
branchdata = Branch();
nb = size(busdata,1);
nl = size(branchdata,1);
sd_pinj = 0.004;
sd_pflow = 0.008;
[Ybus, branch] = buildYbus(nb, branchdata);
Vm0 = busdata(:,7);
ang0 = deg2rad(busdata(:,8));
P_G = busdata(:,3); 
Q_G = busdata(:,4);
P_d = busdata(:,5); 
Q_d = busdata(:,6);
P_spec = P_G - P_d; 
Q_spec = Q_G - Q_d;
[V_1, ang_true, convpf] = newton_raphson_pf(Ybus, P_spec, Q_spec, busdata, Vm0, ang0);
if ~convpf, warning('PF did not converge; check data'); 
end
Sbus = V_1 .* exp(1j*ang_true) .* conj(Ybus * (V_1 .* exp(1j*ang_true)));
Pinj_true = real(Sbus)
Pflow_true = zeros(nl,1);
for k = 1:nl
    f = branch(k).f; t = branch(k).t;
    Vf = V_1(f)*exp(1j*ang_true(f));
    Vt = V_1(t)*exp(1j*ang_true(t));
    y = branch(k).y; bsh = branch(k).b_sh; tap = branch(k).tap;
    I_ft = (Vf - Vt/tap)*y + 1j*(bsh/2)*Vf;
    Pflow_true(k) = real(Vf * conj(I_ft));
end
z_true = [Pinj_true; Pflow_true]
m = length(z_true);
R = diag([ (sd_pinj^2)*ones(nb,1); (sd_pflow^2)*ones(nl,1) ]);
invR = inv(R);
z = z_true + mvnrnd(zeros(m,1), R)'
[V_est, ang_est, h_final, H_final] = run_wls_estimation(Ybus, branch, z, R, busdata);
res_wls = z - h_final;
mah_wls = res_wls' * (R \ res_wls);
nstate = size(H_final,2);

fprintf('Initial WLS Mah (single sample) = %.6g\n', mah_wls);
fprintf('nstate = %d, m = %d\n', nstate, m);
alpha = 0.0005;
thresh_kf = chi2inv(1 - alpha, m);
fprintf('\nKF threshold formula: thresh_kf = chi2inv(1 - alpha, m)\n');
fprintf(' -> m = %d, alpha = %.5g, thresh_kf = %.6g\n', m, alpha, thresh_kf);
dof_wls = max(1, 0);  
thresh_wls = chi2inv(1 - alpha, dof_wls);
fprintf('\nWLS threshold formula: thresh_wls = chi2inv(1 - alpha, dof_wls)\n');
fprintf(' -> dof_wls = m - nstate = %d, alpha = %.5g, thresh_wls = %.6g\n', dof_wls, alpha, thresh_wls);

fprintf('\nSUMMARY: KF thresh = %.6g (χ²_%d),  WLS thresh = %.6g (χ²_%d)\n', thresh_kf, m, thresh_wls, dof_wls);
c_user = [0 0 0.577 -0.577 0 0 -0.577 0 0 0]';
if length(c_user) ~= nstate
    error('Provided c_user length (%d) != nstate (%d). Adjust c_user or system.', length(c_user), nstate);
end
c0 = c_user / (norm(c_user) + eps);
fprintf('\nUsing given c_user (norm=%.6g)\n', norm(c_user));
a_dir = H_final * c0;
norm_a_dir = norm(a_dir);
fprintf('Measurement-space direction norm = %.6g\n', norm_a_dir);
if norm_a_dir < 1e-12
    warning('H_final*c0 ~ 0 -> provided state direction has negligible measurement effect.');
end

HtRinv = H_final' * invR;
M = HtRinv * H_final + 1e-12 * eye(nstate);
Minv = pinv(M);


mah_wls_orig = mah_wls;

s_vals = logspace(-5,2,800);  
rel_tol_wls = 0.25;            
best = struct('s',[],'mah_kf',-inf,'mah_wls',inf,'rel',inf);
best_max_kf = -inf;


xhat0 = [V_est; ang_est];
P_pred0 = 1e-3 * eye(nstate);
S0 = H_final * P_pred0 * H_final' + R;

fprintf('Searching for s in [%.1e, %.1e] to try cross KF thresh (rel_tol_wls=%.3g)...\n', s_vals(1), s_vals(end), rel_tol_wls);
for i = 1:length(s_vals)
    s_try = s_vals(i);
    a_try = s_try * a_dir;
    z_try = z + a_try;
    xls_try = Minv * (HtRinv * z_try);
    r_wls_try = z_try - H_final * xls_try;
    mah_wls_try = r_wls_try' * (R \ r_wls_try);
    rel_change = abs(mah_wls_try - mah_wls_orig) / (abs(mah_wls_orig) + 1e-12);
    innov_try = z + a_try - H_final * xhat0;
    mah_kf_try = innov_try' * (S0 \ innov_try);

    if rel_change <= rel_tol_wls && mah_kf_try > thresh_kf
        best.s = s_try; best.mah_kf = mah_kf_try; best.mah_wls = mah_wls_try; best.rel = rel_change;
        fprintf(' -> Found s=%.3g : KF Mah=%.6g > thresh_kf (WLS rel=%.3g)\n', s_try, mah_kf_try, rel_change);
        break;
    end
    if rel_change <= rel_tol_wls && mah_kf_try > best_max_kf
        best_max_kf = mah_kf_try;
        best.last = struct('s', s_try, 'mah_kf', mah_kf_try, 'mah_wls', mah_wls_try, 'rel', rel_change);
    end
end

if ~isempty(best.s)
    s_final = best.s;
elseif isfield(best,'last')
    s_final = best.last.s;
    fprintf('No s crossed KF thresh within WLS tol; using best s=%.3g (KF Mah=%.6g, WLS rel=%.3g)\n', s_final, best.last.mah_kf, best.last.rel);
else
    s_final = s_vals(end);
    fprintf('No candidate found; using largest s=%.3g\n', s_final);
end

a_stealth = s_final * a_dir;
fprintf('Final s_final = %.3g, norm(a_stealth) = %.6g\n', s_final, norm(a_stealth));

T = 120;
attack_start = 10;
time = linspace(0,1,T);

attack_idx = min(max(attack_start,1), T);
attack_time_norm = time(attack_idx);
fprintf('Attack starts at index = %d/%d, normalized time = %.6f\n', attack_idx, T, attack_time_norm);

Qkf = 1e-8 * eye(nstate);
P0 = 1e-3 * eye(nstate);

xtrue_state = [V_1; ang_true];


xhat_noatt = zeros(nstate, T);
xhat_att = zeros(nstate, T);
rstat_kf_noatt = zeros(1,T);
rstat_kf_att = zeros(1,T);
rstat_wls_noatt = zeros(1,T);
rstat_wls_att = zeros(1,T);


P = P0; xhat = [V_est; ang_est];
for k = 1:T
    xhat_pred = xhat;
    P_pred = P + Qkf;
    zt = H_final * xtrue_state + mvnrnd(zeros(m,1), R)';
    S_k = H_final * P_pred * H_final' + R;
    innov = zt - H_final * xhat_pred;
    rstat_kf_noatt(k) = innov' * (S_k \ innov);
    Kk = P_pred * H_final' / S_k;
    xhat = xhat_pred + Kk * innov;
    P = (eye(nstate) - Kk * H_final) * P_pred;
    xhat_noatt(:,k) = xhat;
 
    xls_k = Minv * (HtRinv * zt);
    r_wls_k = zt - H_final * xls_k;
    rstat_wls_noatt(k) = r_wls_k' * (R \ r_wls_k);
end

P = P0; xhat = [V_est; ang_est];
for k = 1:T
    xhat_pred = xhat;
    P_pred = P + Qkf;
    zt = H_final * xtrue_state + mvnrnd(zeros(m,1), R)';
    if k >= attack_start
        zt = zt + a_stealth;
    end
    S_k = H_final * P_pred * H_final' + R;
    innov = zt - H_final * xhat_pred;
    rstat_kf_att(k) = innov' * (S_k \ innov);
    Kk = P_pred * H_final' / S_k;
    xhat = xhat_pred + Kk * innov;
    P = (eye(nstate) - Kk * H_final) * P_pred;
    xhat_att(:,k) = xhat;
    % WLS stat
    xls_k = Minv * (HtRinv * zt);
    r_wls_k = zt - H_final * xls_k;
    rstat_wls_att(k) = r_wls_k' * (R \ r_wls_k);
end

% safe automatic scaling detection for voltages (first nb rows of xhat)
maxV_sample = max(max(abs(xhat_noatt(1:nb,:))));
if maxV_sample < 5
    base_kV = 1; % change if your system uses different base
    V_noatt_kV = xhat_noatt(1:nb,:) * base_kV;
    V_att_kV = xhat_att(1:nb,:) * base_kV;
    fprintf('Detected voltage in p.u. (max ~ %.3g). Converting using base_kV = %.3g.\n', maxV_sample, base_kV);
else
    V_noatt_kV = xhat_noatt(1:nb,:);
    V_att_kV = xhat_att(1:nb,:);
    fprintf('Detected voltages already large (max ~ %.3g). No p.u.->kV conversion.\n', maxV_sample);
end

kf_diff = rstat_kf_att - rstat_kf_noatt;

fprintf('\n--- Diagnostics ---\n');
fprintf('KF peak (noatt)=%.6g, KF peak (att)=%.6g\n', max(rstat_kf_noatt), max(rstat_kf_att));
fprintf('WLS peak (noatt)=%.6g, WLS peak (att)=%.6g\n', max(rstat_wls_noatt), max(rstat_wls_att));
exceed_kf_idx = find(rstat_kf_att > thresh_kf);
exceed_wls_idx = find(rstat_wls_att > thresh_wls);
if isempty(exceed_kf_idx)
    fprintf('KF: no exceed of threshold %.6g\n', thresh_kf);
else
    fprintf('KF exceeds at indices: '); fprintf('%d ', exceed_kf_idx); fprintf('\n');
end
if isempty(exceed_wls_idx)
    fprintf('WLS: no exceed of threshold %.6g\n', thresh_wls);
else
    fprintf('WLS exceeds at indices: '); fprintf('%d ', exceed_wls_idx); fprintf('\n');
end
attack_time_norm = attack_time_norm;
figure('Name','Per-bus voltages','Position',[200 50 700 900]);
for b = 1:nb
    ax = subplot(nb,1,b);
    plot(time, V_noatt_kV(b,:), '--', 'LineWidth', 1.4); hold on;
    plot(time, V_att_kV(b,:), '-', 'LineWidth', 1.6);
    xline(attack_time_norm,'k:','LineWidth',1);
    if b==1, legend({'No attack','With attack'}, 'Location','best'); end
    medv = median([V_noatt_kV(b,:), V_att_kV(b,:)]);
    vmin = min(min(V_noatt_kV(b,:)), min(V_att_kV(b,:)));
    vmax = max(max(V_noatt_kV(b,:)), max(V_att_kV(b,:)));
    margin = max(0.05*abs(medv), 0.2);
    ylim([vmin - margin, vmax + margin]);
    ylabel(sprintf('Bus %d [kV]', b));
    if b==nb, xlabel('time (normalized)'); else set(gca,'XTickLabel',[]); end
    grid on; set(gca,'FontSize',10);
end
sgtitle('Per-bus voltages (dashed=no attack, solid=with attack)');

figure('Name','WLS Mah vs time','Position',[950 100 760 300]);
plot(time, rstat_wls_noatt, '--o', 'LineWidth',1.2, 'MarkerSize',4); hold on;
plot(time, rstat_wls_att,   '-s',  'LineWidth',1.4, 'MarkerSize',5);
hline_wls = yline(thresh_wls, 'r--', 'LineWidth',1.6);
hline_wls.Label = sprintf('WLS thresh = %.3g (chi2inv(1-\\alpha, dof=%d))', thresh_wls, dof_wls);
hline_wls.LabelHorizontalAlignment = 'right';
xline(attack_time_norm, 'k:', 'LineWidth',1.2);
xlabel('time'); ylabel('WLS Mah (r'' R^{-1} r)'); title(sprintf('WLS residual vs time (dof=%d, alpha=%.5g)', dof_wls, alpha));
legend({'WLS no attack','WLS with attack','Location','northeastoutside'});
ylim([0, max([max(rstat_wls_noatt), max(rstat_wls_att), thresh_wls])*1.2]);
grid on; set(gca,'FontSize',10);


figure('Name','KF Mah vs time','Position',[950 420 760 300]);
plot(time, rstat_kf_noatt, '--o', 'LineWidth',1.2, 'MarkerSize',4); hold on;
plot(time, rstat_kf_att,   '-s',  'LineWidth',1.4, 'MarkerSize',5);
hline_kf = yline(thresh_kf, 'r--', 'LineWidth',1.6);
hline_kf.Label = sprintf('KF thresh = %.3g (chi2inv(1-\\alpha, m=%d))', thresh_kf, m);
hline_kf.LabelHorizontalAlignment = 'right';
xline(attack_time_norm, 'k:', 'LineWidth',1.2);
xlabel('time'); ylabel('KF Mah (innovation'' S^{-1} innovation)'); title(sprintf('KF innovation Mahalanobis vs time (m=%d, alpha=%.5g)', m, alpha));
legend({'KF no attack','KF with attack','Location','northeastoutside'});
ylim([0, max([max(rstat_kf_noatt), max(rstat_kf_att), thresh_kf])*1.2]);
grid on; set(gca,'FontSize',10);

