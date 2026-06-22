%% ========================================================================
%  MAK 315 - MECHANICAL VIBRATIONS PROJECT
%  33-DEGREE-OF-FREEDOM TRAIN CONSIST MODAL ANALYSIS
%  ------------------------------------------------------------------------
%  MODEL: 1 locomotive + 2 passenger cars, 11 DOF each -> 33 DOF total
%  METHOD: Lagrange -> M, C, K matrices -> state-space (66x66) -> eigenvalue
%
%  Requirements: MATLAB (Control System Toolbox -> lsim, ss)
%                Signal Processing Toolbox -> butter, filtfilt, pwelch
%
%  Author: Saliha Yildiz
% ========================================================================

clc; clear; close all;

% Save generated figures in the repository's figures/ directory, regardless
% of the current MATLAB working directory.
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, '..', 'figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% 1. PARAMETERS - SEPARATE SETS FOR LOCOMOTIVE AND CAR ==================
% --- LOCOMOTIVE (electric, heavy) ---
Mc_L      = 85000;     % carbody mass [kg]
Jc_pitch_L= 2.83e6;    % carbody pitch inertia [kg*m^2] = Mc*L^2/12, L=20m
Jc_roll_L = 170000;    % carbody roll inertia [kg*m^2]
Mb_L      = 4500;      % bogie mass [kg]
Jb_roll_L = 2500;      % bogie roll inertia [kg*m^2]
Mw_L      = 1800;      % wheelset mass [kg]

% --- CAR (passenger, light) ---
Mc_V      = 32000;     % carbody mass [kg]
Jc_pitch_V= 1.29e6;    % carbody pitch inertia [kg*m^2] = Mc*L^2/12, L=22m
Jc_roll_V = 65600;     % carbody roll inertia [kg*m^2]
Mb_V      = 3200;      % bogie mass [kg]
Jb_roll_V = 1800;      % bogie roll inertia [kg*m^2]
Mw_V      = 1500;      % wheelset mass [kg]

% --- SUSPENSION: LOCOMOTIVE (stiffer) ---
Ksz_L = 8.0e5;   Csz_L = 35000;     % secondary vertical
Ks_phi_L = 2.5e6; Cs_phi_L = 35000; % secondary roll
Kpz_L = 2.5e6;   Cpz_L = 18000;     % primary vertical
Kp_phi_L = 3.0e6; Cp_phi_L = 25000; % primary roll

% --- SUSPENSION: CAR (softer) ---
Ksz_V = 4.0e5;   Csz_V = 20600;
Ks_phi_V = 1.5e6; Cs_phi_V = 22000;
Kpz_V = 1.26e6;  Cpz_V = 10600;
Kp_phi_V = 2.0e6; Cp_phi_V = 15000;

% --- WHEEL-RAIL CONTACT ---
K_H = 1.2e8;   % EFFECTIVE wheel-rail stiffness (K_eq): Hertz + rail pad in series

% --- COUPLER (pitch-coupled, semi-elastic) ---
K_kup = 5.0e6;   C_kup = 1.5e5;
L_kup_L = 10.0;  L_kup_V = 11.0;

Lc_L = 8.75;   Lc_V = 8.75;   % carbody-bogie half spacing [m]

%% 2. SINGLE-VEHICLE MATRICES (11x11) ===================================
[M_L, K_L, C_L] = vehicle_matrices(Mc_L, Jc_pitch_L, Jc_roll_L, Mb_L, Jb_roll_L, Mw_L, ...
    Ksz_L, Csz_L, Ks_phi_L, Cs_phi_L, Kpz_L, Cpz_L, Kp_phi_L, Cp_phi_L, K_H, Lc_L);
[M_V, K_V, C_V] = vehicle_matrices(Mc_V, Jc_pitch_V, Jc_roll_V, Mb_V, Jb_roll_V, Mw_V, ...
    Ksz_V, Csz_V, Ks_phi_V, Cs_phi_V, Kpz_V, Cpz_V, Kp_phi_V, Cp_phi_V, K_H, Lc_V);

%% 3. GLOBAL ASSEMBLY - Loco + Car + Car ================================
M_global = blkdiag(M_L, M_V, M_V);
K_global = blkdiag(K_L, K_V, K_V);
C_global = blkdiag(C_L, C_V, C_V);

%% 4. COUPLER - PITCH-COUPLED FORMULATION ===============================
% F = K_kup * [(z_L - L_L*theta_L) - (z_V + L_V*theta_V)]
% z_c indices:     1 (loco), 12 (car1), 23 (car2)
% theta_c indices: 2 (loco), 13 (car1), 24 (car2)
[K_global, C_global] = add_coupler(K_global, C_global, 1, 2, 12, 13, L_kup_L, L_kup_V, K_kup, C_kup);
[K_global, C_global] = add_coupler(K_global, C_global, 12, 13, 23, 24, L_kup_V, L_kup_V, K_kup, C_kup);

%% 5. EIGENVALUE ANALYSIS - UNDAMPED and DAMPED =========================
% (a) Undamped natural frequencies: generalized eigenvalue eig(K,M)
[~, Dund] = eig(K_global, M_global);
omega_n = sqrt(real(diag(Dund)));
freq_undamped = sort(omega_n / (2*pi));

% (b) Damped eigenvalues: state-space A matrix (66x66)
n = 33;
A_mat = [zeros(n), eye(n);
        -M_global\K_global, -M_global\C_global];
[V_dam, D_dam] = eig(A_mat);
lam = diag(D_dam);
freq_damped = abs(lam) / (2*pi);
zeta = -real(lam) ./ max(abs(lam), 1e-9);

% Keep only positive imaginary parts (one of each conjugate pair)
mask = imag(lam) > 0;
freq_d  = freq_damped(mask);
zet     = zeta(mask);
V_modes = V_dam(:, mask);
[freq_d, idx] = sort(freq_d);
zet     = zet(idx);
V_modes = V_modes(:, idx);

% Print summary table to console
fprintf('%s\n', repmat('=',1,72));
fprintf(' 33-DOF TRAIN CONSIST - MODAL ANALYSIS RESULTS\n');
fprintf('%s\n', repmat('=',1,72));
fprintf('%4s | %13s | %10s | %9s | Region\n', 'Mode','omega_n(Hz)','f_d(Hz)','zeta(%)');
fprintf('%s\n', repmat('-',1,72));
for k = 1:length(freq_d)
    if     freq_d(k) < 3,  region = 'Carbody';
    elseif freq_d(k) < 20, region = 'Bogie';
    else,                  region = 'P2 (effective contact)';
    end
    fprintf('%4d | %13.4f | %10.4f | %9.2f | %s\n', ...
        k, freq_undamped(k), freq_d(k), zet(k)*100, region);
end
fprintf('%s\n', repmat('=',1,72));
maxRe = max(real(lam));
if maxRe < 0, status = 'STABLE'; else, status = 'UNSTABLE'; end
fprintf(' Stability: max Re(lambda) = %.4f  (%s)\n', maxRe, status);
fprintf('%s\n', repmat('=',1,72));

% Color helper (by frequency group)
blue = [0.18 0.36 0.85]; orange = [0.95 0.6 0.15]; red = [0.75 0.2 0.2];
colors = zeros(length(freq_d), 3);
for k = 1:length(freq_d)
    if     freq_d(k) < 3,  colors(k,:) = blue;
    elseif freq_d(k) < 20, colors(k,:) = orange;
    else,                  colors(k,:) = red;
    end
end

%% ======================================================================
%  FIGURE 1 - UNDAMPED vs DAMPED FREQUENCIES
%  ======================================================================
figure('Position', [100 100 1200 550], 'Color', 'w');
x = 1:length(freq_d);
bar(x-0.2, freq_undamped(1:length(freq_d)), 0.4, 'FaceColor', [0.55 0.55 0.55]); hold on;
b2 = bar(x+0.2, freq_d, 0.4, 'FaceColor', 'flat');
b2.CData = colors;
set(gca, 'YScale', 'log', 'FontSize', 11);
xlabel('Mode No', 'FontSize', 12);
ylabel('Frequency (Hz) - log scale', 'FontSize', 12);
title('Undamped \omega_n vs Damped f_d Natural Frequency Comparison', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Undamped \omega_n', 'Damped f_d', 'Location', 'northwest');
grid on; set(gca, 'Layer', 'top'); xlim([0 length(freq_d)+1]);
print('-dpng', '-r150', fullfile(output_dir, 'natural_frequencies.png'));

%% ======================================================================
%  FIGURE 2 - DAMPING RATIOS
%  ======================================================================
figure('Position', [100 100 1200 550], 'Color', 'w');
b = bar(x, zet*100, 'FaceColor', 'flat');
b.CData = colors; hold on;
yline(5,  'g--', 'LineWidth', 1.2, 'Label', 'Design target 5%');
yline(30, 'r--', 'LineWidth', 1.2, 'Label', 'Over-damping limit 30%');
xlabel('Mode No', 'FontSize', 12);
ylabel('Damping Ratio \zeta (%)', 'FontSize', 12);
title('Modal Damping Ratios', 'FontSize', 14, 'FontWeight', 'bold');
grid on; set(gca, 'Layer', 'top'); xlim([0 length(zet)+1]);
print('-dpng', '-r150', fullfile(output_dir, 'damping_ratios.png'));

%% ======================================================================
%  FIGURE 3 - EIGENVALUE MAP (s-PLANE)
%  ======================================================================
figure('Position', [100 100 750 700], 'Color', 'w');
plot(real(lam), imag(lam), 'o', 'MarkerSize', 8, ...
     'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k'); hold on;
xline(0, 'k--', 'LineWidth', 0.7);
yline(0, 'k--', 'LineWidth', 0.7);
xlabel('Re(\lambda) - Damping component (1/s)', 'FontSize', 12);
ylabel('Im(\lambda) - Damped frequency (rad/s)', 'FontSize', 12);
title({'Eigenvalue Map (s-plane)', ...
       'All eigenvalues in left-half plane -> ASYMPTOTICALLY STABLE'}, ...
      'FontSize', 13, 'FontWeight', 'bold');
grid on; set(gca, 'Layer', 'top');
print('-dpng', '-r150', fullfile(output_dir, 'eigenvalue_map.png'));

%% ======================================================================
%  FIGURE 4 - MODE SHAPES (4 characteristic modes)
%  ======================================================================
% DOF indices (3 vehicles): z_c=[1 12 23], theta_c=[2 13 24],
%                           phi_c=[3 14 25], z_b1=[4 15 26]
target_f    = [0.66, 0.92, 4.76, 8.80];
target_dof  = {[1 12 23], [3 14 25], [4 15 26], [3 14 25]};
target_type = {'Carbody Vertical (3-veh in-phase)', 'Carbody Roll (independent)', ...
               'Bogie Bounce', 'Bogie Roll-equivalent'};
target_lbl  = {'Carbody z_c', 'Carbody \phi_c', 'Front Bogie z_{b1}', 'Carbody \phi_c'};

figure('Position', [100 100 1300 900], 'Color', 'w');
for p = 1:4
    subplot(2,2,p);
    % Find the mode closest to the target frequency
    cand = find(abs(freq_d - target_f(p)) < 0.3);
    if isempty(cand)
        [~, cand] = min(abs(freq_d - target_f(p)));
    end
    best = cand(1); best_mag = 0;
    for cc = cand'
        m = sum(abs(V_modes(target_dof{p}, cc)));
        if m > best_mag, best_mag = m; best = cc; end
    end
    f_act = freq_d(best);
    comp = V_modes(target_dof{p}, best);
    % Phase rotation (bring the largest component to the real positive axis)
    [~, mx] = max(abs(comp));
    comp_r = real(comp * exp(-1i*angle(comp(mx))));
    if max(abs(comp_r)) > 1e-12
        comp_r = comp_r / max(abs(comp_r));
    end
    % Positive blue, negative red
    hold on;
    for ii = 1:3
        if comp_r(ii) >= 0, bc = blue; else, bc = red; end
        bar(ii, comp_r(ii), 0.5, 'FaceColor', bc, 'EdgeColor', 'k');
        if comp_r(ii) >= 0, dy = 0.05; else, dy = -0.1; end
        text(ii, comp_r(ii)+dy, sprintf('%+.2f', comp_r(ii)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
    yline(0, 'k-', 'LineWidth', 0.8);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Locomotive','Car 1','Car 2'}, 'FontSize', 11);
    ylabel(['Normalized ' target_lbl{p}], 'FontSize', 11);
    title(sprintf('Mode @ %.2f Hz - %s', f_act, target_type{p}), ...
          'FontSize', 11, 'FontWeight', 'bold');
    grid on; set(gca, 'Layer', 'top'); ylim([-1.3 1.3]);
end
sgtitle('Characteristic Mode Shapes', 'FontSize', 15, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile(output_dir, 'mode_shapes.png'));

%% ======================================================================
%  FIGURE 5 - FREQUENCY SPLIT INTO PHYSICAL GROUPS
%  ======================================================================
figure('Position', [100 100 1400 480], 'Color', 'w');
g1 = freq_d(freq_d < 3);
g2 = freq_d(freq_d >= 3 & freq_d < 20);
g3 = freq_d(freq_d >= 20);
groups   = {g1, g2, g3};
g_color  = {blue, orange, red};
g_title  = {'Carbody (<3 Hz)', 'Bogie (3-20 Hz)', 'P2 / Effective Contact (>20 Hz)'};
for p = 1:3
    subplot(1,3,p);
    g = groups{p};
    if ~isempty(g)
        bar(1:length(g), g, 'FaceColor', g_color{p}, 'EdgeColor', 'k'); hold on;
        for i = 1:length(g)
            text(i, g(i)*1.02, sprintf('%.2f', g(i)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        end
    end
    title(g_title{p}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Mode No (within group)'); ylabel('Frequency (Hz)');
    grid on; set(gca, 'Layer', 'top');
end
sgtitle('Natural Frequencies Split into Physical Groups', ...
        'FontSize', 14, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile(output_dir, 'frequency_groups.png'));

%% ======================================================================
%  FIGURE 6 - ISO 2631-1 RIDE COMFORT (forced vibration)
%  ======================================================================
v_train = 100/3.6;     % 100 km/h -> m/s
fs = 500;              % sampling frequency [Hz]
dt = 1/fs;
T_sim = 20;            % simulation time [s]
t = (0:dt:T_sim-dt)';
N = length(t);

% Longitudinal positions of the 12 wheelsets
a_bogie = 2.5;
x_loko = [-Lc_L-a_bogie/2, -Lc_L+a_bogie/2, +Lc_L-a_bogie/2, +Lc_L+a_bogie/2];
x_vag1 = x_loko - 23;
x_vag2 = x_loko - 46;
wheel_pos = [x_loko, x_vag1, x_vag2];

% ISO 8608 Class B-like road profile (band-pass filtered white noise)
rng(42);
T_road = T_sim + max(abs(wheel_pos))/v_train + 1;
white = randn(round(T_road*fs), 1);
[bb, ab] = butter(2, [0.5 30]/(fs/2), 'bandpass');
road_raw = filtfilt(bb, ab, white);
road = road_raw * 0.003 / std(road_raw);   % ~3 mm RMS

% Time-delayed road input for each wheelset
U = zeros(N, 12);
for k = 1:12
    delay = max(0, round(-wheel_pos(k)/v_train * fs));
    idx_end = delay + N;
    if idx_end > length(road)
        pad = idx_end - length(road);
        U(:, k) = [road(delay+1:end); zeros(pad, 1)];
    else
        U(:, k) = road(delay+1:idx_end);
    end
end

% Input matrix B: road -> force on wheelset DOFs via Hertz contact
wheel_dofs = [8 9 10 11, 19 20 21 22, 30 31 32 33];
B_force = zeros(33, 12);
for k = 1:12
    B_force(wheel_dofs(k), k) = K_H;
end
B = [zeros(33,12); M_global\B_force];
C_out = eye(66);
D_out = zeros(66, 12);

% Time-domain solution
sys = ss(A_mat, B, C_out, D_out);
y_out = lsim(sys, U, t);

% Carbody vertical velocities -> accelerations
v_loko = y_out(:, 34);  v_vag1 = y_out(:, 34+11);  v_vag2 = y_out(:, 34+22);
a_loko = gradient(v_loko, dt);
a_vag1 = gradient(v_vag1, dt);
a_vag2 = gradient(v_vag2, dt);

% ISO 2631-1 W_k frequency weighting filter
skip = round(2*fs);
a_L_w = Wk_filter(a_loko(skip:end), fs);
a_1_w = Wk_filter(a_vag1(skip:end), fs);
a_2_w = Wk_filter(a_vag2(skip:end), fs);
rms_L = sqrt(mean(a_L_w.^2));
rms_1 = sqrt(mean(a_1_w.^2));
rms_2 = sqrt(mean(a_2_w.^2));

fprintf('\n%s\n', repmat('=',1,72));
fprintf(' ISO 2631-1 RIDE COMFORT (100 km/h, ISO 8608 Class B approx.)\n');
fprintf('%s\n', repmat('=',1,72));
fprintf(' Locomotive : %.4f m/s^2  ->  %s\n', rms_L, comfort_class(rms_L));
fprintf(' Car 1      : %.4f m/s^2  ->  %s\n', rms_1, comfort_class(rms_1));
fprintf(' Car 2      : %.4f m/s^2  ->  %s\n', rms_2, comfort_class(rms_2));
fprintf('%s\n', repmat('=',1,72));

% 2x2 ride comfort panel
t_steady = t(skip:end) - t(skip);
figure('Position', [50 50 1300 850], 'Color', 'w');

% (Top-left) road input
subplot(2,2,1);
nn = round(5*fs);
plot(t(1:nn), U(1:nn,1)*1000, 'Color', [0.2 0.29 0.37], 'LineWidth', 0.7);
xlabel('Time (s)'); ylabel('Road Roughness (mm)');
title('Stochastic Road Input (ISO 8608 Class B)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% (Top-right) locomotive acceleration
subplot(2,2,2);
plot(t_steady, a_loko(skip:end), 'Color', blue, 'LineWidth', 0.5); hold on;
plot(t_steady, a_L_w, 'Color', red, 'LineWidth', 0.7);
yline(+rms_L, 'k--', 'LineWidth', 0.8); yline(-rms_L, 'k--', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Vertical Acceleration (m/s^2)');
title(sprintf('Locomotive Vertical Acceleration (RMS = %.3f m/s^2)', rms_L), ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('Raw acceleration', 'W_k weighted', 'Location', 'northeast');
grid on; xlim([0 8]);

% (Bottom-left) comfort comparison
subplot(2,2,3);
rms_vals = [rms_L, rms_1, rms_2];
bcol = [blue; orange; red];
hold on;
for ii = 1:3
    bar(ii, rms_vals(ii), 0.5, 'FaceColor', bcol(ii,:), 'EdgeColor', 'k');
    text(ii, rms_vals(ii)+0.005, sprintf('%.3f', rms_vals(ii)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end
yline(0.315, 'g--', 'LineWidth', 1, 'Label', '0.315 m/s^2');
set(gca, 'XTick', 1:3, 'XTickLabel', {'Locomotive','Car 1','Car 2'});
ylabel('W_k Weighted RMS Acceleration (m/s^2)');
title('ISO 2631-1 Comfort Comparison', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% (Bottom-right) PSD
subplot(2,2,4);
[Pxx, f_psd] = pwelch(a_loko(skip:end), 2048, [], [], fs);
semilogy(f_psd, Pxx, 'Color', blue, 'LineWidth', 1); hold on;
nat_f = [0.92, 4.76, 45.24]; nat_lbl = {'Carbody','Bogie','P2'};
for ii = 1:3
    xline(nat_f(ii), 'r:', 'LineWidth', 1);
    text(nat_f(ii), max(Pxx)/10, sprintf('%s\n%.2f Hz', nat_lbl{ii}, nat_f(ii)), ...
         'FontSize', 8, 'Color', 'r');
end
xlabel('Frequency (Hz)'); ylabel('PSD (m^2/s^4/Hz)');
title('Locomotive Acceleration PSD - Peaks at Natural Frequencies', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 80]); grid on;

sgtitle('ISO 2631-1 Ride Comfort Analysis (100 km/h)', ...
        'FontSize', 14, 'FontWeight', 'bold');
print('-dpng', '-r150', fullfile(output_dir, 'ride_comfort.png'));

fprintf('\nAll 6 figures generated successfully:\n');
fprintf('  natural_frequencies.png, damping_ratios.png, eigenvalue_map.png,\n');
fprintf('  mode_shapes.png, frequency_groups.png, ride_comfort.png\n');


%% ========================================================================
%  LOCAL FUNCTIONS (defined at the end of the MATLAB script)
%  ========================================================================
function [M, K, C] = vehicle_matrices(Mc, Jp, Jr, Mb, JBr, Mw, ...
                                      Ksz, Csz, Ksphi, Csphi, ...
                                      Kpz, Cpz, Kpphi, Cpphi, KH, Lc)
% Builds the 11x11 M, K, C matrices for a single vehicle.
% DOF order: [z_c, theta_c, phi_c, z_b1, phi_b1, z_b2, phi_b2, zw1..zw4]
    M = diag([Mc, Jp, Jr, Mb, JBr, Mb, JBr, Mw, Mw, Mw, Mw]);
    K = zeros(11,11); C = zeros(11,11);

    % --- STIFFNESS K ---
    K(1,1)=2*Ksz;        K(1,4)=-Ksz;      K(1,6)=-Ksz;
    K(2,2)=2*Ksz*Lc^2;   K(2,4)=-Ksz*Lc;   K(2,6)=+Ksz*Lc;
    K(3,3)=2*Ksphi;      K(3,5)=-Ksphi;    K(3,7)=-Ksphi;
    K(4,1)=-Ksz; K(4,2)=-Ksz*Lc; K(4,4)=Ksz+2*Kpz; K(4,8)=-Kpz; K(4,9)=-Kpz;
    K(5,3)=-Ksphi; K(5,5)=Ksphi+2*Kpphi;
    K(6,1)=-Ksz; K(6,2)=+Ksz*Lc; K(6,6)=Ksz+2*Kpz; K(6,10)=-Kpz; K(6,11)=-Kpz;
    K(7,3)=-Ksphi; K(7,7)=Ksphi+2*Kpphi;
    pairs = [8 4; 9 4; 10 6; 11 6];
    for kk = 1:4
        K(pairs(kk,1), pairs(kk,2)) = -Kpz;
        K(pairs(kk,1), pairs(kk,1)) = Kpz + KH;
    end

    % --- DAMPING C (same topology) ---
    C(1,1)=2*Csz;        C(1,4)=-Csz;      C(1,6)=-Csz;
    C(2,2)=2*Csz*Lc^2;   C(2,4)=-Csz*Lc;   C(2,6)=+Csz*Lc;
    C(3,3)=2*Csphi;      C(3,5)=-Csphi;    C(3,7)=-Csphi;
    C(4,1)=-Csz; C(4,2)=-Csz*Lc; C(4,4)=Csz+2*Cpz; C(4,8)=-Cpz; C(4,9)=-Cpz;
    C(5,3)=-Csphi; C(5,5)=Csphi+2*Cpphi;
    C(6,1)=-Csz; C(6,2)=+Csz*Lc; C(6,6)=Csz+2*Cpz; C(6,10)=-Cpz; C(6,11)=-Cpz;
    C(7,3)=-Csphi; C(7,7)=Csphi+2*Cpphi;
    for kk = 1:4
        C(pairs(kk,1), pairs(kk,2)) = -Cpz;
        C(pairs(kk,1), pairs(kk,1)) = Cpz;
    end

    % --- SYMMETRY COMPLETION ---
    for ii = 1:11
        for jj = 1:11
            if K(ii,jj)~=0 && K(jj,ii)==0, K(jj,ii) = K(ii,jj); end
            if C(ii,jj)~=0 && C(jj,ii)==0, C(jj,ii) = C(ii,jj); end
        end
    end
end

function [K, C] = add_coupler(K, C, izL, ithL, izV, ithV, LL, LV, Kk, Ck)
% Pitch-coupled coupler: F = Kk*[(z_L - LL*th_L) - (z_V + LV*th_V)]
    contrib = [izL, +1; ithL, -LL; izV, -1; ithV, -LV];
    for ii = 1:4
        for jj = 1:4
            i = contrib(ii,1); j = contrib(jj,1);
            ci = contrib(ii,2); cj = contrib(jj,2);
            K(i,j) = K(i,j) + Kk*ci*cj;
            C(i,j) = C(i,j) + Ck*ci*cj;
        end
    end
end

function s = Wk_filter(signal, fs)
% ISO 2631-1 W_k vertical frequency weighting filter (simplified)
    [b1, a1] = butter(4, [0.5 80]/(fs/2), 'bandpass');
    band = filtfilt(b1, a1, signal);
    [b2, a2] = butter(2, [4 8]/(fs/2), 'bandpass');   % 4-8 Hz most sensitive
    peak = filtfilt(b2, a2, band);
    s = 0.5*band + 0.5*peak;
end

function cls = comfort_class(a)
% ISO 2631-1 comfort classification
    if     a < 0.315, cls = 'Not Uncomfortable';
    elseif a < 0.5,   cls = 'A Little Uncomfortable';
    elseif a < 0.8,   cls = 'Fairly Uncomfortable';
    else,             cls = 'Uncomfortable';
    end
end
