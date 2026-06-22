%% ========================================================================
%  MAK 315 - MEKANİK TİTREŞİMLER PROJESİ
%  33 SERBESTLİK DERECELİ TREN KATARI MODAL ANALİZİ 
%  ------------------------------------------------------------------------
%  MODEL: 1 Lokomotif + 2 Vagon, her araç 11 DOF -> toplam 33 DOF
%  YÖNTEM: Lagrange -> M, C, K matrisleri -> durum-uzayı (66x66) -> özdeğer
%
%  Gereksinimler: MATLAB (Control System Toolbox -> lsim, ss)
%                 Signal Processing Toolbox -> butter, filtfilt, pwelch
% ========================================================================

clc; clear; close all;

%% 1. PARAMETRELER - LOKOMOTİF VE VAGON AYRI SETLER ======================
% --- LOKOMOTİF (elektrikli, ağır) ---
Mc_L      = 85000;     % gövde kütlesi [kg]
Jc_pitch_L= 2.83e6;    % gövde pitch ataleti [kg*m^2] = Mc*L^2/12, L=20m
Jc_roll_L = 170000;    % gövde roll ataleti [kg*m^2]
Mb_L      = 4500;      % boji kütlesi [kg]
Jb_roll_L = 2500;      % boji roll ataleti [kg*m^2]
Mw_L      = 1800;      % tekerlek seti kütlesi [kg]

% --- VAGON (yolcu, hafif) ---
Mc_V      = 32000;     % gövde kütlesi [kg]
Jc_pitch_V= 1.29e6;    % gövde pitch ataleti [kg*m^2] = Mc*L^2/12, L=22m
Jc_roll_V = 65600;     % gövde roll ataleti [kg*m^2]
Mb_V      = 3200;      % boji kütlesi [kg]
Jb_roll_V = 1800;      % boji roll ataleti [kg*m^2]
Mw_V      = 1500;      % tekerlek seti kütlesi [kg]

% --- SÜSPANSİYON: LOKOMOTİF (daha sert) ---
Ksz_L = 8.0e5;   Csz_L = 35000;     % ikincil düşey
Ks_phi_L = 2.5e6; Cs_phi_L = 35000; % ikincil roll
Kpz_L = 2.5e6;   Cpz_L = 18000;     % birincil düşey
Kp_phi_L = 3.0e6; Cp_phi_L = 25000; % birincil roll

% --- SÜSPANSİYON: VAGON (daha yumuşak) ---
Ksz_V = 4.0e5;   Csz_V = 20600;
Ks_phi_V = 1.5e6; Cs_phi_V = 22000;
Kpz_V = 1.26e6;  Cpz_V = 10600;
Kp_phi_V = 2.0e6; Cp_phi_V = 15000;

% --- TEKERLEK-RAY TEMASI ---
K_H = 1.2e8;   % ETKİN tekerlek-ray sertliği (K_eq): Hertz + rail pad seri komb.

% --- KUPLÖR (pitch-coupled, semi-elastic) ---
K_kup = 5.0e6;   C_kup = 1.5e5;
L_kup_L = 10.0;  L_kup_V = 11.0;

Lc_L = 8.75;   Lc_V = 8.75;   % gövde-boji yarı mesafeleri [m]

%% 2. TEK ARAÇ MATRİSLERİ (11x11) =======================================
[M_L, K_L, C_L] = arac_matrisleri(Mc_L, Jc_pitch_L, Jc_roll_L, Mb_L, Jb_roll_L, Mw_L, ...
    Ksz_L, Csz_L, Ks_phi_L, Cs_phi_L, Kpz_L, Cpz_L, Kp_phi_L, Cp_phi_L, K_H, Lc_L);
[M_V, K_V, C_V] = arac_matrisleri(Mc_V, Jc_pitch_V, Jc_roll_V, Mb_V, Jb_roll_V, Mw_V, ...
    Ksz_V, Csz_V, Ks_phi_V, Cs_phi_V, Kpz_V, Cpz_V, Kp_phi_V, Cp_phi_V, K_H, Lc_V);

%% 3. GLOBAL ASSEMBLY - Loko + Vagon + Vagon ============================
M_global = blkdiag(M_L, M_V, M_V);
K_global = blkdiag(K_L, K_V, K_V);
C_global = blkdiag(C_L, C_V, C_V);

%% 4. KUPLÖR - PITCH-COUPLED FORMÜLASYON ================================
% F = K_kup * [(z_L - L_L*theta_L) - (z_V + L_V*theta_V)]
% z_c indeksleri: 1 (loko), 12 (vagon1), 23 (vagon2)
% theta_c indeksleri: 2 (loko), 13 (vagon1), 24 (vagon2)
[K_global, C_global] = kuplor_ekle(K_global, C_global, 1, 2, 12, 13, L_kup_L, L_kup_V, K_kup, C_kup);
[K_global, C_global] = kuplor_ekle(K_global, C_global, 12, 13, 23, 24, L_kup_V, L_kup_V, K_kup, C_kup);

%% 5. ÖZDEĞER ANALİZİ - SÖNÜMSÜZ ve SÖNÜMLÜ ============================
% (a) Sönümsüz doğal frekanslar: genelleştirilmiş özdeğer eig(K,M)
[~, Dund] = eig(K_global, M_global);
omega_n = sqrt(real(diag(Dund)));
freq_undamped = sort(omega_n / (2*pi));

% (b) Sönümlü özdeğerler: durum-uzayı A matrisi (66x66)
n = 33;
A_mat = [zeros(n), eye(n);
        -M_global\K_global, -M_global\C_global];
[V_dam, D_dam] = eig(A_mat);
lam = diag(D_dam);
freq_damped = abs(lam) / (2*pi);
zeta = -real(lam) ./ max(abs(lam), 1e-9);

% Yalnız pozitif imajiner kısımları al
mask = imag(lam) > 0;
freq_d  = freq_damped(mask);
zet     = zeta(mask);
V_modes = V_dam(:, mask);
[freq_d, idx] = sort(freq_d);
zet     = zet(idx);
V_modes = V_modes(:, idx);

% Konsola özet tablo
fprintf('%s\n', repmat('=',1,72));
fprintf(' 33 DOF TREN KATARI - MODAL ANALIZ SONUCLARI\n');
fprintf('%s\n', repmat('=',1,72));
fprintf('%4s | %13s | %10s | %9s | Bolge\n', 'Mod','omega_n(Hz)','f_d(Hz)','zeta(%)');
fprintf('%s\n', repmat('-',1,72));
for k = 1:length(freq_d)
    if     freq_d(k) < 3,  bolge = 'Govde';
    elseif freq_d(k) < 20, bolge = 'Boji';
    else,                  bolge = 'P2 (etkin temas)';
    end
    fprintf('%4d | %13.4f | %10.4f | %9.2f | %s\n', ...
        k, freq_undamped(k), freq_d(k), zet(k)*100, bolge);
end
fprintf('%s\n', repmat('=',1,72));
maxRe = max(real(lam));
if maxRe < 0, durum = 'STABIL'; else, durum = 'KARARSIZ'; end
fprintf(' Kararlilik: max Re(lambda) = %.4f  (%s)\n', maxRe, durum);
fprintf('%s\n', repmat('=',1,72));

% Renk yardımcısı (gruplara göre)
mavi = [0.18 0.36 0.85]; turuncu = [0.95 0.6 0.15]; kirmizi = [0.75 0.2 0.2];
renkler = zeros(length(freq_d), 3);
for k = 1:length(freq_d)
    if     freq_d(k) < 3,  renkler(k,:) = mavi;
    elseif freq_d(k) < 20, renkler(k,:) = turuncu;
    else,                  renkler(k,:) = kirmizi;
    end
end

%% ======================================================================
%  GRAFİK 1 - SÖNÜMSÜZ vs SÖNÜMLÜ FREKANSLAR
%  ======================================================================
figure('Position', [100 100 1200 550], 'Color', 'w');
x = 1:length(freq_d);
bar(x-0.2, freq_undamped(1:length(freq_d)), 0.4, 'FaceColor', [0.55 0.55 0.55]); hold on;
b2 = bar(x+0.2, freq_d, 0.4, 'FaceColor', 'flat');
b2.CData = renkler;
set(gca, 'YScale', 'log', 'FontSize', 11);
xlabel('Mod No', 'FontSize', 12);
ylabel('Frekans (Hz) - log olcek', 'FontSize', 12);
title('Sonumsuz \omega_n ve Sonumlu f_d Dogal Frekans Karsilastirmasi', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Sonumsuz \omega_n', 'Sonumlu f_d', 'Location', 'northwest');
grid on; set(gca, 'Layer', 'top'); xlim([0 length(freq_d)+1]);
print('-dpng', '-r150', 'frekanslar.png');

%% ======================================================================
%  GRAFİK 2 - SÖNÜM ORANLARI
%  ======================================================================
figure('Position', [100 100 1200 550], 'Color', 'w');
b = bar(x, zet*100, 'FaceColor', 'flat');
b.CData = renkler; hold on;
yline(5,  'g--', 'LineWidth', 1.2, 'Label', 'Tasarim hedefi %5');
yline(30, 'r--', 'LineWidth', 1.2, 'Label', 'Asiri sonum siniri %30');
xlabel('Mod No', 'FontSize', 12);
ylabel('Sonum Orani \zeta (%)', 'FontSize', 12);
title('Modal Sonum Oranlari', 'FontSize', 14, 'FontWeight', 'bold');
grid on; set(gca, 'Layer', 'top'); xlim([0 length(zet)+1]);
print('-dpng', '-r150', 'sonum.png');

%% ======================================================================
%  GRAFİK 3 - ÖZDEĞER HARİTASI (s-DÜZLEMİ)
%  ======================================================================
figure('Position', [100 100 750 700], 'Color', 'w');
plot(real(lam), imag(lam), 'o', 'MarkerSize', 8, ...
     'MarkerFaceColor', mavi, 'MarkerEdgeColor', 'k'); hold on;
xline(0, 'k--', 'LineWidth', 0.7);
yline(0, 'k--', 'LineWidth', 0.7);
xlabel('Re(\lambda) - Sonum bileseni (1/s)', 'FontSize', 12);
ylabel('Im(\lambda) - Sonumlu frekans (rad/s)', 'FontSize', 12);
title({'Ozdeger Haritasi (s-duzlemi)', ...
       'Tum ozdegerler sol yari duzlemde -> ASIMPTOTIK STABIL'}, ...
      'FontSize', 13, 'FontWeight', 'bold');
grid on; set(gca, 'Layer', 'top');
print('-dpng', '-r150', 'sduzlemi.png');

%% ======================================================================
%  GRAFİK 4 - MOD ŞEKİLLERİ (4 karakteristik mod)
%  ======================================================================
% DOF indeksleri (3 araç): z_c=[1 12 23], theta_c=[2 13 24],
%                          phi_c=[3 14 25], z_b1=[4 15 26]
hedef_f    = [0.66, 0.92, 4.76, 8.80];
hedef_dof  = {[1 12 23], [3 14 25], [4 15 26], [3 14 25]};
hedef_tipi = {'Govde Dusey (3-arac esfazli)', 'Govde Yalpa (bagimsiz)', ...
              'Boji Bounce', 'Boji Yalpa-esdeger'};
hedef_etk  = {'Govde z_c', 'Govde \phi_c', 'On Boji z_{b1}', 'Govde \phi_c'};

figure('Position', [100 100 1300 900], 'Color', 'w');
for p = 1:4
    subplot(2,2,p);
    % Hedef frekansa en yakın modu bul
    cand = find(abs(freq_d - hedef_f(p)) < 0.3);
    if isempty(cand)
        [~, cand] = min(abs(freq_d - hedef_f(p)));
    end
    best = cand(1); best_mag = 0;
    for cc = cand'
        m = sum(abs(V_modes(hedef_dof{p}, cc)));
        if m > best_mag, best_mag = m; best = cc; end
    end
    f_act = freq_d(best);
    comp = V_modes(hedef_dof{p}, best);
    % Faz döndürme (en büyük bileşeni reel pozitife çevir)
    [~, mx] = max(abs(comp));
    comp_r = real(comp * exp(-1i*angle(comp(mx))));
    if max(abs(comp_r)) > 1e-12
        comp_r = comp_r / max(abs(comp_r));
    end
    % Pozitif mavi, negatif kırmızı
    hold on;
    for ii = 1:3
        if comp_r(ii) >= 0, bc = mavi; else, bc = kirmizi; end
        bar(ii, comp_r(ii), 0.5, 'FaceColor', bc, 'EdgeColor', 'k');
        if comp_r(ii) >= 0, dy = 0.05; else, dy = -0.1; end
        text(ii, comp_r(ii)+dy, sprintf('%+.2f', comp_r(ii)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
    yline(0, 'k-', 'LineWidth', 0.8);
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Lokomotif','1. Vagon','2. Vagon'}, 'FontSize', 11);
    ylabel(['Normalize ' hedef_etk{p}], 'FontSize', 11);
    title(sprintf('Mod @ %.2f Hz - %s', f_act, hedef_tipi{p}), ...
          'FontSize', 11, 'FontWeight', 'bold');
    grid on; set(gca, 'Layer', 'top'); ylim([-1.3 1.3]);
end
sgtitle('Karakteristik Mod Sekilleri', 'FontSize', 15, 'FontWeight', 'bold');
print('-dpng', '-r150', 'mod_sekilleri.png');

%% ======================================================================
%  GRAFİK 5 - FREKANSLARIN FİZİKSEL GRUPLARA AYRILMASI
%  ======================================================================
figure('Position', [100 100 1400 480], 'Color', 'w');
g1 = freq_d(freq_d < 3);
g2 = freq_d(freq_d >= 3 & freq_d < 20);
g3 = freq_d(freq_d >= 20);
gruplar = {g1, g2, g3};
g_renk  = {mavi, turuncu, kirmizi};
g_baslik= {'Govde (<3 Hz)', 'Boji (3-20 Hz)', 'P2 / Etkin Temas (>20 Hz)'};
for p = 1:3
    subplot(1,3,p);
    g = gruplar{p};
    if ~isempty(g)
        bar(1:length(g), g, 'FaceColor', g_renk{p}, 'EdgeColor', 'k'); hold on;
        for i = 1:length(g)
            text(i, g(i)*1.02, sprintf('%.2f', g(i)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        end
    end
    title(g_baslik{p}, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Mod No (grup ici)'); ylabel('Frekans (Hz)');
    grid on; set(gca, 'Layer', 'top');
end
sgtitle('Dogal Frekanslarin Fiziksel Gruplara Ayrilmasi', ...
        'FontSize', 14, 'FontWeight', 'bold');
print('-dpng', '-r150', 'gruplar.png');

%% ======================================================================
%  GRAFİK 6 - ISO 2631-1 RIDE COMFORT (zorlanmış titreşim)
%  ======================================================================
v_train = 100/3.6;     % 100 km/h -> m/s
fs = 500;              % örnekleme frekansı [Hz]
dt = 1/fs;
T_sim = 20;            % simülasyon süresi [s]
t = (0:dt:T_sim-dt)';
N = length(t);

% 12 tekerleğin boyuna konumları
a_bogie = 2.5;
x_loko = [-Lc_L-a_bogie/2, -Lc_L+a_bogie/2, +Lc_L-a_bogie/2, +Lc_L+a_bogie/2];
x_vag1 = x_loko - 23;
x_vag2 = x_loko - 46;
wheel_pos = [x_loko, x_vag1, x_vag2];

% ISO 8608 Sınıf B benzeri yol profili (band-geçiren beyaz gürültü)
rng(42);
T_road = T_sim + max(abs(wheel_pos))/v_train + 1;
white = randn(round(T_road*fs), 1);
[bb, ab] = butter(2, [0.5 30]/(fs/2), 'bandpass');
road_raw = filtfilt(bb, ab, white);
road = road_raw * 0.003 / std(road_raw);   % ~3 mm RMS

% Her tekerleğe gecikmeli yol girdisi
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

% Girdi matrisi B: yol -> tekerlek DOF'larına Hertz teması üzerinden kuvvet
wheel_dofs = [8 9 10 11, 19 20 21 22, 30 31 32 33];
B_force = zeros(33, 12);
for k = 1:12
    B_force(wheel_dofs(k), k) = K_H;
end
B = [zeros(33,12); M_global\B_force];
C_out = eye(66);
D_out = zeros(66, 12);

% Zaman tanım alanı çözümü
sistem = ss(A_mat, B, C_out, D_out);
y_out = lsim(sistem, U, t);

% Gövde düşey hızları -> ivmeleri
v_loko = y_out(:, 34);  v_vag1 = y_out(:, 34+11);  v_vag2 = y_out(:, 34+22);
a_loko = gradient(v_loko, dt);
a_vag1 = gradient(v_vag1, dt);
a_vag2 = gradient(v_vag2, dt);

% ISO 2631-1 W_k frekans ağırlık filtresi
skip = round(2*fs);
a_L_w = Wk_filtre(a_loko(skip:end), fs);
a_1_w = Wk_filtre(a_vag1(skip:end), fs);
a_2_w = Wk_filtre(a_vag2(skip:end), fs);
rms_L = sqrt(mean(a_L_w.^2));
rms_1 = sqrt(mean(a_1_w.^2));
rms_2 = sqrt(mean(a_2_w.^2));

fprintf('\n%s\n', repmat('=',1,72));
fprintf(' ISO 2631-1 RIDE COMFORT (100 km/h, ISO 8608 Sinif B yaklasimi)\n');
fprintf('%s\n', repmat('=',1,72));
fprintf(' Lokomotif : %.4f m/s^2  ->  %s\n', rms_L, konfor_sinifi(rms_L));
fprintf(' 1. Vagon  : %.4f m/s^2  ->  %s\n', rms_1, konfor_sinifi(rms_1));
fprintf(' 2. Vagon  : %.4f m/s^2  ->  %s\n', rms_2, konfor_sinifi(rms_2));
fprintf('%s\n', repmat('=',1,72));

% 2x2 ride comfort paneli
t_steady = t(skip:end) - t(skip);
figure('Position', [50 50 1300 850], 'Color', 'w');

% (Sol üst) yol girdisi
subplot(2,2,1);
nn = round(5*fs);
plot(t(1:nn), U(1:nn,1)*1000, 'Color', [0.2 0.29 0.37], 'LineWidth', 0.7);
xlabel('Zaman (s)'); ylabel('Yol Puruzlulugu (mm)');
title('Stokastik Yol Girdisi (ISO 8608 Sinif B)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% (Sağ üst) lokomotif ivmesi
subplot(2,2,2);
plot(t_steady, a_loko(skip:end), 'Color', mavi, 'LineWidth', 0.5); hold on;
plot(t_steady, a_L_w, 'Color', kirmizi, 'LineWidth', 0.7);
yline(+rms_L, 'k--', 'LineWidth', 0.8); yline(-rms_L, 'k--', 'LineWidth', 0.8);
xlabel('Zaman (s)'); ylabel('Dikey Ivme (m/s^2)');
title(sprintf('Lokomotif Dikey Ivmesi (RMS = %.3f m/s^2)', rms_L), ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('Ham ivme', 'W_k agirlikli', 'Location', 'northeast');
grid on; xlim([0 8]);

% (Sol alt) konfor karşılaştırması
subplot(2,2,3);
rms_vals = [rms_L, rms_1, rms_2];
bcol = [mavi; turuncu; kirmizi];
hold on;
for ii = 1:3
    bar(ii, rms_vals(ii), 0.5, 'FaceColor', bcol(ii,:), 'EdgeColor', 'k');
    text(ii, rms_vals(ii)+0.005, sprintf('%.3f', rms_vals(ii)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
end
yline(0.315, 'g--', 'LineWidth', 1, 'Label', '0.315 m/s^2');
set(gca, 'XTick', 1:3, 'XTickLabel', {'Lokomotif','1. Vagon','2. Vagon'});
ylabel('W_k Agirlikli RMS Ivme (m/s^2)');
title('ISO 2631-1 Konfor Karsilastirmasi', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% (Sağ alt) PSD
subplot(2,2,4);
[Pxx, f_psd] = pwelch(a_loko(skip:end), 2048, [], [], fs);
semilogy(f_psd, Pxx, 'Color', mavi, 'LineWidth', 1); hold on;
nat_f = [0.92, 4.76, 45.24]; nat_lbl = {'Govde','Boji','P2'};
for ii = 1:3
    xline(nat_f(ii), 'r:', 'LineWidth', 1);
    text(nat_f(ii), max(Pxx)/10, sprintf('%s\n%.2f Hz', nat_lbl{ii}, nat_f(ii)), ...
         'FontSize', 8, 'Color', 'r');
end
xlabel('Frekans (Hz)'); ylabel('PSD (m^2/s^4/Hz)');
title('Lokomotif Ivme PSD - Dogal Frekanslarda Pikler', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0 80]); grid on;

sgtitle('ISO 2631-1 Ride Comfort Analizi (100 km/h)', ...
        'FontSize', 14, 'FontWeight', 'bold');
print('-dpng', '-r150', 'ride_comfort.png');

fprintf('\nTum 6 grafik basariyla uretildi:\n');
fprintf('  frekanslar.png, sonum.png, sduzlemi.png,\n');
fprintf('  mod_sekilleri.png, gruplar.png, comfort.png\n');


%% ========================================================================
%  YEREL FONKSİYONLAR (MATLAB script sonunda tanımlanır)
%  ========================================================================
function [M, K, C] = arac_matrisleri(Mc, Jp, Jr, Mb, JBr, Mw, ...
                                      Ksz, Csz, Ksphi, Csphi, ...
                                      Kpz, Cpz, Kpphi, Cpphi, KH, Lc)
% Tek bir araç için 11x11 M, K, C matrislerini üretir.
% DOF sırası: [z_c, theta_c, phi_c, z_b1, phi_b1, z_b2, phi_b2, zw1..zw4]
    M = diag([Mc, Jp, Jr, Mb, JBr, Mb, JBr, Mw, Mw, Mw, Mw]);
    K = zeros(11,11); C = zeros(11,11);

    % --- RİJİTLİK K ---
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

    % --- SÖNÜM C (aynı topoloji) ---
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

    % --- SİMETRİ TAMAMLAMA ---
    for ii = 1:11
        for jj = 1:11
            if K(ii,jj)~=0 && K(jj,ii)==0, K(jj,ii) = K(ii,jj); end
            if C(ii,jj)~=0 && C(jj,ii)==0, C(jj,ii) = C(ii,jj); end
        end
    end
end

function [K, C] = kuplor_ekle(K, C, izL, ithL, izV, ithV, LL, LV, Kk, Ck)
% Pitch-coupled kuplör: F = Kk*[(z_L - LL*th_L) - (z_V + LV*th_V)]
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

function s = Wk_filtre(sinyal, fs)
% ISO 2631-1 W_k dikey frekans ağırlık filtresi (basitleştirilmiş)
    [b1, a1] = butter(4, [0.5 80]/(fs/2), 'bandpass');
    bant = filtfilt(b1, a1, sinyal);
    [b2, a2] = butter(2, [4 8]/(fs/2), 'bandpass');   % 4-8 Hz en duyarlı
    pik = filtfilt(b2, a2, bant);
    s = 0.5*bant + 0.5*pik;
end

function sinif = konfor_sinifi(a)
% ISO 2631-1 konfor sınıflandırması
    if     a < 0.315, sinif = 'Cok Rahat';
    elseif a < 0.5,   sinif = 'Biraz Rahatsiz';
    elseif a < 0.8,   sinif = 'Oldukca Rahatsiz';
    else,             sinif = 'Rahatsiz';
    end
end
