%% MATLAB Script for NIRS Data Processing, MI Computation and Hard‑Coded BGL Estimation
% 1. Read two CSVs (red ~630 nm, NIR ~950 nm)
% 2. Compute ∆OD via Beer–Lambert
% 3. Invert extinction matrix → NHbO2, NHb
% 4. Low‑pass filter (2 Hz) to remove tremors
% 5. Hilbert transform → amplitudes & phases
% 6. Compute SaO2 & phase delay Δθ
% 7. Compute MI
% 8. Apply dummy α=1 correction
% 9. Apply hard‑coded calibration: BGL_est = a*MI_corr + b
% 10. Plot signals, MI, SaO2, and estimated BGL

%% 1. Read Input CSVs
red_data = readmatrix('adc_red1_reading(in).csv');
nir_data = readmatrix('adc_ir1_reading(in).csv');
t     = red_data(:,1);
Iout1 = red_data(:,2);
Iout2 = nir_data(:,2);
I0_1 = Iout1(1);
I0_2 = Iout2(1);

%% 2. Optical Density
deltaOD1 = log10(I0_1 ./ Iout1);
deltaOD2 = log10(I0_2 ./ Iout2);

%% 3. Invert Extinction Matrix
epsilonHbO2 = [942, 1214];     % [ε_HbO2(red), ε_HbO2(NIR)]
epsilonHb   = [6509.6, 693.44]; % [ε_Hb(red), ε_Hb(NIR)]
E     = [epsilonHbO2(1), epsilonHb(1);
         epsilonHbO2(2), epsilonHb(2)];
E_inv = inv(E);
N_HbO2 = zeros(size(t));
N_Hb   = zeros(size(t));
for i = 1:length(t)
    v = E_inv * [deltaOD1(i); deltaOD2(i)];
    N_HbO2(i) = v(1);
    N_Hb(i)   = v(2);
end

%% 4. Low‑Pass Filter (2 Hz)
fs = 1/(t(2)-t(1));
[b,a] = butter(2, 2/(fs/2), 'low');
N_HbO2_f = filtfilt(b,a, N_HbO2);
N_Hb_f   = filtfilt(b,a, N_Hb);

%% 5. Hilbert → Amplitude & Phase
h_HbO2 = hilbert(N_HbO2_f);
h_Hb   = hilbert(N_Hb_f);
A_HbO2 = abs(h_HbO2);
A_Hb   = abs(h_Hb);
p_HbO2 = angle(h_HbO2);
p_Hb   = angle(h_Hb);

%% 6. SaO2 & Phase Delay
SaO2      = A_HbO2 ./ (A_HbO2 + A_Hb + eps);
deltaTheta = abs(p_Hb - p_HbO2);
deltaTheta = mod(deltaTheta, pi);

%% 7. Metabolic Index
MI = SaO2 .* (1 - SaO2) .* deltaTheta;

%% 8. α‑Correction (dummy = 1)
alpha = ones(size(MI));
MI_corr = alpha .* MI;

%% 9. Hard‑Coded Calibration → BGL_est
a = 100;   % slope (mg/dL per MI‑unit)
b =  40;   % intercept (mg/dL)
BGL_est = a * MI_corr + b;

%% 10. Plot Everything
figure;
subplot(5,1,1);
plot(t, N_HbO2, 'b', t, N_Hb, 'r');
legend('NHbO2','NHb'); ylabel('NIRS');
title('Raw NIRS Signals');

subplot(5,1,2);
plot(t, N_HbO2_f, 'b', t, N_Hb_f, 'r');
legend('NHbO2\_f','NHb\_f'); ylabel('Filtered');
title('After Low‑Pass (Tremor Removal)');

subplot(5,1,3);
plot(t, MI_corr);
ylabel('MI'); title('Metabolic Index (corrected)');

subplot(5,1,4);
plot(t, SaO2, 'k');
ylabel('SaO2'); title('Arterial Oxygen Saturation');

subplot(5,1,5);
plot(t, BGL_est, 'm');
xlabel('Time (s)');
ylabel('BGL (mg/dL)');
title('Estimated Blood Glucose (hard‑coded calib)');
