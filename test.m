%% MATLAB Script for NIRS Data Processing and Metabolic Index Calculation
% This script:
% - Reads two CSV files:
%       red_nirs_data.csv for the red wavelength (~630 nm)
%       nir_nirs_data.csv for the near-infrared (~950 nm)
%   Each file is assumed to have two columns: time (or sample index) and measured intensity.
%
% - Computes optical density differences using the modified Beer–Lambert law.
% - Solves for oxyhemoglobin (NHbO2) and deoxyhemoglobin (NHb) signals using the extinction
%   coefficient matrix inversion.
% - Applies a 2nd order Butterworth band-pass filter to extract the pulsatile (AC) components.
% - Uses the Hilbert transform to extract instantaneous amplitude and phase.
% - Computes arterial oxygen saturation (SaO2), phase difference (Δθ), and the metabolic index (MI).
% - Optionally applies an amplitude correction.

%% 1. Read Input CSV Files for Red and NIR Data
red_data = readmatrix('adc_red1_reading(in).csv');   % Red data (~630 nm)
nir_data = readmatrix('adc_ir1_reading(in).csv');     % NIR data (~950 nm)

t = red_data(:, 1);         % Time vector (assumes same time base for both files)
Iout1 = red_data(:, 2);     % Measured intensity for red wavelength (~630 nm)
Iout2 = nir_data(:, 2);     % Measured intensity for NIR (~950 nm)

I0_1 = Iout1(1);            % Baseline intensity for red (first sample)
I0_2 = Iout2(1);            % Baseline intensity for NIR (first sample)

%% 2. Compute Optical Density Differences (ΔOD)
deltaOD1 = log10(I0_1 ./ Iout1);  % For red (~630 nm)
deltaOD2 = log10(I0_2 ./ Iout2);  % For NIR (~950 nm)

%% 3. Solve for NIRS Signals (NHbO2 and NHb)
% Set extinction coefficients (example values; adjust as needed for your wavelengths)
epsilonHbO2 = [368, 1222];  % [ε_HbO2(red), ε_HbO2(NIR)]
epsilonHb   = [3750.12,  763.84];  % [ε_Hb(red),   ε_Hb(NIR)]

E = [epsilonHbO2(1), epsilonHb(1);
     epsilonHbO2(2), epsilonHb(2)];
E_inv = inv(E);

N_HbO2 = zeros(size(t));
N_Hb   = zeros(size(t));

for i = 1:length(t)
    deltaOD_i = [deltaOD1(i); deltaOD2(i)];
    N_temp = E_inv * deltaOD_i;
    N_HbO2(i) = N_temp(1);
    N_Hb(i)   = N_temp(2);
end

%% 4. Extract AC (Pulsatile) Components via Band-Pass Filtering
fs = 1 / (t(2) - t(1));  % Estimate the sampling frequency from time vector differences
Wn = [0.8, 10] / (fs/2); % Normalized passband frequencies (0.8 Hz to 10 Hz)
[b, a] = butter(2, Wn, 'bandpass');

N_HbO2_AC = filtfilt(b, a, N_HbO2);
N_Hb_AC   = filtfilt(b, a, N_Hb);

%% 5. Hilbert Transform for Amplitude and Phase Extraction
analytic_HbO2 = hilbert(N_HbO2_AC);
analytic_Hb   = hilbert(N_Hb_AC);

A_HbO2 = abs(analytic_HbO2);  % Instantaneous amplitude of oxyhemoglobin signal
A_Hb   = abs(analytic_Hb);    % Instantaneous amplitude of deoxyhemoglobin signal

phase_HbO2 = angle(analytic_HbO2);  % Instantaneous phase of oxyhemoglobin signal
phase_Hb   = angle(analytic_Hb);    % Instantaneous phase of deoxyhemoglobin signal

%% 6. Compute Arterial Oxygen Saturation and Phase Delay
SaO2 = A_HbO2 ./ (A_HbO2 + A_Hb + eps);  % SaO2 = A_HbO2/(A_HbO2+A_Hb)
deltaTheta = abs(phase_Hb - phase_HbO2);   % Phase difference Δθ
deltaTheta = mod(deltaTheta, pi);          % Wrap phase difference into [0, π]

%% 7. Compute the Metabolic Index (MI)
MI = SaO2 .* (1 - SaO2) .* deltaTheta;


%% 8. Optional: Amplitude Correction (set correction factor α = 1 if no correction)
alpha = ones(size(MI));  % Dummy correction factor (replace with your own computation if needed)
%alpha = 1;
MI_corrected = alpha .* MI;

%% 9. Plot Results
figure;
subplot(3,1,1);
plot(t, N_HbO2, 'b', t, N_Hb, 'r');
xlabel('Time (s)');
ylabel('NIRS Signal');
legend('NHbO2', 'NHb');
title('Computed NIRS Signals');

subplot(3,1,2);
plot(t, SaO2, 'k');
xlabel('Time (s)');
ylabel('SaO2');
title('Arterial Oxygen Saturation');

subplot(3,1,3);
plot(t, MI, 'b', t, MI_corrected, 'r--');
xlabel('Time (s)');
ylabel('Metabolic Index');
legend('MI', 'MI Corrected');
title('Metabolic Index Computation');
