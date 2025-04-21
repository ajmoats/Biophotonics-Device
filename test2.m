%% MATLAB Script for NIRS Data Processing and Metabolic Index Calculation
% This script:
%   1. Reads two CSV files:
%         red_nirs_data.csv (Red wavelength ~630 nm)
%         nir_nirs_data.csv (Near infrared ~950 nm)
%      Each file is assumed to have two columns: time (or sample index) and measured intensity.
%   2. Computes optical density (OD) differences using the modified Beer–Lambert law.
%   3. Solves for NIRS signals (NHbO2 and NHb) using an extinction coefficient matrix inversion.
%   4. Removes body tremor noise using a low-pass filter (cutoff = 2 Hz) to eliminate components
%      at frequencies above ~2 Hz (including tremors at 3, 5, 9, and 13 Hz).
%   5. Applies the Hilbert transform to obtain amplitude and instantaneous phase.
%   6. Computes arterial oxygen saturation (SaO2) and the phase difference (Δθ).
%   7. Calculates the metabolic index (MI).
%   8. Optionally applies an amplitude correction.
%
% Adjust file names, filter cutoff frequencies, extinction coefficients, and other parameters as needed.

%% 1. Read Input CSV Files for Red and NIR Data
red_data = readmatrix('adc_red1_reading(in).csv');   % Red data (~630 nm)
nir_data = readmatrix('adc_ir1_reading(in).csv');     % NIR data (~950 nm)

t = red_data(:, 1);         % Time vector (assumes both files share the same time base)
Iout1 = red_data(:, 2);     % Measured intensity for red wavelength (~630 nm)
Iout2 = nir_data(:, 2);     % Measured intensity for NIR wavelength (~950 nm)

I0_1 = Iout1(1);            % Baseline intensity for red (first sample)
I0_2 = Iout2(1);            % Baseline intensity for NIR (first sample)

%% 2. Compute Optical Density Differences (∆OD)
% Using the modified Beer–Lambert law:
%   ∆OD(λ,t) = log10( Iout(λ,t0) / Iout(λ,t) )
deltaOD1 = log10(I0_1 ./ Iout1);  % For red
deltaOD2 = log10(I0_2 ./ Iout2);  % For NIR

%% 3. Solve for NIRS Signals (NHbO2 and NHb)
% Set extinction coefficients (example values; modify these as appropriate for your wavelengths)
epsilonHbO2 = [942, 1214];  % [ε_HbO2(red), ε_HbO2(NIR)]
epsilonHb   = [6509.6,  693.44];  % [ε_Hb(red),   ε_Hb(NIR)]

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

%% 4. Remove Body Tremor Noise using a Low-Pass Filter
% Here we remove frequencies above ~2 Hz (which include tremor components at 3, 5, 9, and 13 Hz)
% because the desired pulsatile (e.g., heartbeat) signal is at a lower frequency.
fs = 1 / (t(2) - t(1));  % Estimate the sampling frequency from t
cutoff = 2;              % Low-pass cutoff frequency in Hz
Wn = cutoff / (fs/2);    % Normalize cutoff frequency to the Nyquist frequency
[b, a] = butter(2, Wn, 'low');  % 2nd order Butterworth low-pass filter

N_HbO2_filtered = filtfilt(b, a, N_HbO2);
N_Hb_filtered   = filtfilt(b, a, N_Hb);

%% 5. Hilbert Transform for Amplitude and Phase Extraction
analytic_HbO2 = hilbert(N_HbO2_filtered);
analytic_Hb   = hilbert(N_Hb_filtered);

A_HbO2 = abs(analytic_HbO2);     % Instantaneous amplitude of oxyhemoglobin signal
A_Hb   = abs(analytic_Hb);       % Instantaneous amplitude of deoxyhemoglobin signal

phase_HbO2 = angle(analytic_HbO2);  % Instantaneous phase of oxyhemoglobin signal
phase_Hb   = angle(analytic_Hb);    % Instantaneous phase of deoxyhemoglobin signal

%% 6. Compute Arterial Oxygen Saturation (SaO2) and Phase Difference (Δθ)
SaO2 = A_HbO2 ./ (A_HbO2 + A_Hb + eps);  % Compute SaO2
deltaTheta = abs(phase_Hb - phase_HbO2);   % Compute phase difference Δθ
deltaTheta = mod(deltaTheta, pi);          % Wrap phase differences into [0, π]

%% 7. Compute the Metabolic Index (MI)
% According to Eq. (32) from the paper:
%   MI(t) = SaO2(t) · (1 - SaO2(t)) · |Δθ(t)|
MI = SaO2 .* (1 - SaO2) .* deltaTheta;

%% 8. Optional: Amplitude Correction (if applicable)
alpha = ones(size(MI));      % Dummy correction factor (set to 1 if no correction)
MI_corrected = alpha .* MI;

%% 9. Plot the Results
figure;
subplot(4,1,1);
plot(t, N_HbO2, 'b', t, N_Hb, 'r');
xlabel('Time (s)'); ylabel('NIRS Signal');
legend('NHbO2', 'NHb');
title('Computed NIRS Signals (Before Tremor Filtering)');

subplot(4,1,2);
plot(t, N_HbO2_filtered, 'b', t, N_Hb_filtered, 'r');
xlabel('Time (s)'); ylabel('Filtered NIRS Signal');
legend('NHbO2 filtered', 'NHb filtered');
title('NIRS Signals after Low-Pass Filtering (Tremor Removal)');

subplot(4,1,3);
plot(t, MI, 'b', t, MI_corrected, 'r--');
xlabel('Time (s)'); ylabel('Metabolic Index');
legend('MI', 'MI Corrected');
title('Metabolic Index Computation');

subplot(4,1,4);
plot(t, SaO2, 'k');
xlabel('Time (s)');
ylabel('SaO2');
title('Arterial Oxygen Saturation');
