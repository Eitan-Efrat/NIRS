function [dHbR, dHbO, fig] = CalcNIRS1(dataFile, SDS, tissueType, plotChannelIdx, ...
    extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
% CalcNIRS1 - Calculate and plot dHbR and dHbO changes from NIRS intensity data.
%
% [dHbR, dHbO, fig] = CalcNIRS1(dataFile, SDS, tissueType, plotChannelIdx, ...
%     extinctionCoefficientsFile, DPFperTissueFile, relDPFfile)
%
% Inputs:
%   dataFile: .mat file containing:
%              SD.Lambda   - a 2-element vector with wavelengths (nm)
%              t           - time vector (sec)
%              d           - intensity data; 40 columns (first 20 for wavelength1,
%                            last 20 for wavelength2)
%   SDS:        Source-detector separation (cm, positive scalar)
%   tissueType: Tissue type as a string (e.g., 'adult_forearm', 'adult_head', etc.)
%   plotChannelIdx: Vector (1 to 20) indicating channels to plot (if empty, no plot)
%   extinctionCoefficientsFile: CSV file with extinction coefficients (columns: wavelength, HbO2, HHb, etc.)
%   DPFperTissueFile: TXT file containing tissue and DPF values.
%   relDPFfile: CSV file with columns: wavelength and relDPFcoeff.
%
% Outputs:
%   dHbR: Time series of deoxyhemoglobin changes (n x 20 matrix)
%   dHbO: Time series of oxyhemoglobin changes (n x 20 matrix)
%   fig:  Handle to time-course plot figure (empty if no plotChannelIdx provided)
%
% This function calculates the optical density (OD) using:
%       OD = log10(I0./I(t))
% (where I0 is the baseline intensity at the first time point). It then uses a
% two-wavelength linear system based on extinction coefficients and effective path
% lengths to compute the changes in hemoglobin concentrations.
%
% Example:
%   [dHbR, dHbO, fig] = CalcNIRS1('FN_031_V2_Postdose2_Nback.mat', 3, 'adult_forearm', [1 2], ...
%                    'ExtinctionCoefficientsData.csv', 'DPFperTissue.txt', 'RelativeDPFCoefficients.csv');

%% Input Checks
% Check that file inputs are valid strings and exist.
if ~(ischar(dataFile) || isstring(dataFile))
    error('dataFile must be a string.');
end
if ~exist(dataFile, 'file')
    error('dataFile does not exist.');
end

if ~(isscalar(SDS) && isnumeric(SDS) && SDS > 0)
    error('SDS must be a positive numeric scalar.');
end

if ~(ischar(tissueType) || isstring(tissueType))
    error('tissueType must be a string.');
end

if nargin < 4 || isempty(plotChannelIdx)
    plotChannelIdx = [];
end

if ~(ischar(extinctionCoefficientsFile) || isstring(extinctionCoefficientsFile))
    error('extinctionCoefficientsFile must be a string.');
end
if ~exist(extinctionCoefficientsFile, 'file')
    error('extinctionCoefficientsFile does not exist.');
end

if ~(ischar(DPFperTissueFile) || isstring(DPFperTissueFile))
    error('DPFperTissueFile must be a string.');
end
if ~exist(DPFperTissueFile, 'file')
    error('DPFperTissueFile does not exist.');
end

if ~(ischar(relDPFfile) || isstring(relDPFfile))
    error('relDPFfile must be a string.');
end
if ~exist(relDPFfile, 'file')
    error('relDPFfile does not exist.');
end

%% Load Data File
% Load variables SD, t, and d from the MAT file.
data = load(dataFile);
if ~isfield(data, 'SD') || ~isfield(data, 't') || ~isfield(data, 'd')
    error('The MAT file must contain fields SD, t, and d.');
end
SD = data.SD;
t  = data.t;
d  = data.d;
if ~isfield(SD, 'Lambda') || numel(SD.Lambda) ~= 2
    error('SD must contain a 2-element vector Lambda.');
end
if size(d,2) ~= 40
    error('Intensity data d must have 40 columns (20 for each wavelength).');
end

% Split intensity data into two wavelengths:
d1 = d(:, 1:20);  % Intensity for wavelength 1
d2 = d(:, 21:40); % Intensity for wavelength 2

%% Load Supporting Tables
extTable = readtable(extinctionCoefficientsFile);
DPFTable = readtable(DPFperTissueFile, 'FileType', 'text');
relDPFTable = readtable(relDPFfile);

%% Compute Effective Path Lengths
% Extract measurement wavelengths.
lambda1 = SD.Lambda(1);
lambda2 = SD.Lambda(2);

% Get base DPF for the specified tissue type.
idx = strcmpi(DPFTable.Tissue, tissueType);
if ~any(idx)
    error('Specified tissue type not found in DPFperTissueFile.');
end
baseDPF = DPFTable.DPF(idx);

% Interpolate the relative DPF for each wavelength.
relDPF1 = interp1(relDPFTable.wavelength, relDPFTable.relDPFcoeff, lambda1, 'linear');
relDPF2 = interp1(relDPFTable.wavelength, relDPFTable.relDPFcoeff, lambda2, 'linear');

% Compute effective path lengths:
L_eff1 = SDS * baseDPF * relDPF1;
L_eff2 = SDS * baseDPF * relDPF2;

%% Calculate Optical Density (OD)
% Use the first time point as the baseline intensity.
I0_1 = d1(1, :);
I0_2 = d2(1, :);
OD1 = log10(repmat(I0_1, size(d1, 1), 1) ./ d1);
OD2 = log10(repmat(I0_2, size(d2, 1), 1) ./ d2);

%% Get Extinction Coefficients and Build Extinction Matrix
% Interpolate extinction coefficients for HbO2 and HHb at each wavelength.
extCoeff1_HbO = interp1(extTable.wavelength, extTable.HbO2, lambda1, 'linear');
extCoeff1_HbR = interp1(extTable.wavelength, extTable.HHb,  lambda1, 'linear');
extCoeff2_HbO = interp1(extTable.wavelength, extTable.HbO2, lambda2, 'linear');
extCoeff2_HbR = interp1(extTable.wavelength, extTable.HHb,  lambda2, 'linear');

% Build the extinction matrix (rows: wavelengths; columns: [HbR, HbO]).
A = [extCoeff1_HbR, extCoeff1_HbO; extCoeff2_HbR, extCoeff2_HbO];
A_inv = inv(A);

%% Solve for dHbR and dHbO
% Initialize output matrices.
nChannels = 20;             % There are 20 channels.
nTime = length(t);          % Number of time points.
dHbR = zeros(nTime, nChannels);
dHbO = zeros(nTime, nChannels);

% For each channel, solve the system:
%    [OD1/L_eff1; OD2/L_eff2] = A * [dHbR; dHbO]
for ch = 1:nChannels
    Y = [OD1(:, ch) / L_eff1, OD2(:, ch) / L_eff2]';  % 2 x nTime matrix
    X = A_inv * Y;
    dHbR(:, ch) = X(1, :)';
    dHbO(:, ch) = X(2, :)';
end

%% Plot Time Courses for Specified Channels (if requested)
if ~isempty(plotChannelIdx)
    fig = figure;
    for i = 1:length(plotChannelIdx)
        ch = plotChannelIdx(i);
        subplot(length(plotChannelIdx), 1, i);
        plot(t, dHbR(:, ch), 'b', t, dHbO(:, ch), 'r');
        xlabel('Time (s)');
        ylabel('d Concentration');
        title(sprintf('Channel %d', ch));
        legend('dHbR', 'dHbO');
    end
else
    fig = [];
end

%% FFT and SNR Calculation for Channel 1 dHbO
x = dHbO(:, 1);
N = length(x);
if N < 2
    error('Not enough data points for FFT.');
end
dt = t(2) - t(1);  % Assumes uniform sampling
Fs = 1 / dt;
X_fft = fft(x);
f = (0:N-1) * (Fs / N);
% Retain only positive frequencies
idxPositive = f <= Fs/2;
f = f(idxPositive);
X_mag = abs(X_fft(idxPositive));

% Define signal as the maximum magnitude below 2.5 Hz and noise as the mean above 2.5 Hz.
signal_level = max(X_mag(f < 2.5));
noise_level = mean(X_mag(f >= 2.5));
SNR = signal_level / noise_level;

figure;
plot(f, X_mag, 'k');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(sprintf('FFT of Ch1 dHbO (SNR = %.2f)', SNR));
end
