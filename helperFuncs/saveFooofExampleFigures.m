function hFigs = saveFooofExampleFigures(chanDat, saveDir, bgCol, labCol)
% saveFooofExampleFigures
%
% Creates and saves three example figures:
%   1) log-scaled power spectrum (flat + aperiodic)
%   2) log-scaled aperiodic
%   3) flat alone
%
% INPUTS
%   chanDat : struct containing
%       chanDat.fooof.spectra_flat_log10
%       chanDat.fooof.aperiodic_log10
%       chanDat.tf.frex
%
%   saveDir : folder to save figures into
%   bgCol   : background color, e.g. [1 1 1]
%   labCol  : axis/label color, e.g. [0 0 0]
%
% OUTPUT
%   hFigs   : struct with handles to the three figures
%
% Notes:
% - Uses channel index 2, matching your scratch code.
% - Uses only frequency indices 120:end.
% - Assumes the spectra are already in log10 units.
% - Frequency axis is displayed on a log scale for cleaner spectrum plotting.
%   If you want linear x-axis, remove: set(ax,'XScale','log')

if nargin < 2 || isempty(saveDir)
    saveDir = pwd;
end
if nargin < 3 || isempty(bgCol)
    bgCol = [1 1 1];
end
if nargin < 4 || isempty(labCol)
    labCol = [0 0 0];
end

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% -------------------------
% Pull relevant data
% -------------------------
freqIdx  = 120:numel(chanDat.tf.frex);
plotFrex = chanDat.tf.frex(freqIdx);

flatRaw = squeeze(chanDat.fooof.spectra_flat_log10(:,2,freqIdx));
aperRaw = squeeze(chanDat.fooof.aperiodic_log10(:,2,freqIdx));

% Make sure spectra are [nFreq x nSpectra] for plotting
if isvector(flatRaw)
    flatRaw = reshape(flatRaw, 1, []);
end
if isvector(aperRaw)
    aperRaw = reshape(aperRaw, 1, []);
end

if size(flatRaw,1) == numel(plotFrex)
    flatSpec = flatRaw;
elseif size(flatRaw,2) == numel(plotFrex)
    flatSpec = flatRaw';
else
    error('spectra_flat_log10 dimensions do not match chanDat.tf.frex(120:end).');
end

if size(aperRaw,1) == numel(plotFrex)
    aperSpec = aperRaw;
elseif size(aperRaw,2) == numel(plotFrex)
    aperSpec = aperRaw';
else
    error('aperiodic_log10 dimensions do not match chanDat.tf.frex(120:end).');
end

fullSpec = flatSpec + aperSpec;

% -------------------------
% Colors for traces / means
% -------------------------
allTraceCol = 0.75 .* [1 1 1];  % light gray for individual spectra
fullMeanCol = [0 0 0];          % black
aperMeanCol = [0.85 0.33 0.10]; % orange-ish
flatMeanCol = [0 0.45 0.74];    % blue-ish

% -------------------------
% Make + save figures
% -------------------------
hFigs = struct();

hFigs.full = local_makeAndSaveFig( ...
    plotFrex, fullSpec, ...
    'Log-scaled power spectrum (flat + aperiodic)', ...
    fullMeanCol, allTraceCol, ...
    fullfile(saveDir, 'example_logScaledPower_flatPlusAperiodic'), ...
    bgCol, labCol);

hFigs.aperiodic = local_makeAndSaveFig( ...
    plotFrex, aperSpec, ...
    'Log-scaled aperiodic', ...
    aperMeanCol, allTraceCol, ...
    fullfile(saveDir, 'example_logScaledAperiodic'), ...
    bgCol, labCol);

hFigs.flat = local_makeAndSaveFig( ...
    plotFrex, flatSpec, ...
    'Flat alone', ...
    flatMeanCol, allTraceCol, ...
    fullfile(saveDir, 'example_flatOnly'), ...
    bgCol, labCol);

end


function hFig = local_makeAndSaveFig(plotFrex, specMat, figTitle, meanCol, traceCol, fileStem, bgCol, labCol)

hFig = figure('Position', [0, 0, 600, 370], 'Color', bgCol);
ax = axes('Parent', hFig);
hold(ax, 'on');

% all spectra
plot(ax, plotFrex, specMat, 'Color', traceCol, 'LineWidth', 0.75);

% mean spectrum
plot(ax, plotFrex, mean(specMat, 2, 'omitnan'), 'Color', meanCol, 'LineWidth', 2.8);

% aesthetics
ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
box(ax, 'off');

set(ax, 'XScale', 'log');
xlim(ax, [plotFrex(1) plotFrex(end)]);
xline(37.5, 'color', [203, 157, 6]./255, 'linestyle', '--', 'linewidth', 6)
xlabel(ax, 'Frequency (Hz)', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);
ylabel(ax, 'log_{10} power', ...
    'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);
% title(ax, figTitle, ...
%     'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Dotum', 'Color', labCol);

set(ax, 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Dotum');

% save
exportgraphics(hFig, [fileStem '.png'], 'Resolution', 300);
% savefig(hFig, [fileStem '.fig']);

end