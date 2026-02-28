function chanDat = addBreathHRV_fromRRint(chanDat, rawDat, fs, varargin)
% addBreathHRV_fromRRint
%   Computes continuous HRV measures from rawDat.RRint (seconds, continuous, cleaned)
%   and stores per-breath means into chanDat.behDat.
%
% Inputs required in workspace/structs:
%   rawDat.RRint                1 x time   (seconds)
%   rawDat.behDat.finalOnset    breaths x 1 sample indices into RRint
%   rawDat.behDat.length        breaths x 1 (seconds)
%   fs                          scalar Hz (e.g., 500)
%
% Outputs added to chanDat.behDat (per breath):
%   HRV_RMSSD30   (seconds)  -- 30s windowed RMSSD (from downsampled RR)
%   HRV_SDNN30    (seconds)  -- 30s windowed SDNN  (from downsampled RR)
%   HRV_RSAamp    (bpm)      -- RSA amplitude (envelope of resp-band HR)
%   HR_mean       (bpm)      -- mean HR in breath
%   RR_mean       (seconds)  -- mean RR in breath
%
% Options (name-value):
%   'WinSec'      (default 30)
%   'HRVfs'       (default 10)         % internal rate for HRV traces (Hz); must divide fs
%   'RSAbandHz'   (default [0.04 0.5]) % respiratory band for RSA on HR
%   'RSAon'       (default 'HR')       % 'HR' or 'RR'
%   'EnvLowpassHz'(default 2)          % anti-alias LP before downsample

% ----------------- options -----------------
p = inputParser;
p.addParameter('WinSec', 30, @(x)isscalar(x)&&x>0);
p.addParameter('HRVfs',  10, @(x)isscalar(x)&&x>0);
p.addParameter('RSAbandHz', [0.04 0.5], @(v)isnumeric(v)&&numel(v)==2&&v(1)>0&&v(2)>v(1));
p.addParameter('RSAon', 'HR', @(s)ischar(s)||isstring(s));
p.addParameter('EnvLowpassHz', 2, @(x)isscalar(x)&&x>0);
p.parse(varargin{:});
opt = p.Results;

WinSec    = opt.WinSec;
hrvFs     = opt.HRVfs;
RSAbandHz = opt.RSAbandHz(:).';
RSAon     = upper(string(opt.RSAon));
lpHz      = opt.EnvLowpassHz;

% ----------------- pull inputs -----------------
RR = double(rawDat.RRint(:).');  % 1 x T (seconds)
T  = numel(RR);

onset = rawDat.behDat.finalOnset(:);
blen  = rawDat.behDat.length(:);
nBreaths = numel(onset);

if height(chanDat.behDat) ~= nBreaths
    error('chanDat.behDat height (%d) must match number of breaths in rawDat.behDat (%d).', ...
        height(chanDat.behDat), nBreaths);
end

% ----------------- downsample to HRVfs -----------------
ds = fs / hrvFs;
if abs(ds - round(ds)) > 1e-9
    error('HRVfs must evenly divide fs. fs=%g, HRVfs=%g', fs, hrvFs);
end
ds = round(ds);
fs_ds = fs / ds;

% anti-alias LP before decimation (RR variations are slow)
lpHz_use = min(lpHz, fs/2 - 1);
lpFilt = designfilt('lowpassiir','FilterOrder',4, ...
    'HalfPowerFrequency', lpHz_use, 'SampleRate', fs, 'DesignMethod','butter');
RR(isnan(RR)) = median(RR, 'omitnan');
RR_lp = filtfilt(lpFilt, RR);
RR_ds = RR_lp(1:ds:end);              % 1 x N
N     = numel(RR_ds);

% ----------------- continuous HRV traces (30s window) -----------------
wlen = max(3, round(WinSec * fs_ds));
if mod(wlen,2)==0, wlen = wlen + 1; end

% SDNN over RR (seconds)
sdnn_t = movstd(RR_ds, wlen, 0, 'omitnan', 'Endpoints','shrink');

% RMSSD over successive RR differences (seconds)
dRR = diff(RR_ds);
rmssd_core = sqrt(movmean(dRR.^2, max(2,wlen-1), 'omitnan', 'Endpoints','shrink'));
rmssd_t = [NaN, rmssd_core];          % align roughly to RR_ds time base (1 x N)

% ----------------- RSA (amplitude of resp-band oscillation) -----------------
switch RSAon
    case "HR"
        sig = 60 ./ RR_ds;            % bpm
    case "RR"
        sig = RR_ds;                  % seconds
    otherwise
        error('RSAon must be ''HR'' or ''RR''.');
end

nyq = fs_ds/2;
b1 = max(RSAbandHz(1), 0.001);
b2 = min(RSAbandHz(2), nyq-0.001);
if b1 >= b2
    error('Invalid RSAbandHz after clamping to Nyquist. fs_ds=%g, band=[%g %g].', fs_ds, b1, b2);
end

bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1', b1, 'HalfPowerFrequency2', b2, ...
    'SampleRate', fs_ds, 'DesignMethod','butter');

sig_bp = filtfilt(bpFilt, sig);
rsa_t  = abs(hilbert(sig_bp));        % amplitude envelope (bpm if RSAon='HR')

% ----------------- per-breath means (map indices to downsampled grid) -----------------
rmssd_b = nan(nBreaths,1);
sdnn_b  = nan(nBreaths,1);
rsa_b   = nan(nBreaths,1);
hr_b    = nan(nBreaths,1);
rr_b    = nan(nBreaths,1);

for b = 1:nBreaths
    s0 = onset(b);
    if ~isfinite(s0) || s0 < 1 || s0 > T, continue; end

    e0 = s0 + round(blen(b)*fs) - 1;
    e0 = min(max(e0,1), T);

    % map to downsampled indices
    sD = floor((s0-1)/ds) + 1;
    eD = floor((e0-1)/ds) + 1;
    sD = max(1, min(N, sD));
    eD = max(1, min(N, eD));
    if eD < sD, continue; end

    rmssd_b(b) = mean(rmssd_t(sD:eD), 'omitnan');
    sdnn_b(b)  = mean(sdnn_t(sD:eD),  'omitnan');
    rsa_b(b)   = mean(rsa_t(sD:eD),   'omitnan');

    % optional breath means from RR/HR
    rr_b(b)    = mean(RR_ds(sD:eD),   'omitnan');
    hr_b(b)    = mean(60./RR_ds(sD:eD), 'omitnan');
end

% ----------------- store into chanDat.behDat -----------------
chanDat.behDat.HRV_RMSSD30 = rmssd_b;     % seconds
chanDat.behDat.HRV_SDNN30  = sdnn_b;      % seconds
chanDat.behDat.HRV_RSAamp  = rsa_b;       % bpm if RSAon='HR', else seconds
chanDat.behDat.HR_mean     = hr_b;        % bpm
chanDat.behDat.RR_mean     = rr_b;        % seconds

% stash some meta (optional)
chanDat.hrv.meta.WinSec    = WinSec;
chanDat.hrv.meta.HRVfs     = fs_ds;
chanDat.hrv.meta.RSAbandHz = [b1 b2];
chanDat.hrv.meta.RSAon     = char(RSAon);

end