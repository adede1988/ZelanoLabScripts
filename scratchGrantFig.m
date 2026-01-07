

curDat = load('R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\251006_OBE_NWU_RY_1\raw\raw_focusedBreathing_echem\raw_focusedBreathing_echem.mat');
curDat = curDat.curDat; 

data = curDat.rawData.trial{1};


outDat = struct; 
outDat.data = data;
outDat.labels = curDat.outLabs; 
outDat.fs = 2000; 
outDat.figs = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\251006_OBE_NWU_RY_1';
outDat.sessID =  '251006_OBE_NWU_RY_1';
P = struct; 
P.macroRemove = [];
P.spikeClean = false; 
P.spikeWin = 44; 
P.rspIDX = 1; 
P.rspFlip = 1; 
outDat.rspIDX = P.rspIDX;
outDat.rspFlip = P.rspFlip; 
outDat = preprocess_macros(outDat,P); 
outDat = process_respiration_breathing(outDat, P); 
bmObj = outDat.bmObj; 
 outDat.behDat = table();

    % Build time axis in seconds for index conversion
    tim = (1:size(outDat.data, 2)) / outDat.fs;

    % Column 2: onset time → sniffOnset and finalOnset (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 2));
    outDat.behDat.sniffOnset = idx;
    outDat.behDat.finalOnset = idx;

    % Column 12: condition
    outDat.behDat.condition = bmObj(:, 12);

    % Column 1: onset Y value
    outDat.behDat.Yonset = bmObj(:, 1);

    % Column 3: peak Y value
    outDat.behDat.inhaleMax = bmObj(:, 3);

    % Column 4: peak time → inMaxTim (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 4));
    outDat.behDat.inMaxTim = idx;

    % Column 5: end Y value
    outDat.behDat.Yend = bmObj(:, 5);

    % Column 6: end time → endTim (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 6));
    outDat.behDat.endTim = idx;

    % Column 7: length (end - onset)
    outDat.behDat.length = bmObj(:, 7);

    % Column 8: amp (peak Y - avg of two ends)
    outDat.behDat.amp = bmObj(:, 8);

    % Column 10: exhale peak Y value
    outDat.behDat.exhaleMin = bmObj(:, 10);

    % Column 11: exhale peak time → exMinTim (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 11));
    outDat.behDat.exMinTim = idx;

    % Column 14: index
    outDat.behDat.index = bmObj(:, 14);

     R = preprocess_respiration_wholetrace(outDat);
    plot_sniff_epochs(outDat, R);

  splitFreq = 2;          % Hz cutoff between "low" and "high" components
    hpOrder   = 4;           % 4th-order Butterworth for high-pass
    
    % Design high-pass for the spike-y part (> splitFreq)
    [b_hp, a_hp] = butter(hpOrder, splitFreq/(outDat.fs/2), 'high');
    
    % High-frequency component of the IC
    x_high = filtfilt(b_hp, a_hp, outDat.data(42:46,:).').';   % column

    for ii = 1:length(outDat.bmObj(:,1))
        % mask =  outDat.behDat.sniffOnset(ii)-500:outDat.behDat.sniffOnset(ii)+1000; 
        % x = outDat.data(46, mask);
        % fs = outDat.fs; 
        % % x: raw vector (Nx1 or 1xN)
        % % fs: sampling rate (Hz)
        % 
        % x = x(:);  % column
        % N = numel(x);
        % 
        % pad_ms = 20;
        % pad = round(pad_ms/1000 * fs);   % samples
        % 
        % thr = 200;  % <-- set this to something sensible for your units
        % above = x > thr;
        % 
        % % Find suprathreshold segments
        % d = diff([false; above; false]);
        % seg_starts = find(d == 1);
        % seg_ends   = find(d == -1) - 1;
        % 
        % % Peak index (max abs) within each segment
        % stim_idx = zeros(numel(seg_starts),1);
        % for k = 1:numel(seg_starts)
        %     s = seg_starts(k);
        %     e = seg_ends(k);
        %     [~, rel] = max(abs(x(s:e)));
        %     stim_idx(k) = s + rel - 1;
        % end
        % 
        % % Zero padded window around each peak
        % x_clean = x;
        % for k = 1:numel(stim_idx)
        %     c  = stim_idx(k);
        %     lo = max(1, c - pad);
        %     hi = min(N, c + pad);
        %     x_clean(lo:hi) = NaN;
        % end
        % 
        % figure
        % plot(x, '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
        % 
        % hold on
        %  % light grey dashed
        % plot(x_clean, 'k-', 'LineWidth', 1)
        % yyaxis right
        % plot(outDat.data(42, mask ))
        % xlim([0,1500])
        % legend({'x\_clean','x (raw)', 'respiration'}, 'Location', 'best')
        % title(ii)
        % 
        % xticks([0:200:1600])
        % xticklabels([-1000:400:2200])
        % 



% multi channel


mask = outDat.behDat.sniffOnset(ii)-500 : outDat.behDat.sniffOnset(ii)+1000;
fs   = outDat.fs;

chans = 46:48;
X = outDat.data(chans, mask);   % 3 x T
T = size(X,2);

% --- detect artifacts using channel 46 only ---
x_det = X(1,:).';               % channel 46 for detection (column)
N = numel(x_det);

pad_ms = 20;
pad = round(pad_ms/1000 * fs);

thr   = 200;                    % set for your units
above = x_det > thr;            % (use abs(x_det) > thr if polarity can flip)

d = diff([false; above; false]);
seg_starts = find(d == 1);
seg_ends   = find(d == -1) - 1;

stim_idx = zeros(numel(seg_starts),1);
for k = 1:numel(seg_starts)
    s = seg_starts(k);
    e = seg_ends(k);
    [~, rel] = max(abs(x_det(s:e)));
    stim_idx(k) = s + rel - 1;
end

% --- build a single NaN mask and apply to all 3 channels ---
nan_mask = false(N,1);
for k = 1:numel(stim_idx)
    c  = stim_idx(k);
    lo = max(1, c - pad);
    hi = min(N, c + pad);
    nan_mask(lo:hi) = true;
end

X_clean = X;                    % 3 x T
X_clean(:, nan_mask) = NaN;

% --- plot: all 3 channels with vertical offsets, plus grey originals ---
off = 20;
offsets = (0:numel(chans)-1) * off;

fig = figure('Color','w', 'visible', false); hold on

% --- left axis: ephys (offset) ---
yyaxis left
for iC = 1:numel(chans)
    y_raw   = X(iC,:)       + offsets(iC);
    y_clean = X_clean(iC,:) + offsets(iC);

    plot(y_raw,  '--', 'Color', [0.85 0.85 0.85], 'LineWidth', 1); % lighter grey original
    plot(y_clean, 'k-',  'LineWidth', 1);                          % black cleaned
end

xlim([0 1500])
ylim([-20 70])

% X axis labels in ms relative to inhale
xticks(0:200:1600)
xticklabels(-1000:400:2200)
xlabel('Time Relative to Inhale (ms)')

ylabel('Voltage (\muV)')

% Remove y tick labels (keep ticks off entirely)
yticks([])
% If you prefer to keep ticks but hide labels:
% yticklabels({})

% Make axes bolder + white panel
ax = gca;
ax.LineWidth = 1.5;
ax.Box = 'off';
ax.XColor = 'k';
ax.YColor = 'k';

% --- right axis: respiration (smoothed) ---
yyaxis right
resp = outDat.data(42, mask);
resp_s = smoothdata(resp, 'movmean', max(1, round(0.05*fs))); % 50 ms smoothing window
plot(resp_s, 'LineWidth', 3, 'Color', [0.13 0.55 0.13])       % forest green (RGB)
ax = gca;
ax.YColor = [0.13 0.55 0.13];
ax.LineWidth = 1.5;
yticks([])
ylabel('Respiration (au; inhale up)')

% Inhale time marker
% xline(500, '--', 'Color', 'r', 'LineWidth', 2);

% title(ii)

% --- scale bar: 20 uV in bottom-right quadrant (left axis units) ---
yyaxis left
bar_uV = 20;
xL = xlim; yL = ylim;

% place it inset from bottom-right
x_bar = xL(1) + 200;                 % 80 samples left from right edge
y_bar0 = yL(2) - 0.25*range(yL);     % a bit above bottom
line([x_bar x_bar], [y_bar0 y_bar0 + bar_uV], 'Color','k', 'LineWidth', 2);
text(x_bar + 10, y_bar0 + bar_uV/2, sprintf('%d \\muV', bar_uV), ...
    'Color','k', 'FontSize', 10, 'VerticalAlignment','middle');


% Set physical size: 3 in x 3 in
set(fig, 'Units','inches', 'Position',[1 1 3 3]);      % on-screen size (optional)
set(fig, 'PaperUnits','inches', 'PaperPosition',[0 0 3 3]);
set(fig, 'PaperSize',[3 3]);

% Save JPG at 300 dpi
tmp = ['sniffAligned_ephys_' num2str(ii) '.jpg'];
outFile = fullfile(outDat.figs, tmp);
print(fig, outFile, '-djpeg', '-r300');  % 300 DPI JPEG
% legend({'raw', 'clean', 'respiration'}, 'Location', 'northwest')



    end