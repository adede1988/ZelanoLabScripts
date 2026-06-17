load('R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\251008_EEG_NWU_JC\preProc\251008_EEG_NWU_JC_breathingPreProc.mat')

rspDat = chanDat.data(36,:);

rspDat = bandpass(rspDat, [.05, 1], 500); 
fs = 500; % sampling rate
fs_new = 25;
x_ds = resample(rspDat, fs_new, fs);

win_size = fs_new * 60;  % 60 sec
step = fs_new * 5;       % 5 sec



m = 2;
r = 0.2 * std(rspDat);     % standard choice

[sampen_vals, t_idx] = windowed_sampen(x_ds, win_size, step, m, r);

% plot
time = (1:length(x_ds)) / fs_new;

% --- colors / formatting ---
bgCol  = [26 24 56]/255;       % #1A1838
labCol = [255 234 177]/255;    % #FFEAB1
rspCol = [223 230 218]/255;    % light sage green

fs = 500;

plotSpecs = {
    4440,   135180,  'audio book',              'audio_book_respiration.jpg'
    142080, 302520,  'mindfulness meditation',  'mindfulness_meditation_respiration.jpg'
    917350, 1076540, 'fast meditation',         'fast_meditation_respiration.jpg'
    1078960,1208330, 'slow meditation',         'slow_meditation_respiration.jpg'
};
figSaveDir = 'G:\My Drive\cZelano\presentations\sharkTank2026';
for ii = 1:size(plotSpecs,1)

    % --- data ---
    seg = rspDat(plotSpecs{ii,1}:plotSpecs{ii,2});
    tim = (1:length(seg)) ./ fs;

    % --- figure / axes ---
    fig = figure('Color', bgCol, 'Position', [0, 0, 1000, 600]);
    ax = axes;
    hold(ax, 'on');

    ax.Color      = bgCol;
    ax.LineWidth  = 2.2;
    ax.FontSize   = 16;
    ax.FontWeight = 'bold';
    ax.FontName   = 'Dotum';
    ax.XColor     = labCol;
    ax.YColor     = labCol;
    ax.TickDir    = 'out';
    ax.TickLength = [0.018 0.018];
    box(ax, 'off');

    % --- plot respiration ---
    plot(tim, seg, '-', 'Color', rspCol, 'LineWidth', 3.2);

    xlim([0 max(tim)]);
    ylim([-400 400]);

    % --- labels ---
    xlabel('time (seconds)', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    ylabel('Respiration', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    title(plotSpecs{ii,3}, ...
        'FontSize', 22, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    set(gca, 'Layer', 'top');

    % --- save ---
    outFile = fullfile(figSaveDir, plotSpecs{ii,4});
    exportgraphics(fig, outFile, ...
        'BackgroundColor', bgCol, ...
        'Resolution', 300);
end

figure;
plot(rspDat);
xlim([1078960 1208330])


figure
plot(time, x_ds)
xlim([])

xline(time(ceil(chanDat.TTL ./ (500/25))+1 ) )
yyaxis right
plot(t_idx/fs_new, sampen_vals, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Signal / SampEn');
legend('Respiration', 'Sample Entropy');