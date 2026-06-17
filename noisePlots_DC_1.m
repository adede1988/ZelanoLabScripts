%% scratch DC_1 analysis: 

names = {'raw_task1', 'raw_task2', 'raw_task3', 'raw_task4'};
set(0, 'defaultfigurewindowstyle', 'docked')

allDat = cell(4,1); 
for ss = 1:4
    allDat{ss} = load(['R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\260122_OBE_NWU_DC_1\raw\'...
        names{ss} '\' names{ss} '.mat']); 
    allDat{ss} = allDat{ss}.curDat; 
end

for ss = 1:4
    
    fs = 2000;                          % Hz
    nChan = 38;
    nSamp = size(allDat{ss}.rawData.trial{1}, 2);
    
    % time in minutes
    tMin = (0:nSamp-1) / fs / 60;
    
    % channel labels
    chanLabels = { ...
        'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9','CP5','CP1', ...
        'Pz','P3','P7','O1','Oz','O2','P4','P8','TP10','CP6','CP2','Cz', ...
        'C4','T8','FT10','FC6','FC2','F4','F8','Fp2', ...
        'macro1','macro2','macro3','macro4','macro5','macro6'};
    
    % offsets for stacked plotting
    spacing = 200;                       % uV between traces
    offsets = (1:nChan) * spacing;
    
    % colors
    macroGreen = [0.20 0.60 0.20];     % leaf green-ish
    
    figure; hold on
    
    for ii = 1:nChan
        if ii <= 32
            thisColor = [0 0 0];       % scalp EEG = black
        else
            thisColor = macroGreen;    % macro channels = green
        end
        
        plot(tMin, allDat{ss}.rawData.trial{1}(ii,:) + offsets(ii), ...
            'Color', thisColor, 'LineWidth', 0.8);
    end
    
    % y-axis channel labels only
    set(gca, 'YTick', offsets, 'YTickLabel', chanLabels, ...
        'TickDir', 'out', 'Box', 'off')
    
    ylabel('Channel')
    xlabel('Time (min)')
    
    % leave blank area to the right for scale bar
    xPad = 0.08 * range(tMin);
    xlim([tMin(1), tMin(end) + xPad])
    ylim([-100 ii*spacing+100])
    % vertical 100 uV scale bar in blank area
    barX = tMin(end) + 0.04 * range(tMin);
    barY0 = offsets(3);                 % place near lower traces
    barH = 200;                         % 100 uV
    
    plot([barX barX], [barY0 barY0+barH], 'k', 'LineWidth', 2)
    plot([barX-0.003*range(tMin) barX+0.003*range(tMin)], [barY0 barY0], 'k', 'LineWidth', 2)
    plot([barX-0.003*range(tMin) barX+0.003*range(tMin)], [barY0+barH barY0+barH], 'k', 'LineWidth', 2)
    
    text(barX + 0.008*range(tMin), barY0 + barH/2, '200 \muV', ...
        'Rotation', 90, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left')
    
    % optional: put channel 1 at top
    set(gca, 'YDir', 'reverse')
    title(['Recording number: ' num2str(ss)])

end