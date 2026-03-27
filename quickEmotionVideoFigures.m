%% scratch analysis of emotional Video task
datPre = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\CHANDAT_processed\';
CP = load(fullfile(datPre, '251009_OBE_NWU_CP_1_macro_EmotionalMovieTask_49.mat'));
TI = load(fullfile(datPre, '250904_OBE_NWU_TI_macro_EmotionalMovieTask_48.mat'));
AS = load(fullfile(datPre, '250225_OBE_NWU_AS_4_macro_EmotionalMovieTask_13.mat'));

CP = CP.chanDat; 
TI = TI.chanDat; 
AS = AS.chanDat; 

CP.use(69) = 0; 
CP.use(117) = 0; 
TI.use(506) = 0; 
TI.use(439) = 0; 

allDat   = {CP, TI, AS};
subjName = {'CP', 'TI', 'AS'};
condName = {'neutral', 'happy', 'sad'};

%% TF power over respiration 

for s = 1:3
    
    chanDat = allDat{s};
    figure('Name', subjName{s}, 'Color', 'w');
    
    for c = 1:3
        
        subplot(3,1,c)
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        
        imagesc( squeeze(mean(chanDat.tf.powZ(useIdx,:,:), 1, 'omitnan'))' )
        set(gca, 'YDir', 'normal')
        clim([-10 15])
        colorbar
        
        % frequency labels
        yticks(50:50:300)
        yticklabels(round(chanDat.tf.frex(50:50:300),1))
        
        % overlay mean breath trace on right axis
        yyaxis right
        plot(squeeze(mean(chanDat.tf.breathSeg(useIdx,:), 1, 'omitnan')), ...
            'LineWidth', 2)
        
        title([subjName{s} ' ' condName{c} ' video TF power across respiration'])
        
        if c == 3
            xlabel('Time')
        end
        
    end
end



%% TF ITPC over breath: 
figure('Name', 'itpc at resp phase 18', 'Color', 'w');
for c = 1:3
    
    subplot(3,1,c)
    hold on 
    for s = 1:3
        chanDat = allDat{s};
        
        
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        itpc = squeeze(abs(mean(exp(1i*chanDat.tf.phase(useIdx,:,:)), 1, 'omitnan')))' ;
        fx = chanDat.tf.frex; 

        plot(fx(fx>1 & fx<20), itpc(fx>1 & fx<20,18))
       
        
   
        
    end
    legend(subjName)
end



%% TF power locked to gamma peak
% --- theme colors ---
bgCol  = [26 24 56]/255;      % #1A1838
labCol = [255 234 177]/255;   % #FFEAB1

for s = 1:3
    
    chanDat = allDat{s};
    fig = figure('Name', subjName{s}, ...
        'Color', bgCol, ...
        'Position', [80 80 950 900], ...
        'InvertHardcopy', 'off');
    
    for c = 1:3
        
        ax = subplot(3,1,c);
        hold(ax, 'on')
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        
        trials = chanDat.trial.data(useIdx,:); 
        gamMed = median(chanDat.gammaBurst.freqHz(useIdx)); 
        
        trials = bandpass(trials.', [gamMed-5 gamMed+5], 500);   % [time x trials]
        phase  = angle(hilbert(trials)); 
        idxVals = chanDat.gammaBurst.t0_idx(useIdx) + 1000; 
        
        nTrials = numel(idxVals);
        win = -250:250;
        
        for ii = 1:nTrials
            searchIdx = idxVals(ii)-10 : idxVals(ii)+10;
            [~, tmp] = min(abs(phase(searchIdx, ii)));
            idxVals(ii) = searchIdx(tmp);
        end
        
        erpi = idxVals + win;   % [nTrials x nWin]
        trialNum = repmat((1:nTrials)', 1, numel(win));
        linIdx   = sub2ind(size(trials), erpi, trialNum);
        erp      = trials(linIdx);
        
        % ---------------- LEFT AXIS ----------------
        yyaxis left
        
        imagesc(ax, -500:50:500, [], ...
            squeeze(mean(chanDat.gammaLockTF.powZ_primary(useIdx,:,:), 1, 'omitnan'))')
        set(ax, 'YDir', 'normal')
        clim(ax, [-10 15])
        colormap(ax, turbo)
        
        yticks(ax, 10:10:100)
        yticklabels(ax, string(round(chanDat.gammaLockTF.frexSel(10:10:100),1)))
        
        ylabel(ax, 'frequency (Hz)', ...
            'FontSize', 18, ...
            'FontWeight', 'bold', ...
            'FontName', 'Dotum', ...
            'Color', labCol)
        ylim([0 100])
        % ---------------- RIGHT AXIS ----------------
        yyaxis right
        plot(ax, -500:2:500, mean(erp,1), 'k', 'LineWidth', 1.5)
        
        % optional right-axis label
        ylabel(ax, 'bandpassed signal', ...
            'FontSize', 18, ...
            'FontWeight', 'bold', ...
            'FontName', 'Dotum', ...
            'Color', labCol)
        
        % ---------------- COMMON STYLING ----------------
        ax.Color      = bgCol;
        ax.LineWidth  = 2.2;
        ax.FontSize   = 16;
        ax.FontWeight = 'bold';
        ax.FontName   = 'Dotum';
        ax.XColor     = labCol;
        ax.TickDir    = 'out';
        ax.TickLength = [0.018 0.018];
        box(ax, 'off');
        
        % explicitly style left and right y-axes
        ax.YAxis(1).Color = labCol;     % left axis color
        ax.YAxis(2).Color = labCol;    % right axis color
        xlim([-500 500])
        % colorbar styling
        cb = colorbar;
        cb.Color      = labCol;
        cb.LineWidth  = 2;
        cb.FontSize   = 14;
        cb.FontWeight = 'bold';
        cb.FontName   = 'Dotum';
        cb.Label.String = 'power (z)';
        cb.Label.Color = labCol;
        cb.Label.FontSize = 16;
        cb.Label.FontWeight = 'bold';
        cb.Label.FontName = 'Dotum';
        
        title([subjName{s} ' ' condName{c} ' video TF power at \gamma peak'], ...
            'FontSize', 18, ...
            'FontWeight', 'bold', ...
            'FontName', 'Dotum', ...
            'Color', labCol)

        if c == 3
            xlabel(ax, 'Time relative to \gamma peak (ms)', ...
                'FontSize', 18, ...
                'FontWeight', 'bold', ...
                'FontName', 'Dotum', ...
                'Color', labCol)
        end
        
        set(ax, 'Layer', 'top');
    end
end

%% TF ITPC locked to gamma peaks

for s = 1:3
    
    chanDat = allDat{s};
    figure('Name', subjName{s}, 'Color', 'w');
    
    for c = 1:3
        
        subplot(3,1,c)
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        
        imagesc( squeeze(abs(mean(exp(1i*chanDat.gammaLockTF.phase_primary(useIdx,:,:)), 1, 'omitnan')))' )
        set(gca, 'YDir', 'normal')
        clim([.05 .3])
        colorbar
        
        % frequency labels
        yticks(10:10:100)
        yticklabels(round(chanDat.gammaLockTF.frexSel(10:10:100),1))
        xticks(1:2:21)
        xticklabels(-500:100:500)
       
        
        title([subjName{s} ' ' condName{c} ' video ITPC locked to \gamma peak'])
        
        if c == 3
            xlabel('Time relative to \gamma peak (ms)')
        end
        
    end
end


%% gamma to breath phase consistency

for s = 1:3
    
    chanDat = allDat{s};
    figure('Name', subjName{s}, 'Color', 'w');
    
    for c = 1:3
        
        subplot(3,1,c)
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        
        polarhistogram(chanDat.pac.diag.fb_phase(useIdx,:,:), 12 )
      
       
        
        title([subjName{s} ' ' condName{c} ' video \gamma peak respiratory phase'])
        
       
        
    end
end

%% gamma peak ERP
figure('Name', 'ERP to gammaPeak', 'Color', 'w');

win = -500:500;
condName = {'neutral','happy','sad'};

for c = 1:3
    
    subplot(3,1,c)
    hold on
    
    for s = 1:3
        chanDat = allDat{s};
        
        useIdx   = chanDat.use == 1 & chanDat.behDat.condition == c;
        trialIdx = find(useIdx);
        peakIdx  = chanDat.pac_peaks.peakIDX_local(useIdx);
        
        erpMat = nan(numel(trialIdx), numel(win));
        
        for t = 1:numel(trialIdx)
            thisTrial = trialIdx(t);
            thisPeak  = peakIdx(t);
            sampIdx   = thisPeak + win;
            
            if all(sampIdx >= 1) && all(sampIdx <= size(chanDat.trial.data,2))
                erpMat(t,:) = chanDat.trial.data(thisTrial, sampIdx);
            end
        end
        
        plot(win, mean(erpMat,1,'omitnan'), 'LineWidth', 2)
    end
    
    xline(0,'k--')
    xlabel('Samples from gamma peak')
    ylabel('Amplitude')
    title(['ERP to gamma peak: ' condName{c}])
    legend(subjName, 'Location', 'best')
end


%% Theta gamma coupling strength

for s = 1:3
    
    chanDat = allDat{s};
    figure('Name', subjName{s}, 'Color', 'w');
    
    for c = 1:3
        
        subplot(3,1,c)
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        
        imagesc( squeeze(mean(chanDat.pacTheta.pac(useIdx,1:50,:, 10), 1, 'omitnan'))' )
        set(gca, 'YDir', 'normal')
        clim([-.5 .5])
        colorbar
        
        % frequency labels
        yticks(5:5:50)
        yticklabels(round(chanDat.pacTheta.PACfrex(5:5:50),1))
       
        ylabel('phase freq (Hz)')
       
        
        title([subjName{s} ' ' condName{c} ' video theta gamma coupling'])
        
        if c == 3
            xlabel('breath template sampled time')
        end
        
    end
end


%% gamma envelope power

for s = 1:3
    
    chanDat = allDat{s};
    figure('Name', subjName{s}, 'Color', 'w');
    
    for c = 1:3
        
        subplot(3,1,c)
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        % [~, order] = sort(chanDat.)
        imagesc( squeeze(abs(mean(exp(1i*chanDat.gammaEnv.phase(useIdx,:,:)), 1, 'omitnan')))' )
        set(gca, 'YDir', 'normal')
        clim([.05 .3])
        colorbar
        
        % frequency labels
        yticks(2:2:20)
        yticklabels(round(chanDat.gammaEnv.frex(2:2:20),1))
        
        % overlay mean breath trace on right axis
        yyaxis right
        plot(squeeze(mean(chanDat.tf.breathSeg(useIdx,:), 1, 'omitnan')), ...
            'LineWidth', 2)
        
        title([subjName{s} ' ' condName{c} ' video gamma envelope ITPC'])
        
        if c == 3
            xlabel('Time')
        end
        
    end
end


%% sorted gamma envelope: 


for s = 1:3
    
    chanDat = allDat{s};
    figure('Name', subjName{s}, 'Color', 'w');
    
    for c = 1:3
        
        subplot(3,1,c)
        
        useIdx = chanDat.use == 1 & chanDat.behDat.condition == c;
        [~, order] = sort(chanDat.gammaBurst.t0_idx(useIdx));
        tmp = chanDat.gammaEnv.gamEnv(useIdx,:);
        overallM = mean(tmp, 'all', 'omitnan'); 
        overallS = std(tmp, [], 'all', 'omitnan'); 
        tmp = (tmp - overallM) ./ overallS; 
        imagesc( tmp(order,:))
        set(gca, 'YDir', 'normal')
        clim([-1 3])
        colorbar
        
        % frequency labels
        % yticks(2:2:20)
        % yticklabels(round(chanDat.gammaEnv.frex(2:2:20),1))
        
        % overlay mean breath trace on right axis
        % yyaxis right
        % plot(squeeze(mean(chanDat.tf.breathSeg(useIdx,:), 1, 'omitnan')), ...
        %     'LineWidth', 2)
        
        title([subjName{s} ' ' condName{c} ' video gamma envelope'])
        
        if c == 3
            xlabel('Time')
        end
        
    end
end