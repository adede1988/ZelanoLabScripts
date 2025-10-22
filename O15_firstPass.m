clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ... 
           'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\',...
           'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

%prefix index for data folder: 
datPrei = [1,1,1,2,2,1,1]; 

sessionIDs = {'250818_Dupi_NMH_JH_1', ... %has behavior
               '250623_DUPI_NMH_KS_2',... %has behavior
               '250623_Dupi_NMH_KS_1',... %has behavior
               '250908_OBE_NWU_AS', ...   %has behavior
                '250904_OBE_NWU_TI', ...  %has behavior
                '250818_Dupi_NMH_JH_2',...
                '250811_Dupi_NMH_TPB_1'};   %has behavior

%there are multiple respiration channels in many recordings
%which one is right for each session: 
rspIDX = [1,1,1,1,1,1,1]; 
rspFlip = [1,1,1,1,1,1,1]; %hard code flip

addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
addpath([codePre 'ZelanoLabScripts'])
addpath(genpath([codePre 'myFrequentUse']))
addpath([codePre 'myFrequentUse/export_fig_repo'])
addpath(genpath("C:\Users\dtf8829\Documents\eeglab2025.0.0"))

addpath([codePre 'fieldtrip-20230118'])
addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])

saveFigs = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\OBE_FIGURES'; 

set(0, 'defaultfigurewindowstyle', 'normal')
ft_defaults

for sessi = 4:length(sessionIDs)

    outDat = load([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_O15preproc.mat']).outDat; 


    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.trialDat(idx, :, :); 
    
      
        
    idx = rspIDX(sessi); 
    rspDat = squeeze(rspDat(idx, :, :)); 


    idx = cellfun(@(x) contains(x, 'macBP'), outDat.labels);
    macOut = outDat.trialDat(idx, :, :); 
    macLabs = outDat.labels(idx); 
    
  

    [pow, f] = trial_psd(macOut, outDat.fs); 

   
    [badIDX] = covMatClean_breathing(macOut, ones(size(macOut,3),1), ...
                                             ones(size(macOut,3),1)*3001);

    % for ii = 1:length(badIDX)
    %     figure
    % 
    %     hold on 
    %     for chani = 1:5
    %         plot(macOut(chani,:,badIDX(ii)) + chani*100)
    %     end
    %     title(ii)
    % 
    % end

    macOut(:, :, badIDX) = []; 
    tmpBeh =  outDat.sniffInfo;
    tmpRsp = rspDat; 
    tmpRsp(:, badIDX) = []; 
    tmpBeh(badIDX, :) = []; 
    outDat.badIDX = badIDX; 

    




    frex = logspace(log10(.1),log10(200),200);
    numfrex = length(frex); 
    stds = logspace(log10(3),log10(10),numfrex);
    figure('position', [0,0,1000, 1000])
    cndNames = {'start', 'free', 'confirm'};
    
    for chani = 1:5
        chani
        if ~exist(['a' datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_O15_macPowPhase_' num2str(chani)...
                    '.mat'], 'file')
        pow = getChanTrialTF(squeeze(macOut(chani,:,:)), frex, ...
            numfrex, stds, outDat.fs);
        phase = angle(pow); 
        pow = abs(pow).^2; 
        powz = arrayfun(@(x) myChanZscore(pow(:,:,x)), ...
                        1:size(pow,3), 'UniformOutput',false ); %z-score
        powz = cell2mat(powz); %organize
        powz = reshape(powz, size(powz,1), size(powz,2)/numfrex, []); 

       

        pow = powz; 

        powOut = struct; 
        powOut.channel = macLabs{chani}; 
        powOut.sessID = outDat.sessID; 
        powOut.sniffInfo = tmpBeh; 
        powOut.rspDat = tmpRsp;
        powOut.tim = outDat.tim; 
        powOut.fs = outDat.fs; 
        powOut.behDat = outDat.behDat; 
        powOut.covCleaningDone = 1; 
        powOut.pow = pow; 
        powOut.phase = phase; 
        powOut.frex = frex; 
        powOut.stds = stds; 
        idx = cellfun(@(x) contains(x, 'blinkIndicator'), outDat.labels);
        powOut.blinkDat = outDat.trialDat(idx,:,:); 
        powOut.blinkDat(:, :, badIDX) = []; 

        save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                sessionIDs{sessi} '_O15_macPowPhase_' num2str(chani)...
                '.mat'], ...
                'powOut', "-v7.3")

        else
        powOut = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                '\preProc\' sessionIDs{sessi} '_O15_macPowPhase_' ...
                num2str(chani) '.mat']).powOut;
        % blinkMask = squeeze(outDat.trialDat(52,:,:));
        % blinkMask = blinkMask==2;
        % blinkMask(:,badIDX) = []; 
        % [Xlock, onsetSamples, onsetTrials] = ...
        %     lock_power_to_blinks(blinkMask, powz);
        % figure
        % imagesc(squeeze(median(Xlock, 2))')
        % set(gca, 'ydir', 'normal')
        % yticks(10:20:numfrex)
        % yticklabels(round(frex(10:20:numfrex)))
        end
    for cnd = 1:3
        subplot(5,3,(chani-1)*3+cnd)
        plotMat = squeeze(mean(powOut.pow(:, ...
            powOut.sniffInfo(:,6)==cnd, :), 2, 'omitnan'))';
        tStart = powOut.tim(1);
        tEnd   = powOut.tim(end);
        
        imagesc([powOut.tim(1) powOut.tim(end)], ...
        [1 size(plotMat,1)], ...
        plotMat)
        set(gca, 'ydir', 'normal')

        yticks(10:20:numfrex)
        yticklabels(round(frex(10:20:numfrex)))
        caxis([-5, 10])
        
        yyaxis right
        plot(powOut.tim, mean(powOut.rspDat(:, ...
                            powOut.sniffInfo(:,6)==cnd), 2), ...
                             'color', 'green', 'linewidth', 2)
       
        title([sessionIDs{sessi} ' ' cndNames{cnd}],...
            'interpreter', 'none')
    end

        


    end
    export_fig([saveFigs '/' sessionIDs{sessi} '_O15_sniff_tf_clean' '.jpg'], '-r300')


    

end


  % idx = cellfun(@(x) contains(x, 'Fp'), outDat.labels); 
    % blinkDat = squeeze(outDat.trialDat(idx, :, :)); 
    % for triali = 1:size(macroDat,3)
    %     figure
    %     hold on 
    %     for chani = 1:6
    %         plot(macroDat(chani,:,triali) + chani*100)
    % 
    %     end
    %     plot(blinkDat(:,triali), 'color', 'k', 'linewidth', 2)
    %     yyaxis right
    %     plot(rspDat(:,triali), 'color', 'green', 'linewidth', 2)
    %     title(triali)
    % end


    % macOut = zeros([5, size(macroDat, [2,3])]);
    % %do bipolar rereferencing 
    % for chani = 1:5
    %     macOut(chani, :, :) = squeeze(macroDat(chani, :, :) -...
    %                                   macroDat(chani+1,:,:)); 
    % end
    % 
 % switch sessionIDs{sessi}
    %     case '250904_OBE_NWU_TI'
    %         spikeThresh = 50;
    %         spikeWin = 11; 
    %         knownIC = 1; 
    %     case '250623_DUPI_NMH_KS_2'
    %         spikeThresh = 15;
    %         spikeWin = 11; 
    %         knownIC = 4; 
    %     case '250623_Dupi_NMH_KS_1'
    %         spikeThresh = 15;
    %         spikeWin = 11; 
    %         knownIC = 10; 
    %     case '250818_Dupi_NMH_JH_2'
    %         spikeThresh = 20;
    %         spikeWin = 11; 
    %         knownIC = 4; 
    %     case '250818_Dupi_NMH_JH_1'
    %         spikeThresh = 20;
    %         spikeWin = 11; 
    %         knownIC = 1; 
    %     case '250908_OBE_NWU_AS'
    %         spikeThresh = 20;
    %         spikeWin = 11; 
    %         knownIC = 10; 
    %     otherwise
    %         spikeThresh = 20;
    %         spikeWin = 11; 
    %         knownIC = []; 
    % end

    % covariance matrix based cleaning doesn't really do a good job here
 % % test = sum(powz(:,:,frex>50), 3);
        % % 
        % % test2 = squeeze(macOut(chani,:,:)); 
        % % test2(test>prctile(test(:), 95)) = NaN; 
        % % 
        % % % for triali = 1:size(powz,2)
        % % %     figure
        % % %     plot(squeeze(macOut(chani,:,triali)))
        % % %     hold on 
        % % %     plot(squeeze(test2(:,triali)))
        % % %     title(triali)
        % % % 
        % % % end
        % % 
        % % for fi = 1:numfrex
        % %     tmp = squeeze(pow(:,:,fi));
        % %     tmp(test>prctile(test(:), 95)) = NaN;
        % %     pow(:,:,fi) = tmp; 
        % % end
        % % powz2 = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
        % % powz2 = cell2mat(powz2); %organize
        % % powz2 = reshape(powz2, size(powz2,1), size(powz2,2)/numfrex, []); %organize

        % %remove bad trials for this channel individually 
        % test = permute(pow, [3,1,2]); 
        % badIDX = covMatClean_breathing(test, ones(size(macOut,3),1), ...
        %                                      ones(size(macOut,3),1)*3001);
        % 
        % pow(:, badIDX, :) = []; 
        % tmpRsp(:, badIDX) = []; 
        % tmpBeh(badIDX, :) = []; 
       
        
        % meanVals = squeeze(mean(powz, 2)); 
        % 
        % 
        % %%%%%%% targeted cleaning: 
        % [timeIdx, freqIdx] = find(meanVals > 5);
        % badTrials =[]; 
        % for k = 1:numel(timeIdx)
        %     t = timeIdx(k);
        %     f = freqIdx(k);
        % 
        %     % Extract values across trials at this bin
        %     vals = squeeze(powz(t, :, f)); 
        %     thresh = prctile(squeeze(powz(t, :, f)), 95);   % 56 x 1
        %     vals(vals>thresh) = []; 
        %     meanVal = mean(vals); 
        %     sdVal = std(vals); 
        %     % Step 2: compute within-bin z-scores across trials
        %     zvals = (squeeze(powz(t, :, f)) - meanVal) ./ sdVal;
        % 
        %     % Step 3: flag outlier trials
        %     badTrials = [badTrials, find(zvals > 5)];  % threshold can be tuned
        % 
        % 
        % end
        % trialCounts = histcounts(badTrials, 0.5:1:(size(powz,2)+0.5));
        % badIdx = find(trialCounts / numel(timeIdx) > .1); 
        
        % pow(:, badIDX, :) = []; 
        % tmpRsp(:, badIDX) = []; 
        % tmpBeh(badIDX, :) = []; 
        % pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
        % pow = cell2mat(pow); %organize
        % pow = reshape(pow, size(pow,1), size(pow,2)/numfrex, []); %organize
        % 
         
                

        %%%%%%%%%%%
    % [b,a] = butter(4, [5,150]/(outDat.fs/2), 'bandpass');
    % gammaSig = arrayfun(@(x) ...
    %         hilbert(filtfilt(b,a, squeeze(macOut(x,:,:)))),...
    %             1:5, 'uniformoutput', false); 
    % gammaSig = cat(3, gammaSig{:}); 
    % gammaSig = real(permute(gammaSig, [3,1,2])); 
    % [test, prominence] = detect_spikes(macOut,spikeThresh, spikeWin,...
    %     false, gammaSig); 
    % 
    % out = ica_flag_spikes(macOut, test, prominence, 'Fs', 500, ...
    %     'knownIC', knownIC);
    % 
    % 
    % 


    % for ii = 1:20
    % % ii = 54
    % 
    %     figure
    % 
    %     tmp = squeeze(real(macOut(2,:,ii)));
    %     plot(outDat.tim, tmp)
    %     hold on 
    %     plot(outDat.tim,out.data_clean(2,:,ii))
    % 
    % 
    % 
    % end

   % macOut = out.data_clean; 
    % covariance matrix based cleaning doesn't really do a good job here
    % [badIDX] = covMatClean_breathing(macOut, ones(size(macOut,3),1), ...
    %                                          ones(size(macOut,3),1)*3001);
    % 
    % for bi = 1:length(badIDX)
    %     figure
    %     hold on 
    %     for ci = 1:5
    % 
    %         plot(macOut(1,:,badIDX(bi))+ ci*100)
    % 
    %     end
    % 
    %     title(badIDX(bi))
    % end
    % macOut(:,:,badIDX) = []; 
    % outDat.sniffInfo(badIDX, :) = []; 
    % rspDat(:, badIDX) = []; 


% 
% 
% 
% 
% for ii = 1:82
% figure
% plot(squeeze(macOut(1,:,ii)))
% title(ii)
% end
% 
% 
% 
% 
% 
% for ii = 1:79
% figure
% imagesc(squeeze(pow(:,ii,:)))
% title(ii)
% caxis([-20, 50])
% end



% % [b,a] = butter(4, [5,150]/(outDat.fs/2), 'bandpass');
    % % gammaSig = arrayfun(@(x) ...
    % %         hilbert(filtfilt(b,a, squeeze(macOut(x,:,:)))),...
    % %             1:5, 'uniformoutput', false); 
    % % gammaSig = cat(3, gammaSig{:}); 
    % % gammaSig = real(permute(gammaSig, [3,1,2])); 
    % % % limVals = prctile(gammaSig, [5 95], [2,3]); 
    % % % meanVal = mean(gammaSig(gammaSig > limVals(1) & ...
    % % %                         gammaSig < limVals(2)), "all");
    % % % sdVal = std(gammaSig(gammaSig > limVals(1) & ...
    % % %                         gammaSig < limVals(2)),[], "all");
    % % % gammaSig = (gammaSig - meanVal) ./ ...
    % % %                         sdVal; 
    % % [b,a] = butter(4, [35,150]/(outDat.fs/2), 'bandpass');
    % % highSig = arrayfun(@(x) ...
    % %         hilbert(filtfilt(b,a, squeeze(macOut(x,:,:)))),...
    % %             1:5, 'uniformoutput', false); 
    % % highSig = cat(3, highSig{:}); 
    % % highSig = abs(permute(highSig, [3,1,2]).^2); 
    % % 
    % % % limVals = prctile(highSig(:), [5 95]); 
    % % % meanVal = mean(highSig(highSig > limVals(1) & ...
    % % %                         highSig < limVals(2)), "all");
    % % % sdVal = std(highSig(highSig > limVals(1) & ...
    % % %                         highSig < limVals(2)),[], "all");
    % % % highSig = (highSig - meanVal) ./ ...
    % % %                         sdVal; 
    % % 
    % % 
    % % 
    % % [test, ~, surroundVals, surroundVals2] = detect_spikes(macOut, 20, 41, false, gammaSig); 

    % test = detect_spikes_refined(macOut, 20, ...
    % 'Fs',2000, 'UseGammaMask',true, 'GammaBand',[40 90], ...
    % 'MinCycles',3, 'EnvThreshK',2.5, 'ZeroCrossMin',4);