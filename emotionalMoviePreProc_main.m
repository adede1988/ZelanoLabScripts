

clear
codePre = 'G:\My Drive\GitHub\';
datPre  = {'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ...
         'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\', ...
         'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

sessionIDs = {'251009_OBE_NWU_CP_1', ... %preProcessed
           '250225_OBE_NWU_AS_4', ...  
            '250904_OBE_NWU_TI'};           
datPrei = [2,2,2];       

addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))
addpath(genpath([codePre 'ZelanoLabScripts']))
addpath([codePre 'slowBreathing'])

figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';
% Load once, reuse
EEGLOC = readtable([codePre 'ZelanoLabScripts/myEEGcoords_thetaPhi.csv']);

for s = 1:numel(sessionIDs)
  S.id   = sessionIDs{s};
  S.root = datPre{datPrei(s)};
  S.figPath = figPath; 
  % <<< all subject-specifics here
  [raw, P]      = getSessionParams_emotionTask(S);

  outDat = assemble_outDat_EmotionalMovie(raw, S, P);

  TTL = detect_ttls_EmotionalMovie(raw, P);
  outDat.TTL = TTL; 

  outDat = downsample_data(outDat, P.fs_target);

  if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end


  outDat = preprocess_macros(outDat, P);



  outDat = process_respiration_breathing(outDat, P); 
  outDat.moreThan1 = 1; 
  outDat.rspIDX = P.rspIDX;
  outDat.rspFlip = P.rspFlip; 
  

  outDat = build_behavior_table_EmotionalMovie(outDat);

  R = preprocess_respiration_wholetrace(outDat);
  plot_sniff_epochs(outDat, R);

  if ~isfolder(fullfile(outDat.OGdataDir,  'preProc'))
      mkdir(fullfile(outDat.OGdataDir,  'preProc'))
  end

  save(fullfile(outDat.OGdataDir,  'preProc', ...
                [outDat.sessID '_' outDat.task 'preproc.mat']), ...
                'outDat', "-v7.3")

end

% 
% %% time alignment: 
% 
% idx = cellfun(@(x) contains(x, 'event'), dat.outLabs);
% photoDiode = dat.rawData.trial{1}(idx, :); 
% 
% 
% % photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
% figure; plot(photoDiode)
% % tim = dat.rawdata_emotionClips.time{1}; 
% %TTLs in sample indices 
% TTLs = find(photoDiode(1:length(photoDiode)-1)>-800 &...
%      photoDiode(2:length(photoDiode))<-800);
% diffVals = diff(TTLs); 
% TTLs(diffVals<350) = []; 
% xline(TTLs)
% %task is so simple, there's no need for behavioral alignment: 
% %can just decode from the TTL pulses directly 
% TTLout = struct;
% TTLout.timStamp = zeros(50,1); 
% TTLout.type = zeros(50,1); 
% ii = 1;
% ti = 1; 
% while true
% 
%     curTTL = TTLs(ii); 
%     TTLout.timStamp(ti) = curTTL;
%     if TTLs(ii+2) - curTTL < 1800
%         TTLout.type(ti) = 3; %negative
%         ii = ii + 3; 
%     elseif TTLs(ii+1) - curTTL < 1000
%         TTLout.type(ti) = 2; %positive
%         ii = ii + 2; 
%     else
%         TTLout.type(ti) = 1; %neutral 
%         ii = ii + 1; 
%     end
%     ti = ti + 1; 
%     if ii > 349
%         break
%     end
% 
% end
% fs = dat.rawData.fsample; 
% typeidx = zeros(length(dat.rawData.time{1}),1); 
% for ii = 1:length(TTLout.timStamp)-1
%     TTLout.timStampEnd(ii) = TTLout.timStamp(ii+1) - fs*(1.6);
%     typeidx(TTLout.timStamp(ii):TTLout.timStamp(ii+1)) = TTLout.type(ii); 
% end
% TTLout.timStampEnd(end) = TTLout.timStamp(end)+5*fs; 
% 
% 
% %rereference the 6 macro contacts: 
% 
% idx = cellfun(@(x) contains(x, 'macro'), dat.outLabs);
% ephysDat = dat.rawData.trial{1}(idx,:);
% idx = cellfun(@(x) contains(x, 'rsp'), dat.outLabs);
% rspDat = dat.rawData.trial{1}(idx,:); 
% rspDat = rspDat(1,:); 
% 
% macOut = zeros([5,size(ephysDat,2)]);
% %do bipolar rereferencing 
% for chani = 1:5
%     macOut(chani, :) = squeeze(ephysDat(chani, :) -...
%                                   ephysDat(chani+1,:)); 
% end
% 
% 
% ephysDat = macOut; 
% 
% %get the breath info
% %col 1: onset Y value
% %col 2: onset tim
% %col 3: peak Y value
% %col 4: peak tim
% %col 5: end Y value
% %col 6: end tim
% %col 7: length (end tim - onset tim)
% %col 8: amp (peak Y - avg of two ends)
% %col 9: idx of peak in rspSig2
% %col10: exhale peak Y value
% %col11: exhale peak tim
% %col12: empty
% %col13: empty
% %col14: index
% bmObj = breathTemplates4(rspDat, fs);
% tim = dat.rawData.time{1};
% %col12: video type: 
% for ii = 1:size(bmObj,1)
%     curTim = bmObj(ii, 2); 
%     bmObj(ii, 12) = typeidx(find(curTim <= tim, 1));
% 
% 
% end
% 
% figure
% scatter(bmObj(bmObj(:,12)==1,7), bmObj(bmObj(:,12)==1,8), 'blue', 'filled')
% hold on 
% scatter(bmObj(bmObj(:,12)==2,7), bmObj(bmObj(:,12)==2,8), 'green', 'filled')
% scatter(bmObj(bmObj(:,12)==3,7), bmObj(bmObj(:,12)==3,8), 'red', 'filled')
% 
% %use bulb gamma to create "trials" of gamma peaks
% chanOfInterest = 3; 
% gamFrex = logspace(log10(10), log10(70), 100);
% 
% [phase, pow] = multiphasevec3(gamFrex, ephysDat(chanOfInterest,:), fs, 4,0); 
% 
% %what is the peak gamma frequency: 
% figure
% plot(gamFrex, squeeze(mean(pow, 3)), 'linewidth', 3)
% title('channel 3 power spectrum')
% [~, maxidx] = max(squeeze(mean(pow,3))); 
% xline(gamFrex(maxidx))
% %gamma burst occurences
% figure
% plot(squeeze(pow(:,maxidx,:)))
% gamPow = squeeze(pow(:,maxidx,:)); 
% gamPow = (gamPow - mean(gamPow)) / std(gamPow); 
% gamPow = smoothdata(gamPow, 'gaussian', 100); 
% plot(gamPow)
% yline(3)
% 
% %get gamma burst onsets: 
% gamBursts = struct; 
% gamBursts.onset = find(gamPow(1:length(gamPow)-1) < 3 & ...
%                 gamPow(2:length(gamPow)) > 3);
% 
% %adjust onsets to nearest oscillatory peak: 
% gamBursts.onset = arrayfun(@(x) x - 25 + find(abs(...
%                 atan2(sin(phase(x-25:x+25) - (pi/2)), ...
%                 cos(phase(x-25:x+25) - (pi/2)))) == ...
%             min(abs(atan2(sin(phase(x-25:x+25) - (pi/2)), ...
%                 cos(phase(x-25:x+25) - (pi/2))))) ), ...
%                 gamBursts.onset);
% 
% %trim off the last couple bursts to avoid hitting the end of the recording
% gamBursts.onset(length(gamBursts.onset)-2:length(gamBursts.onset)) = [];
% gamBursts.onset(1:2) = [];
% 
% %gamma burst lengths: 
% gamBursts.lengths = arrayfun(@(x) find(gamPow(x:length(gamPow)-1) >3 & ...
%                    gamPow(x+1:length(gamPow)) <3, 1), ...
%                    gamBursts.onset) ./ fs; 
% 
% gamBursts.onset(gamBursts.lengths>.2) = [];
% gamBursts.lengths(gamBursts.lengths>.2) = [];
% 
% %get the video type that each burst occurred in: 
% gamBursts.type = typeidx(gamBursts.onset); 
% 
% %phase of respiration 
% rspPhase = angle(hilbert(smoothdata(rspDat, 'gaussian', fs)));
% gamBursts.rspPhase = rspPhase(gamBursts.onset); 
% 
% figure
% subplot 311
% histogram(gamBursts.rspPhase(gamBursts.type == 1), 20,...
%     'normalization', 'probability');
% ylim([0, .15])
% subplot 312
% histogram(gamBursts.rspPhase(gamBursts.type == 2), 20,...
%     'normalization', 'probability');
% ylim([0, .15])
% subplot 313
% histogram(gamBursts.rspPhase(gamBursts.type == 3), 20,...
%     'normalization', 'probability');
% ylim([0, .15])
% 
% 
% 
% 
% 
% %plot all breath traces for each condition AND gamma burst timing
% figure
% tmpX = linspace(1,100,25); 
% tmpRes = zeros(25, 3); 
% for ii = 1:size(bmObj,1)
%     if bmObj(ii, 12) > 0
%     if bmObj(ii, 12) == 1
%         subplot 311
%         hold on 
%     elseif bmObj(ii, 12) == 2
%         subplot 312
%         hold on 
%     elseif bmObj(ii, 12) == 3
%         subplot 313
%         hold on 
%     end
%     curstart = find(tim>=bmObj(ii,2),1); 
%     curend = find(tim>=bmObj(ii,6),1); 
%     curBursts = (gamBursts.onset(find(gamBursts.onset>curstart & ...
%                      gamBursts.onset<curend)) - curstart) / fs; 
%     len = (curend - curstart) / fs; 
%     curDat = smoothdata(rspDat(curstart:curend), ...
%                     'gaussian', fs); 
%     curTim = 1/fs : 1/fs : length(curDat)/fs; 
% 
%     newTim = linspace(1/fs, length(curDat)/fs, 100);
%     curBurstIDX = arrayfun(@(x) find(x <= newTim, 1), ...
%                                 curBursts); 
%     interpDat = interp1(curTim, curDat, newTim, 'linear'); 
%     plot(interpDat, 'color', [74/255,103/255,65/255,.1])
%     ylim([-200, 200])
% 
%     tmpidx = arrayfun(@(x) find(x <= tmpX, 1), curBurstIDX); 
% 
%     tmpRes(tmpidx, bmObj(ii, 12)) = tmpRes(tmpidx, bmObj(ii, 12)) + 1; 
% 
%     % scatter(curBurstIDX, zeros(length(curBurstIDX),1)-150, 30, ...
%     %      [220/255,121/255,58/255], 'filled', ...
%     %       'MarkerFaceAlpha', 0.1)
%     end
% 
% end
% subplot 311
% yyaxis right
% plot(tmpX, tmpRes(:,1), 'linewidth', 2, 'color', [220/255,121/255,58/255])
% subplot 312
% yyaxis right
% plot(tmpX, tmpRes(:,2), 'linewidth', 2, 'color', [220/255,121/255,58/255])
% subplot 313
% yyaxis right
% plot(tmpX, tmpRes(:,3), 'linewidth', 2, 'color', [220/255,121/255,58/255])
% 
% 
% 
% 
% %split raw ephys data into trials
% trialDat = zeros(5, fs*2, length(gamBursts.onset));
% trialRsp = zeros(fs*2, length(gamBursts.onset)); 
% for trial = 1:length(gamBursts.onset)
%     curTim = gamBursts.onset(trial); 
%     curTrial = ephysDat(:, curTim-(fs):curTim+(fs)-1);
%     curRsp = rspDat(curTim-(fs):curTim+(fs)-1); 
%     trialDat(:,:,trial) = curTrial; 
%     trialRsp(:,trial) = curRsp; 
% end
% 
% %plot an average gamma burst
% figure
% plot(squeeze(mean(trialDat(3, :, gamBursts.type==1), 3)), ...
%     'color', 'blue', 'linewidth', 2)
% hold on 
% plot(squeeze(mean(trialDat(3, :, gamBursts.type==2), 3)), ...
%     'color', 'green', 'linewidth', 2)
% plot(squeeze(mean(trialDat(3, :, gamBursts.type==3), 3)), ...
%     'color', 'red', 'linewidth', 2)
% 
% figure
% histogram(gamBursts.lengths(gamBursts.type ==1), 'facecolor', 'blue')
% hold on 
% histogram(gamBursts.lengths(gamBursts.type ==2), 'facecolor', 'green')
% histogram(gamBursts.lengths(gamBursts.type ==3), 'facecolor', 'red')
% 
% 
% 
% %get timefrequency decomposition for trialwise data: 
% powDat = cell(5,2); 
% %loop over channels
% for ii = 3
%     ii
%     downDat = resample(squeeze(trialDat(ii,:,:)), 1, 2);
%     highFrex = logspace(log10(2), log10(150), 100);
%     [~, pow] = multiphasevec3(highFrex, downDat', 1000, 6, 0); 
%     pow = permute(pow, [3,1,2]); 
%     powDat{ii,1} = pow; 
%     % getChanMultiTF(...
%         % downDat, highFrex, 1000, .001:.001:2, 5);
%     pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
%     pow = cell2mat(pow); %organize
%     pow = reshape(pow, size(pow,1), size(pow,2)/length(highFrex), []); %organize
% 
% 
%     powDat{ii,2} = pow; 
% 
% 
% end
% 
% %create a 2d matrix of time X bursts X chans/frex
%  allPow = zeros([size(powDat{1,2},[1,2]), size(powDat{1,2},3)*length(powDat)]);
%     for ii = 1:5
%         allPow(:,:,(ii-1)*length(highFrex)+1:ii*length(highFrex)) = powDat{ii,2}; 
%     end
% 
% 
% %calculate the cosine similarity at time points across bursts
% % simMat = zeros([size(allPow,1), size(allPow,2), size(allPow,2)]);
% % 
% % for tt = 1:size(simMat,1)
% % 
% %     %cosine sim of neural representation: 
% %     simMat(tt,:,:) = reshape(cell2mat(arrayfun(@(ii) ...
% %         arrayfun(@(jj) ...
% %             sum(allPow(tt,ii,:) .* allPow(tt,jj,:))/(norm(squeeze(allPow(tt,ii,:))) * norm(squeeze(allPow(tt,jj,:)))), ...
% %         1:size(simMat,2)), ...
% %     1:size(simMat,2), 'uniformoutput', false)), [size(allPow,2), size(allPow,2)]); 
% % 
% % end
% 
% 
% 
% 
% figure('Color', 'w')  % white background
% 
% timeAxis = 800:1200;                % original time vector
% relTime = timeAxis - 1000;          % convert to relative time: -200 to +200
% xTickVals = [800 1000 1200];
% 
% test = sum(pow, 1);
% badIDX = squeeze(sum(test,3)) > 0;
% 
% tmp = powDat{3,2};
% tmpbursts = gamBursts; 
% tmp = tmp(:,100:240, :); 
% tmpbursts.type = tmpbursts.type(100:240); 
% for i = 1:3
%     subplot(1,3,i)
% 
%     % Get average power for condition i
%     plotMat = squeeze(mean(tmp(:, tmpbursts.type==i, :), 2));
% 
%     imagesc(plotMat')  % use relTime on x-axis
%     set(gca, 'YDir', 'normal')
%     caxis([-5, 5])
% 
%     % Set custom ticks and labels
%     % xticks(xTickVals)
%     % xticklabels([-200, 0, 200])
%     % xlim([800 1200])
%     yticks([])
% 
%     % Dashed, thick vertical line at 0
%     xline(1000, 'k--', 'LineWidth', 2)
% 
%     % Title
%     switch i
%         case 1, title('Neutral')
%         case 2, title('Happy')
%         case 3, title('Sad')
%     end
% 
%     % Common labels
%     if i == 1
%         % ylabel('speed of electrical brain waves (lower to higher)')
%     end
%     % xlabel('Time relative to olfactory bulb gamma oscillation activity (milliseconds)')
% end
% 
% 
% 
% figure('Color', 'w')  % white background
% 
% timeAxis = 800:1200;                % original time vector
% relTime = timeAxis - 1000;          % convert to relative time: -200 to +200
% xTickVals = [800 1000 1200];
% 
% tmp = powDat{3,2};
% 
% for i = 1:18
%     subplot(6,3,i)
% 
%     % Get average power for condition i
%     typei = mod(i,3);
%     if typei == 0
%         typei = 3; 
%     end
%     plotMat = squeeze(mean( ...
%         trialDat(floor((i-1)/3)+1,:, gamBursts.type==typei), [1,3]));
% 
%     plot(plotMat)  % use relTime on x-axis
% 
%     yyaxis right
%     plotDat = squeeze(mean(trialRsp(:, gamBursts.type==typei),2)); 
%     plot(plotDat, 'color', 'green')
% 
% 
%     % Dashed, thick vertical line at 0
%     xline(2000, 'k--', 'LineWidth', 2)
% 
%     % Title
%     switch i
%         case 1, title('Neutral')
%         case 2, title('Happy')
%         case 3, title('Sad')
%     end
% 
%     % Common labels
%     if i == 1
%         % ylabel('speed of electrical brain waves (lower to higher)')
%     end
%     % xlabel('Time relative to olfactory bulb gamma oscillation activity (milliseconds)')
% end
% 
% 
% for ii = 1:size(trialDat,3)
%     if rand() < .05
%     figure
%     plot(squeeze(trialDat(3,:,ii)))
%      switch gamBursts.type(ii)
%         case 1, title('Neutral')
%         case 2, title('Happy')
%         case 3, title('Sad')
%     end
%     end
% 
% end
% 
% 
% 
% %stats: 
% 
% %stack the power data from all channels: 
% 
% allPow = cat(4, powDat{:,2});
% nt = size(allPow,2); 
% nChan = size(allPow, 4); 
% chanCode = arrayfun(@(x) ones(nt,1)*x, 1:nChan, 'uniformoutput', false);
% chanCode = cat(1, chanCode{:});
% cndCode = repmat(gamBursts.type, [nChan,1]); 
% 
% tStatHappy = zeros(size(allPow, [1,2])); 
% tStatSad = tStatHappy; 
% 
% %loop time points
% for tt = 1:size(allPow,1)
%     tt
%     %loop frequencies
%     for fi = 1:length(highFrex)
%             curdat = squeeze(allPow(tt, :, fi, :));
%             curdat = curdat(:); 
% 
%             modDat = table(curdat, categorical(cndCode), ...
%                 categorical(chanCode), ...
%                 'VariableNames', {'pow', 'cnd', 'chan'}); 
% 
%             modDat(modDat.cnd == "0", :) = []; 
%             modDat.cnd = removecats(modDat.cnd);
%             lme = fitlme(modDat, ...
%                 'pow ~ cnd + (1|chan) ');
%             tStatHappy(tt,fi) = lme.Coefficients.tStat(2); 
%             tStatSad(tt,fi) = lme.Coefficients.tStat(3); 
% 
% 
%     end
% end
% 
% tStatHappy = zeros(size(allPow, [1,2])); 
% tStatSad = tStatHappy; 
% 
% % Sample every 10th time point
% timeIdx = 1:10:size(allPow, 1);
% 
% % Loop over selected time points
% for ttIdx = 1:length(timeIdx)
%     tt = timeIdx(ttIdx);
%     fprintf('Processing time point %d\n', tt);
% 
%     % Loop over frequencies
%     for fi = 1:length(highFrex)
%         curdat = squeeze(allPow(tt, :, fi, :));
%         curdat = curdat(:); 
% 
%         modDat = table(curdat, categorical(cndCode), ...
%                        categorical(chanCode), ...
%                        'VariableNames', {'pow', 'cnd', 'chan'}); 
% 
%         modDat(modDat.cnd == "0", :) = []; 
%         modDat.cnd = removecats(modDat.cnd);
% 
%         lme = fitlme(modDat, 'pow ~ cnd + (1|chan)');
% 
%         % Store results only at sampled time points
%         tStatHappy(tt, fi) = lme.Coefficients.tStat(2); 
%         tStatSad(tt, fi) = lme.Coefficients.tStat(3); 
%     end
% end
% 
% happy = tStatHappy(timeIdx, 1:100);
% sad = tStatSad(timeIdx, 1:100);
% 
% 
% figure
% subplot 231
% tmp = mean(allPow, 4); 
% plotMat = squeeze(mean(tmp(:, gamBursts.type==1, :), 2)); 
% imagesc(plotMat')
% set(gca, 'ydir', 'normal')
% caxis([-5, 5])
% xlim([200, 800])
% yticks(1:10:100)
% yticklabels(round(highFrex(1:10:100), 2) )
% xline(500)
% yyaxis right 
% curRsp = trialRsp(:,gamBursts.type ==1); 
% curRsp = mean(curRsp, 2); 
% curRsp = resample(curRsp, 1, 2); 
% plot(curRsp, 'color', 'red')
% title('neutral')
% 
% subplot 232
% plotMat = squeeze(mean(tmp(:, gamBursts.type==2, :), 2)); 
% imagesc(plotMat')
% set(gca, 'ydir', 'normal')
% caxis([-5, 5])
% xlim([200, 800])
% yticks(1:10:100)
% yticklabels(round(highFrex(1:10:100), 2) )
% xline(500)
% yyaxis right 
% curRsp = trialRsp(:,gamBursts.type ==2); 
% curRsp = mean(curRsp, 2); 
% curRsp = resample(curRsp, 1, 2); 
% plot(curRsp, 'color', 'red')
% title('happy')
% 
% subplot 233
% plotMat = squeeze(mean(tmp(:, gamBursts.type==3, :), 2)); 
% imagesc(plotMat')
% set(gca, 'ydir', 'normal')
% caxis([-5, 5])
% xlim([200, 800])
% yticks(1:10:100)
% yticklabels(round(highFrex(1:10:100), 2) )
% xline(500)
% yyaxis right 
% curRsp = trialRsp(:,gamBursts.type ==3); 
% curRsp = mean(curRsp, 2); 
% curRsp = resample(curRsp, 1, 2); 
% plot(curRsp, 'color', 'red')
% title('sad')
% 
% subplot 235
% 
% imagesc(abs(happy'))
% set(gca, 'ydir', 'normal')
% caxis([3, 6])
% xlim([20, 80])
% yticks(1:10:100)
% yticklabels(round(highFrex(1:10:100), 2) )
% xline(50)
% title('happy - neutral t stat')
% 
% subplot 236
% 
% imagesc(abs(sad'))
% set(gca, 'ydir', 'normal')
% caxis([3, 6])
% xlim([20, 80])
% yticks(1:10:100)
% yticklabels(round(highFrex(1:10:100), 2) )
% xline(50)
% title('sad - neutral t stat')
% 
% 
% 
% 
% 
% numTime = size(allPow, 1);
% numTrials = size(allPow, 2);
% numFeatures = size(allPow, 3);
% 
% timeIndices = 1:20:numTime;
% simMat = zeros(length(timeIndices), numTrials, numTrials);  % reduced size
% 
% for kk = 1:length(timeIndices)
%     tt = timeIndices(kk);
%     X = squeeze(allPow(tt, :, :));  % [trials x features]
% 
%     % Normalize each row to unit norm
%     Xnorm = sqrt(sum(X.^2, 2));
%     Xnorm(Xnorm == 0) = eps;  % avoid division by zero
%     Xunit = X ./ Xnorm;
% 
%     % Cosine similarity via dot product
%     simMat(kk, :, :) = Xunit * Xunit';
% end
% 
% 
% 
% 
% numTrials = size(simMat, 2);
% sameMask = false(numTrials, numTrials);
% diffMask = false(numTrials, numTrials);
% 
% for i = 1:numTrials
%     for j = 1:numTrials
%         if i ~= j
%             if gamBursts.type(i) == gamBursts.type(j)
%                 sameMask(i, j) = true;
%             else
%                 diffMask(i, j) = true;
%             end
%         end
%     end
% end
% 
% 
% numTimePoints = size(simMat, 1);
% withinCondSim = zeros(1, numTimePoints);
% betweenCondSim = zeros(1, numTimePoints);
% 
% for tt = 1:numTimePoints
%     simSlice = squeeze(simMat(tt, :, :));
% 
%     % Extract upper triangle without diagonal (optional: to avoid duplicate pairs)
%     simSliceTri = triu(simSlice, 1);
% 
%     withinCondSim(tt) = mean(simSliceTri(sameMask), 'omitnan');
%     betweenCondSim(tt) = mean(simSliceTri(diffMask), 'omitnan');
% end
% 
% figure
% t = linspace(-0.5, 0.5, numTimePoints);  % assuming 1 sec window
% plot(t, withinCondSim, 'b', 'LineWidth', 2); hold on;
% plot(t, betweenCondSim, 'r', 'LineWidth', 2);
% yline(0, '--k');
% legend('Within Condition', 'Between Condition');
% xlabel('Time (s) relative to burst onset');
% ylabel('Mean Cosine Similarity');
% title('Neural similarity within vs. between condition');
% 
% 
% 
% 
% 
% [~, order] = sort(TTLout.type); 
% 
% figure
% imagesc(squeeze(simMat(1,order, order)))
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %chans X time X trials
% trialDat = zeros(6, fs*2, length(TTLout.timStamp));
% trialRsp = zeros(fs*2, length(TTLout.timStamp)); 
% for trial = 1:length(TTLout.timStamp)
%     curTim = TTLout.timStampEnd(trial); 
%     curTrial = ephysDat(1:6, curTim:curTim+fs*2-1);
%     curRsp = ephysDat(8,curTim:curTim+fs*2-1); 
%     trialDat(:,:,trial) = curTrial; 
%     trialRsp(:,trial) = curRsp; 
% end
% 
% 
% powDat = cell(6,2); 
% %loop over channels
% for ii = 1:6
%     ii
%     downDat = resample(squeeze(trialDat(ii,:,:)), 1, 2);
%     highFrex = logspace(log10(2), log10(150), 100);
%     [~, pow] = multiphasevec3(highFrex, downDat', 1000, 6, 0); 
%     pow = permute(pow, [3,1,2]); 
%     powDat{ii,1} = pow; 
%     % getChanMultiTF(...
%         % downDat, highFrex, 1000, .001:.001:2, 5);
%     pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
%     pow = cell2mat(pow); %organize
%     pow = reshape(pow, size(pow,1), size(pow,2)/length(highFrex), []); %organize
% 
% 
%     powDat{ii,2} = pow; 
% 
% 
% end
% 
%     allPow = zeros([size(powDat{1},[1,2]), size(powDat{1},3)*length(powDat)]);
%     for ii = 1:6
%         allPow(:,:,(ii-1)*length(highFrex)+1:ii*length(highFrex)) = powDat{ii}; 
%     end
% 
%     %get the cosine similarity matrix for time X trials X trials 
% 
%     % t1 * t2 / magnitude(t1)*magnitude(t2)
% 
%     %loop on time
%     allPowSave = allPow; 
%     allPow = mean(allPow, 1); 
%     simMat = zeros([size(allPow,1), size(allPow,2), size(allPow,2)]);
%     simMatStim = simMat; 
%     for tt = 1:size(simMat,1)
% 
%         %cosine sim of neural representation: 
%         simMat(tt,:,:) = reshape(cell2mat(arrayfun(@(ii) ...
%             arrayfun(@(jj) ...
%                 sum(allPow(tt,ii,:) .* allPow(tt,jj,:))/(norm(squeeze(allPow(tt,ii,:))) * norm(squeeze(allPow(tt,jj,:)))), ...
%             1:size(simMat,2)), ...
%         1:size(simMat,2), 'uniformoutput', false)), [size(allPow,2), size(allPow,2)]); 
% %         [~, order] = sort(subDat.ToRem); 
%         %angle dif sim : 
%         % simMatStim(tt,:,:) = reshape(cell2mat(arrayfun(@(ii) ...
%         %     arrayfun(@(jj) ...
%         %         pi - abs(angdiff(subDat.ToRem(ii)*2/180*pi , subDat.ToRem(jj)*2/180*pi)), ...
%         %     1:size(simMat,2)), ...
%         % 1:size(simMat,2), 'uniformoutput', false)), [size(allPow,2), size(allPow,2)]); 
% 
%     end
% 
% 
%     %plot all trials with respiration: 
%     for ii = 1:175
%         ii
%         figure
%         imagesc(squeeze(allPowSave(:,ii, :))')
%         caxis([-10 100])
%         yyaxis right
%         plot(normalize(resample(trialRsp(:,ii), 1, 2)), 'color', 'green')
%         hold on 
%         plot(normalize(resample(squeeze(trialDat(3,:,ii)), 1, 2)), ...
%             'color', 'red', 'linestyle', '-', 'linewidth', 2)
%         title([num2str(ii) 'type: ' num2str(trialTypes(ii))])
%     end
% 
%     [~, order] = sort(TTLout.type); 
%     figure
%     plot(squeeze(mean(allPow(1, TTLout.type==1, :), 2)))
%     hold on 
% 
%     plot(squeeze(mean(allPow(1, TTLout.type==2, :), 2)))
% 
%     plot(squeeze(mean(allPow(1, TTLout.type==3, :), 2)))
%     yline(0)
%     xline([100.5 200.5 300.5 400.5 500.5])
%     legend({'normal', 'happy', 'sad'})
% 
% 
% 
%     figure
%     imagesc(squeeze(simMat(1,order, order)))
% 
% 
%     simMat = squeeze(simMat); 
%     withinSim = [];   % to store within-type similarities
%     betweenSim = [];  % to store between-type similarities
%     trialTypes = TTLout.type; 
%     numTrials = length(trialTypes);
% 
%     for i = 1:numTrials
%         for j = i+1:numTrials  % avoid duplicates and diagonal
%             if trialTypes(i) == trialTypes(j)
%                 withinSim(end+1) = simMat(i, j);
%             else
%                 betweenSim(end+1) = simMat(i, j);
%             end
%         end
%     end
% 
% % Optional summary statistics
% meanWithin = mean(withinSim);
% meanBetween = mean(betweenSim);
% 
% figure;
% histogram(withinSim, [0:.1:1], 'normalization', 'probability')
% hold on 
% histogram(betweenSim, [0:.1:1], 'normalization', 'probability')
% legend({'within sim', 'between sim'})
% 
% 
%     %look at neural similarity as a function of stimulus similarity across
%     %time, z-scored against a shuffled null
%     neuToStimSim = zeros(size(simMat,1),1); 
%     neuToStimSimNULL = zeros(size(simMat,1), 1000); 
%     parfor tt = 1:size(simMat,1)
% 
%         neu = squeeze(simMat(tt,:,:)); 
%         stim = squeeze(simMatStim(tt, :,:)); 
% 
% 
%         mask = true(size(neu)); 
%         mask = triu(mask,1); 
%         mask = mask(:); 
% 
%         neu = neu(:); 
%         stim = stim(:); 
% 
%         neu(mask==0) = []; 
%         stim(mask==0) = [];
% 
%         neuToStimSim(tt) = corr(neu, stim, 'type', 'Spearman'); 
% 
%         for shuf = 1:100
%             stim = randsample(stim, length(stim), false); 
%             neuToStimSimNULL(tt,shuf) = corr(neu, stim, 'type', 'Spearman'); 
%         end
% 
%     end
% 
%     test = (neuToStimSim - mean(neuToStimSimNULL,2)) ./ std(neuToStimSimNULL,[],2); 
% 
%     allSubDat(chan).ZneuToStimSim = test; 
%     allSubDat(chan).simMat = simMat; 
%     allSubDat(chan).simMatStim = simMatStim;
%     allSubDat(chan).Error = subDat.Error; 
%     allSubDat(chan).toRem = subDat.ToRem; 
%     allSubDat(chan).mulTim = mulTim; 
% 
% 
% 
%     disp(['sub: ' num2str(chan) ' of ' num2str(length(chanFiles)) '; time: ' num2str(round(toc))])
%     figure
%     plot(mulTim, test)
%     title(chan)
% 
% end
