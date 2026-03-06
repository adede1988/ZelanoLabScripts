%% make some nifty topos: 

codePre = 'G:\My Drive\GitHub\';
datPre = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\';

%% set paths

addpath(genpath([codePre 'ZelanoLabScripts']))
addpath([codePre 'myFrequentUse'])

%%

datFolder = [datPre 'CHANDAT_processed/scalpPow_itpc']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 


allSubIDs = cell(length(chanFiles),4); 
for ii = 1:length(chanFiles)
    curChan = chanFiles(ii).name;
    nameBits = split(curChan, '_'); 
    L = length(nameBits); 
    if L == 7
        subID = [nameBits{1} '_' nameBits{2} '_' ...
                 nameBits{3} '_' nameBits{4} '_' nameBits{5}];
        taskID = nameBits{6}; 
        sessID = str2num(nameBits{5}); 
        chanID = split(nameBits{7}, '.mat');
        chanID = chanID{1}; 
        type   = nameBits{2}; 
    else
        subID = [nameBits{1} '_' nameBits{2} '_' ...
                 nameBits{3} '_' nameBits{4}];
        taskID = nameBits{5}; 
        sessID = 1; 
        chanID = split(nameBits{6}, '.mat');
        chanID = chanID{1}; 
        type   = nameBits{2}; 
    end
   
    allSubIDs{ii,2} = taskID; 
    allSubIDs{ii,1} = subID;
    allSubIDs{ii,3} = [subID '_' taskID];
    allSubIDs{ii,4} = chanID; 
    allSubIDs{ii,5} = sessID; 
    allSubIDs{ii,6} = type; 


end

[comboKey, ~, comboIdx] = unique(allSubIDs(:,3),  'stable');
allSubIDs(:,7) = arrayfun(@(x) x, comboIdx, 'uniformoutput', false); 

eegLocs = readtable("G:\My Drive\GitHub\ZelanoLabScripts\eegLocs_standard_coords.csv"); 



%% pre specify variables to extract 

% evidence that there is a gamma oscillation: 
allFlatSpec = cell([length(chanFiles), 1]); 


% other stuff: 
allPowBandMax = nan([length(chanFiles), 3, 3, 50]);
allitpcBandMax= nan([length(chanFiles), 3, 3, 50]);
allPowShuf    = nan([length(chanFiles), 3, 3, 50]); %convert values to percentile
allitpcShuf= nan([length(chanFiles), 3, 3, 50]); 

subERP_PAC_peak = nan([length(chanFiles), 4001]);
subERP_PAC_noPeak=nan([length(chanFiles), 4001]); 
subERP_noPAC_peak=nan([length(chanFiles), 4001]); 

subERP_PAC_HRV  = nan([length(chanFiles), 4001]);
subERP_PAC_noHRV  = nan([length(chanFiles), 4001]);

allgamEnv = cell([length(chanFiles), 1]);
allPACgamPeakidx50    = cell([length(chanFiles), 1]); 
allTaskList =  cell([length(chanFiles), 1]); 



allBehDat = cell([length(chanFiles),1]); 

taskERP = nan([3, length(chanFiles), 4001]);
cndList = {'audio', 'focus', 'shadow'};


%% extraction loop: 
parfor start = 1:length(chanFiles)
    disp(start)
    try
        %load: 
        outSum = load(fullfile(chanFiles(start).folder, chanFiles(start).name));
        outSum = outSum.outSum; 



        % Evidence that there is a gamma oscillation: 
        allFlatSpec{start} = outSum.flatSpec; 



        
        allpowBandMax(start,:,:,:)  = outSum.powBandMax; 
        allitpcBandMax(start,:,:,:) = outSum.itpcBandMax; 
        perValP = nan(size(outSum.powBandMax)); 
        perValI = nan(size(outSum.powBandMax)); 
        [nBand, nCnd, nTim] = size(outSum.powBandMax); 
        nSh = size(outSum.powShuf, 1); 
        for b = 1:nBand
            for c = 1:nCnd
                for t = 1:nTim

                    perValP(b,c,t) = squeeze(sum(outSum.powBandMax(b,c,t) > ...
                        outSum.powShuf(:,b,c,t), 1)) / nSh; 
                    perValI(b,c,t) = squeeze(sum(outSum.itpcBandMax(b,c,t) > ...
                        outSum.itpcShuf(:,b,c,t),1)) / nSh; 
                end
            end
        end


    allPowShuf(start,:,:,:) = perValP; 
    allitpcShuf(start,:,:,:) = perValI; 


    

    % get theta power: 
    outSum.behDat.useVec = outSum.useVec; 
    allBehDat{start, 1} = outSum.behDat;


    useVec = outSum.useVec; 

    allGamEnv{start} = outSum.gamEnv(outSum.useVec==1, :); 
    allPACgamPeakidx50{start}    = outSum.behDat.PACgamPeakidx50(outSum.useVec==1);
    allTaskList{start} = string(outSum.behDat.task(outSum.useVec==1)); 


    HRV_RMS  = outSum.behDat.HRV_RMSSD30(useVec);
    HRV_SDNN = outSum.behDat.HRV_SDNN30(useVec);
    HRV_RSA  = outSum.behDat.HRV_RSAamp(useVec); 
    HRV_mm   = outSum.behDat.RR_max_min(useVec); 
    HR       = outSum.behDat.HR_mean(useVec); 

    minmax01 = @(v) local_minmax01(v);
    rms01  = minmax01(HRV_RMS);
    sdnn01 = minmax01(HRV_SDNN);
    rsa01  = minmax01(HRV_RSA);
    
    % --- overall state = mean of the 3 scaled measures ---
    overall = mean([rms01 sdnn01 rsa01], 2, 'omitnan');
    
    
    HRVidx = overall > median(overall); 

    taskVec = string(outSum.behDat.task); 

    %get ERP for PAC HRV+ breaths
    test = outSum.ERP_pacPeak(outSum.useVec, :); 
    subERP_PAC_HRV(start, :) = mean(test(HRVidx,:), 1, 'omitnan');

    %get ERP for PAC HRV- breaths
    test = outSum.ERP_pacPeak(outSum.useVec, :); 
    subERP_PAC_noHRV(start, :) = mean(test(~HRVidx,:), 1, 'omitnan');

    %get ERP for PAC_peak
    useVec = outSum.behDat.goodBreath == 1 & abs(outSum.behDat.PACgamPeakidx50 - outSum.behDat.gamPeakidx50)<=3; 
    useVec = useVec & outSum.useVec; 
    subERP_PAC_peak(start, :) = mean(outSum.ERP_pacPeak(useVec,:), 1, 'omitnan');

    %get ERP for PAC max that is not the breath-wise peak
    useVec = outSum.behDat.goodBreath == 1 & abs(outSum.behDat.PACgamPeakidx50 - outSum.behDat.gamPeakidx50)>3; 
    useVec = useVec & outSum.useVec; 
    subERP_PAC_noPeak(start, :) = mean(outSum.ERP_pacPeak(useVec,:), 1, 'omitnan');
  
    %get ERP for breath-wise max that is not in PAC window
    useVec = outSum.behDat.goodBreath == 1 & abs(outSum.behDat.PACgamPeakidx50 - outSum.behDat.gamPeakidx50)>3; 
    useVec = useVec & outSum.useVec; 
    subERP_noPAC_peak(start, :) = mean(outSum.ERP_allPeak(useVec,:), 1, 'omitnan');
       

    for c = 1:3
        idx = taskVec == cndList{c}; 
        if sum(idx)>5
            taskERP(c, start, :) = mean(outSum.ERP_pacPeak(idx & outSum.useVec, :), 1, 'omitnan');
        end
    end




    catch
        disp(['failure on ' allSubIDs{start, 4}  allSubIDs{start,1} ' ' allSubIDs{start,2} '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    end

end


%% doing things with extracted values: 


%evidence that there is a gamma oscillation:
%get one index into each control (TI and AS only)
conIDX = cellfun(@(chan,type,sess,subID) strcmp('Cz', chan) & ...
                       strcmpi('obe', type) & ...
                       sess==1 &  ...
                       subID~=14, ...
                       allSubIDs(:,4), allSubIDs(:,6), ...
                       allSubIDs(:,5), allSubIDs(:,7));

conFlatSpec = allFlatSpec(conIDX); 
conFlatSpec = cat(1, conFlatSpec{:}); 
conFlatSpec = squeeze(mean(conFlatSpec, 2));
conFlatSpec(isnan(conFlatSpec(:,3)), :) = []; 

figure; 
frex = logspace(log10(.1),log10(200),300);
imagesc(conFlatSpec(:,frex>4))
plotFrex = frex(frex>4); 
xticklabels(plotFrex(20:20:140))
clim([-.5 1.5])









idx = cellfun(@(x) strcmp('P3', x), allSubIDs(:,4));

uniSub = allSubIDs(idx,3); 
types = allSubIDs(idx,6); 
sess  = allSubIDs(idx,5); 
n = length(uniSub); 
idx = find(idx); 
allSubGamTim = nan(n, 4);
prcTilVals = [10 20 30 40];

figure;
hold on 
for ii = 1:length(idx)
    behDat = allBehDat{idx(ii)}; 
    orderedVals = sort(behDat.gamPeakidx50);
    orderedVals = [1; orderedVals(:); 50]; 
    valsUntil = arrayfun(@(x) find(orderedVals>x,1), prcTilVals) ./ ...
        length(orderedVals);
    allSubGamTim(ii,:) = valsUntil; 
    valsUntil = [0 valsUntil 1]; 
    
    if strcmpi(types{ii}, 'obe')
        plot([0 prcTilVals 50] , valsUntil, 'color', 'k', 'linewidth', 2)
    elseif strcmpi(types{ii}, 'dupi')
        if sess{ii} == 1
            plot([0 prcTilVals 50] , valsUntil, 'color', 'red', 'linewidth', 1)
        elseif sess{ii} == 2
            plot([0 prcTilVals 50] , valsUntil, 'color', 'green', 'linewidth', 1)
        end
    end

end





figure;
subplot 131
idx = cellfun(@(x,y,z) strcmp('Cz', x) & ...
                       strcmp('obe', lower(y)) & ...
                       z==1,     allSubIDs(:,4), allSubIDs(:,6), ...
                                 allSubIDs(:,5));

test = allGamEnv(idx);
test2= allPACgamPeakidx50(idx); 
test = cat(1, test{:}); 
test2= cat(1, test2{:}); 
[conidx, order] = sort(test2); 
test = (test - mean(test, 2)) ./ std(test, [], 2); 
imagesc(test(order, :))
clim([-.5 5])
subplot 132
idx = cellfun(@(x,y,z) strcmp('Cz', x) & ...
                       strcmp('dupi', lower(y)) & ...
                       z==1,     allSubIDs(:,4), allSubIDs(:,6), ...
                                 allSubIDs(:,5));

test = allGamEnv(idx);
test2= allPACgamPeakidx50(idx); 
test = cat(1, test{:}); 
test2= cat(1, test2{:}); 
[dupidx1, order] = sort(test2); 
test = (test - mean(test, 2)) ./ std(test, [], 2); 
imagesc(test(order, :))
clim([-.5 5])

subplot 133
idx = cellfun(@(x,y,z) strcmp('Cz', x) & ...
                       strcmp('dupi', lower(y)) & ...
                       z==2,     allSubIDs(:,4), allSubIDs(:,6), ...
                                 allSubIDs(:,5));

test = allGamEnv(idx);
test2= allPACgamPeakidx50(idx); 
test = cat(1, test{:}); 
test2= cat(1, test2{:}); 
[dupidx2, order] = sort(test2); 
test = (test - mean(test, 2)) ./ std(test, [], 2); 
imagesc(test(order, :))
clim([-.5 5])


conVals  = prctile(conidx, [1:99]);
dupVals1 = prctile(dupidx1,[1:99]);
dupVals2 = prctile(dupidx2,[1:99]);


figure;
scatter(conVals, 1:99)
hold on 
scatter(dupVals1, 1:99)
scatter(dupVals2, 1:99)


plot_group_topos(allPowShuf, allSubIDs, eegLocs)

plot_group_topos(allitpcBandMax, allSubIDs, eegLocs)



plotERP(squeeze(taskERP(1,:,:)),squeeze(taskERP(2,:,:)) , eegLocs, allSubIDs, "OBE", 1,9, {'audio', 'focus'})


plotERP(subERP_PAC_peak,subERP_PAC_noPeak , eegLocs, allSubIDs, "OBE", 1,8, {'PAC_peak', 'PAC_noPeak'})


plotERP(subERP_PAC_peak,subERP_noPAC_peak , eegLocs, allSubIDs, "OBE", 1,8, {'PAC_peak', 'Peak_noPAC'})


plotERP(subERP_PAC_peak,subERP_noPAC_peak , eegLocs, allSubIDs, "Dupi", 1, {'PAC_peak', 'Peak_noPAC'})


plotERP(subERP_PAC_HRV,subERP_PAC_noHRV , eegLocs, allSubIDs, "OBE", 1, 14, {'HRV', 'no HRV'})





%% Is there a relationship between timing of gamma peak and timing of theta max power? 




test = cat(1, allBehDat{:}); 
