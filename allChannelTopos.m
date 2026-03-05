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


allPowBandMax = nan([length(chanFiles), 3, 3, 50]);
allitpcBandMax= nan([length(chanFiles), 3, 3, 50]);
allPowShuf    = nan([length(chanFiles), 3, 3, 50]); %convert values to percentile
allitpcShuf= nan([length(chanFiles), 3, 3, 50]); 

subERP_PAC_peak = nan([length(chanFiles), 4001]);
subERP_PAC_noPeak=nan([length(chanFiles), 4001]); 
subERP_noPAC_peak=nan([length(chanFiles), 4001]); 


parfor start = 1:length(chanFiles)
    disp(start)
    try
        outSum = load(fullfile(chanFiles(start).folder, chanFiles(start).name));
        outSum = outSum.outSum; 
        
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



    %get ERP for PAC_peak
    useVec = outSum.behDat.goodBreath == 1 & abs(outSum.behDat.PACgamPeakidx50 - outSum.behDat.gamPeakidx50)<=3; 
    subERP_PAC_peak(start, :) = mean(outSum.ERP_pacPeak(useVec,:), 1, 'omitnan'); 

    %get ERP for PAC max that is not the breath-wise peak
    useVec = outSum.behDat.goodBreath == 1 & abs(outSum.behDat.PACgamPeakidx50 - outSum.behDat.gamPeakidx50)>3; 
    subERP_PAC_noPeak(start, :) = mean(outSum.ERP_pacPeak(useVec,:), 1, 'omitnan'); 
  
    %get ERP for breath-wise max that is not in PAC window
    useVec = outSum.behDat.goodBreath == 1 & abs(outSum.behDat.PACgamPeakidx50 - outSum.behDat.gamPeakidx50)>3; 
    subERP_noPAC_peak(start, :) = mean(outSum.ERP_allPeak(useVec,:), 1, 'omitnan'); 
       
    catch
        disp(['failure on ' allSubIDs{start, 4}  allSubIDs{start,1} ' ' allSubIDs{start,2} '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    end

end


plot_group_topos(allPowShuf, allSubIDs, eegLocs)

plot_group_topos(allitpcBandMax, allSubIDs, eegLocs)



plotERP(subERP_PAC_peak, subERP_PAC_noPeak, eegLocs, allSubIDs, "Dupi", 2,6, {'peak', 'noPeak'})


