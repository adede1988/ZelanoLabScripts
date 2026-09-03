function run_scores()
% RUN_SCORES  One sweep over all Dupi+OBE finals: behavioral olfactory scores
% (cue/thresh/O15) and macBP 25-58 Hz gamma FOOOF selection (all tasks).
% Loads each final once. Writes long-format CSVs incrementally.

here = fileparts(mfilename('fullpath'));
proj = fileparts(here); repo = fileparts(proj);
addpath(repo); addpath(here);
outDir = fullfile(proj,'out','tables'); if ~exist(outDir,'dir'), mkdir(outDir); end
H = lib_scores();
BAND = [25 58];

idx = readtable(fullfile(outDir,'session_index.csv'), 'TextType','string');
% analysis set: Dupi + OBE only (EEG excluded), on disk
keep = (idx.cohort=="Dupi" | idx.cohort=="OBE") & idx.onDisk==1;
idx = idx(keep,:);
fprintf('Scoring %d finals (Dupi+OBE, on disk)\n', height(idx));

behFile = fullfile(outDir,'behavioral_long.csv');
macFile = fullfile(outDir,'macbp_gamma_long.csv');
blkFile = fullfile(outDir,'breathing_blocks.csv');
bestFile= fullfile(outDir,'macbp_best.csv');

% resume: if macbp_best.csv already has rows, keep them and skip done (sessID,task)
done = containers.Map('KeyType','char','ValueType','logical');
resume = exist(bestFile,'file')==2;
if resume
    try
        Bt = readtable(bestFile,'TextType','string');
        for i=1:height(Bt), done([char(Bt.sessID(i)) '|' char(Bt.task(i))])=true; end
    catch, resume=false; end
end
if resume && done.Count>0
    fprintf('RESUME: %d units already done; processing remaining\n', done.Count);
else
    resume=false;
    fid=fopen(behFile,'w'); fprintf(fid,'sessID,task,cohort,group,participant,sessNum,metric,value,nTrials\n'); fclose(fid);
    fid=fopen(macFile,'w'); fprintf(fid,'sessID,task,cohort,group,participant,sessNum,nMacBP,chanIdx,chanLabel,peakFlatDb,gammaDetected,peakFreq,peakHeight,apExponent,spikeFrac,isBest\n'); fclose(fid);
    fid=fopen(bestFile,'w');fprintf(fid,'sessID,task,cohort,group,participant,sessNum,nMacBP,bestChan,bestLabel,bestPeakFlatDb,bestGammaDetected,bestSpikeFrac,nNonMaxDetected,nMacBPDetected\n'); fclose(fid);
    fidBlk = fopen(blkFile,'w'); fprintf(fidBlk,'sessID,cohort,group,blockTask,condition,noseMouth,warp,nBreaths\n'); fclose(fidBlk);
end

for r = 1:height(idx)
    id=char(idx.sessID(r)); task=char(idx.task(r)); coh=char(idx.cohort(r));
    grp=char(idx.group(r)); part=char(idx.participant(r)); sn=idx.sessNum(r);
    fp=char(idx.finalPath(r));
    if isKey(done,[id '|' task]), continue; end   % resume-skip
    fprintf('[%3d/%3d] %-11s %s\n', r, height(idx), task, id);
    if exist(fp,'file')~=2, fprintf('   missing on disk, skip\n'); continue; end
    try
        S=load(fp); fn=fieldnames(S); od=S.(fn{1}); clear S;
    catch e
        fprintf('   LOAD FAIL: %s\n', e.message); continue;
    end
    fs = 500; if isfield(od,'fs'), fs=double(od.fs); end

    % ---------- behavioral ----------
    try
        if isfield(od,'behDat') && istable(od.behDat)
            bd = od.behDat;
            switch task
                case 'cueTask',    Sc=H.cue_score(bd);
                case 'threshTask', Sc=H.thresh_score(bd);
                case 'O15',        Sc=H.o15_score(bd);
                otherwise,         Sc=[];
            end
            if ~isempty(Sc)
                fn2=fieldnames(Sc);
                fid=fopen(behFile,'a');
                nt = NaN; if isfield(Sc,'nTrials'), nt=Sc.nTrials; end
                for k=1:numel(fn2)
                    if strcmp(fn2{k},'nTrials'), continue; end
                    fprintf(fid,'%s,%s,%s,%s,%s,%d,%s,%.6g,%g\n', id,task,coh,grp,part,sn,fn2{k},Sc.(fn2{k}),nt);
                end
                fclose(fid);
            end
        end
    catch e
        fprintf('   BEH FAIL: %s\n', e.message);
    end

    % ---------- breathing block labels ----------
    if strcmp(task,'breathingTask') && isfield(od,'behDat') && istable(od.behDat)
        try
            bd=od.behDat; vn=bd.Properties.VariableNames;
            colstr = @(nm) getcol(bd,vn,nm);
            tcol=colstr('task'); ccol=colstr('condition'); ncol=colstr('noseMouth'); wcol=colstr('warp');
            key = strcat(tcol,'|',ccol,'|',ncol,'|',wcol);
            uk = unique(key,'stable');
            fidBlk=fopen(blkFile,'a');
            for u=1:numel(uk)
                parts=split(uk(u),'|');
                nb = sum(key==uk(u));
                fprintf(fidBlk,'%s,%s,%s,%s,%s,%s,%s,%d\n', id,coh,grp, ...
                    char(parts(1)),char(parts(2)),char(parts(3)),char(parts(4)),nb);
            end
            fclose(fidBlk);
        catch e
            fprintf('   BLK FAIL: %s\n', e.message);
        end
    end

    % ---------- macBP gamma ----------
    try
        labs = od.labels; labs = cellfun(@(x)char(string(x)),labs,'uni',0);
        isMac = cellfun(@(x) ~isempty(regexpi(x,'macbp','once')), labs);
        macIdx = find(isMac); nMac=numel(macIdx);
        if nMac>=1 && isfield(od,'data')
            peakFlatDb=nan(nMac,1); gammaDet=false(nMac,1); peakFreq=nan(nMac,1);
            peakHeight=nan(nMac,1); apExp=nan(nMac,1); spikeF=nan(nMac,1);
            for m=1:nMac
                Rm=H.macbp_gamma(double(od.data(macIdx(m),:)), fs, BAND);
                peakFlatDb(m)=Rm.peakFlatDb; gammaDet(m)=Rm.gammaDetected;
                peakFreq(m)=Rm.peakFreq; peakHeight(m)=Rm.peakHeight;
                apExp(m)=Rm.apExponent; spikeF(m)=Rm.spikeFrac;
            end
            [~,bi]=max(peakFlatDb);  % best = max peak flattened power
            fid=fopen(macFile,'a');
            for m=1:nMac
                fprintf(fid,'%s,%s,%s,%s,%s,%d,%d,%d,%s,%.6g,%d,%.6g,%.6g,%.6g,%.6g,%d\n', ...
                    id,task,coh,grp,part,sn,nMac,m,labs{macIdx(m)},peakFlatDb(m), ...
                    gammaDet(m), peakFreq(m), peakHeight(m), apExp(m), spikeF(m), double(m==bi));
            end
            fclose(fid);
            nDet = sum(gammaDet);
            nNonMaxDet = nDet - double(gammaDet(bi));
            fid=fopen(bestFile,'a');
            fprintf(fid,'%s,%s,%s,%s,%s,%d,%d,%d,%s,%.6g,%d,%.6g,%d,%d\n', ...
                id,task,coh,grp,part,sn,nMac,bi,labs{macIdx(bi)},peakFlatDb(bi), ...
                gammaDet(bi), spikeF(bi), nNonMaxDet, nDet);
            fclose(fid);
        else
            fprintf('   no macBP channels (nMac=%d)\n', nMac);
        end
    catch e
        fprintf('   MAC FAIL: %s\n', e.message);
    end
    clear od;
end
fprintf('\nDONE run_scores\n');
end

function s = getcol(bd, vn, nm)
% return a height(bd)x1 string column with <missing> replaced by "" (or all "" if absent)
    if any(strcmp(vn,nm))
        v = bd.(nm);
        try, s = string(v); catch, s = strings(height(bd),1); end
    else
        s = strings(height(bd),1);
    end
    s = s(:);
    if numel(s) ~= height(bd), s = strings(height(bd),1); end
    s(ismissing(s)) = "";
end
