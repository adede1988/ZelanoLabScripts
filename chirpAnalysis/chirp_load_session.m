function S = chirp_load_session(fp, C)
% CHIRP_LOAD_SESSION  Load a cue final and extract the chirp-analysis substrate.
%   S = chirp_load_session(fp, C)
%   Robust load (outDat -> chanDat -> out). Selects bestMac (OB channel) + the full macBP set,
%   the sniff-locked onsets (behDat.finalOnset), and respiration. Channels by LABEL, never index.
%
%   S fields: sessID fs labs macIdx macLabs macSig[nMac x nSamp] goodMac bestMac bestIdx
%             bestSig onsets behDat rsp ok msg

    S = struct('ok',false,'msg','');
    try
        s = load(fp); fn = fieldnames(s); od = s.(fn{1}); clear s;
    catch ME, S.msg = ['load failed: ' ME.message]; return; end
    if ~isfield(od,'data') || ~isfield(od,'labels') || ~isfield(od,'fs')
        S.msg = 'missing data/labels/fs'; return; end

    labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
    S.labs = labs; S.fs = od.fs;
    [~,nm] = fileparts(fp); S.sessID = erase(nm, '_cueTaskPreproc');

    isMac = cellfun(@(x) contains(x,'macBP'), labs);
    S.macIdx = find(isMac); S.macLabs = labs(S.macIdx);
    S.macSig = double(od.data(S.macIdx, :));
    % channel-level QC: finite + non-flat
    S.goodMac = arrayfun(@(r) all(isfinite(S.macSig(r,:))) && std(S.macSig(r,:))>0, (1:numel(S.macIdx))');

    % bestMac (OB channel)
    bestIdx = [];
    if isfield(od,'bestMac') && ~isempty(od.bestMac)
        bestIdx = find(strcmp(char(string(od.bestMac)), labs), 1);
        S.bestMac = char(string(od.bestMac));
    end
    if isempty(bestIdx)
        if ~isempty(S.macIdx), bestIdx = S.macIdx(1); S.bestMac = labs{bestIdx};
            S.msg = 'bestMac ABSENT -> used first macBP'; else, S.msg='no macBP'; return; end
    end
    S.bestIdx = bestIdx; S.bestSig = double(od.data(bestIdx,:));

    % onsets: cue = one cued sniff/trial; use behDat.finalOnset
    if isfield(od,'behDat') && istable(od.behDat) && ismember('finalOnset', od.behDat.Properties.VariableNames)
        bd = od.behDat;
        keep = true(height(bd),1);
        if ismember('sniffLabel', bd.Properties.VariableNames)
            keep = keep & (string(bd.sniffLabel) == "cued" | all(string(bd.sniffLabel)=="cued"));
        end
        S.behDat = bd(keep,:);
        S.onsets = double(bd.finalOnset(keep));
    else
        S.msg = 'no behDat.finalOnset'; return;
    end

    % respiration (for figures)
    isR = cellfun(@(x) contains(x,'rsp'), labs);
    if any(isR)
        rA = double(od.data(isR,:)); ri = 1;
        if isfield(od,'rspIDX') && od.rspIDX>=1 && od.rspIDX<=size(rA,1), ri = od.rspIDX; end
        fl = 1; if isfield(od,'rspFlip'), fl = od.rspFlip; end
        S.rsp = rA(ri,:).*fl;
    else, S.rsp = []; end

    S.ok = true;
end
