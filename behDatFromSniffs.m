function behDat = behDatFromSniffs(sniffs, sniffTypes)
%BEHDATFROMSNIFFS  Shared first six columns of a per-sniff behavior table.
%
%   behDat = behDatFromSniffs(sniffs, sniffTypes)
%
%   Used by build_behavior_table_{O15,cueTask,threshTask}; each then broadcasts
%   its own task-specific raw-behavior fields onto these rows. (breathingTask
%   builds a per-breath table instead and does not use this.)
%
%   sniffs columns: 1=onset, 2=trialNum, 3=withinTrialIdx, 4=TTLoffset,
%   6=sniffType (integer index into sniffTypes).
%   sniffTypes: cell of readable labels, one per sniffType integer.

    behDat = table;
    behDat.sniffOnset = sniffs(:,1); % sniff onset
    behDat.n          = sniffs(:,2); % trial num
    behDat.wiTriali   = sniffs(:,3); % sniff within the trial
    behDat.TTLoffSet  = sniffs(:,4); % how far off from the TTL
    behDat.sniffType  = sniffs(:,6); % what kind of sniff is it?
    behDat.sniffLabel = [sniffTypes{sniffs(:,6)}]'; % readable sniff label
end
