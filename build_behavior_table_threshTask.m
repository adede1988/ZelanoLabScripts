function behDat = build_behavior_table_threshTask(sniffs, rawBeh)
% Integrate threshTask (PEA) behavioral data with the detected sniff onsets.
% Shared front matter (first 6 columns) comes from behDatFromSniffs; the loop
% below is the TASK-SPECIFIC broadcast of thresh raw-behavior fields onto each sniff.

    trialTypes = {'air', 'low', 'med'};
    behDat = behDatFromSniffs(sniffs, {"cued"});

    for ii = 1:length(rawBeh.trialNum)
        idx = find(behDat.n == ii);
        for jj = 1:length(idx)
            behDat.odor(idx(jj))         = rawBeh.Odor(ii);
            behDat.pleasantness(idx(jj)) = rawBeh.pleasantness(ii);
            behDat.intensity(idx(jj))    = rawBeh.intensity(ii);
            behDat.type(idx(jj))         = string(trialTypes{rawBeh.Odor(ii)});
        end
    end
end
