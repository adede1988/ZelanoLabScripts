function behDat = build_behavior_table_O15(sniffs, rawBeh)
% Integrate O15 behavioral data with the detected sniff onsets.
% Shared front matter (first 6 columns) comes from behDatFromSniffs; the loop
% below is the TASK-SPECIFIC broadcast of O15 raw-behavior fields onto each sniff.

    behDat = behDatFromSniffs(sniffs, {"start", "free", "confirm"});

    for ii = 1:length(rawBeh.target)
        idx = find(behDat.n == ii);
        for jj = 1:length(idx)
            behDat.target(idx(jj))   = string(rawBeh.target{ii});
            behDat.response(idx(jj)) = string(rawBeh.response{ii});
            behDat.expScore(idx(jj)) = rawBeh.expScore(ii);
        end
    end
end
