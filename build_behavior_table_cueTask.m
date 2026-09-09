function behDat = build_behavior_table_cueTask(sniffs, rawBeh)
% Integrate cueTask behavioral data with the detected sniff onsets.
% Shared front matter (first 6 columns) comes from behDatFromSniffs; the loop
% below is the TASK-SPECIFIC broadcast of cue raw-behavior fields onto each sniff.

    behDat = behDatFromSniffs(sniffs, {"cued"});

    for ii = 1:length(rawBeh.n)
        idx = find(behDat.n == ii);
        for jj = 1:length(idx)
            behDat.cue(idx(jj))      = rawBeh.cue(ii);
            behDat.odor(idx(jj))     = rawBeh.odor(ii);
            behDat.response(idx(jj)) = rawBeh.response(ii);
            if ~isempty(rawBeh.response_str{ii})
                behDat.respString(idx(jj)) = string(rawBeh.response_str(ii));
            else
                behDat.respString(idx(jj)) = "SKIP";
            end
            behDat.type(idx(jj)) = string(rawBeh.type(ii));
        end
    end
end
