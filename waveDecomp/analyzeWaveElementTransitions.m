function trans = analyzeWaveElementTransitions(result, outDir)
% analyzeWaveElementTransitions
%
% Uses result.events from decomposeOceanWaveElements() to compute:
%   1. event-to-event transition counts
%   2. transition probability matrix
%   3. temporal overlap counts
%   4. overlap probability matrix
%   5. lag distributions
%   6. diagnostic plots
%
% Example:
%   result = decomposeOceanWaveElements("beach.wav", "wave_decomposition");
%   trans = analyzeWaveElementTransitions(result, "wave_decomposition");

if nargin < 2
    outDir = result.outDir;
end

if ~exist(outDir, "dir")
    mkdir(outDir);
end

events = result.events;

if isempty(events)
    error("result.events is empty. No events available for transition analysis.");
end

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------
nElements = max(events.elementID);

maxTransitionLagSec = 8.0;
% Only count B as following A if B starts within this many seconds
% after A starts.

allowSelfTransitions = true;
% If false, A -> A transitions are ignored.

usePeakTimes = false;
% If true, use peakSec order instead of onsetSec order.

% ------------------------------------------------------------
% Sort events in temporal order
% ------------------------------------------------------------
if usePeakTimes
    timeField = "peakSec";
else
    timeField = "onsetSec";
end

events = sortrows(events, timeField);

eventTimes = events.(timeField);
eventOnsets = events.onsetSec;
eventOffsets = events.offsetSec;
eventIDs = events.elementID;

% ------------------------------------------------------------
% Sequential transition matrix
% ------------------------------------------------------------
transitionCounts = zeros(nElements, nElements);
transitionLags = cell(nElements, nElements);

for e = 1:(height(events) - 1)

    fromID = eventIDs(e);

    % Candidate future events
    futureIdx = (e + 1):height(events);

    futureOnsets = eventOnsets(futureIdx);
    futureIDs = eventIDs(futureIdx);

    lags = futureOnsets - eventOnsets(e);

    % Keep only events after current event onset and within lag window
    valid = lags >= 0 & lags <= maxTransitionLagSec;

    if ~allowSelfTransitions
        valid = valid & futureIDs ~= fromID;
    end

    if ~any(valid)
        continue
    end

    % Use the first valid event as the next transition
    validFutureIdx = futureIdx(valid);
    nextIdx = validFutureIdx(1);

    toID = eventIDs(nextIdx);
    lag = eventOnsets(nextIdx) - eventOnsets(e);

    transitionCounts(fromID, toID) = transitionCounts(fromID, toID) + 1;
    transitionLags{fromID, toID} = [transitionLags{fromID, toID}; lag];
end

% Convert counts to row-normalized probabilities
transitionProbability = transitionCounts ./ sum(transitionCounts, 2);
transitionProbability(isnan(transitionProbability)) = 0;

% ------------------------------------------------------------
% Overlap matrix
% ------------------------------------------------------------
overlapCounts = zeros(nElements, nElements);
overlapDurations = cell(nElements, nElements);

for a = 1:height(events)

    idA = eventIDs(a);
    onsetA = eventOnsets(a);
    offsetA = eventOffsets(a);

    for b = 1:height(events)

        if a == b
            continue
        end

        idB = eventIDs(b);
        onsetB = eventOnsets(b);
        offsetB = eventOffsets(b);

        overlapDur = min(offsetA, offsetB) - max(onsetA, onsetB);

        if overlapDur > 0
            overlapCounts(idA, idB) = overlapCounts(idA, idB) + 1;
            overlapDurations{idA, idB} = [overlapDurations{idA, idB}; overlapDur];
        end
    end
end

% Probability that B overlaps A, given an A event occurred
nEventsPerElement = accumarray(eventIDs, 1, [nElements 1], @sum, 0);
overlapProbability = overlapCounts ./ nEventsPerElement;
overlapProbability(isnan(overlapProbability)) = 0;

% ------------------------------------------------------------
% Mean lags and overlap durations
% ------------------------------------------------------------
meanTransitionLag = nan(nElements, nElements);
medianTransitionLag = nan(nElements, nElements);
meanOverlapDuration = nan(nElements, nElements);

for i = 1:nElements
    for j = 1:nElements

        if ~isempty(transitionLags{i,j})
            meanTransitionLag(i,j) = mean(transitionLags{i,j});
            medianTransitionLag(i,j) = median(transitionLags{i,j});
        end

        if ~isempty(overlapDurations{i,j})
            meanOverlapDuration(i,j) = mean(overlapDurations{i,j});
        end
    end
end

% ------------------------------------------------------------
% Get labels
% ------------------------------------------------------------
labels = "E" + string(1:nElements);

if isfield(result, "componentInfo") && ...
        istable(result.componentInfo) && ...
        any(strcmp(result.componentInfo.Properties.VariableNames, "suggestedLabel"))

    labels = strings(nElements, 1);

    for k = 1:nElements
        row = result.componentInfo.elementID == k;
        if any(row)
            labels(k) = sprintf("E%d: %s", ...
                k, result.componentInfo.suggestedLabel(find(row, 1)));
        else
            labels(k) = "E" + string(k);
        end
    end
end

% ------------------------------------------------------------
% Save matrices
% ------------------------------------------------------------
transitionCountTable = array2table(transitionCounts, ...
    "VariableNames", matlab.lang.makeValidName("to_" + labels), ...
    "RowNames", matlab.lang.makeValidName("from_" + labels));

transitionProbTable = array2table(transitionProbability, ...
    "VariableNames", matlab.lang.makeValidName("to_" + labels), ...
    "RowNames", matlab.lang.makeValidName("from_" + labels));

overlapProbTable = array2table(overlapProbability, ...
    "VariableNames", matlab.lang.makeValidName("overlaps_" + labels), ...
    "RowNames", matlab.lang.makeValidName("given_" + labels));

writetable(transitionCountTable, ...
    fullfile(outDir, "transition_counts.csv"), ...
    "WriteRowNames", true);

writetable(transitionProbTable, ...
    fullfile(outDir, "transition_probabilities.csv"), ...
    "WriteRowNames", true);

writetable(overlapProbTable, ...
    fullfile(outDir, "overlap_probabilities.csv"), ...
    "WriteRowNames", true);

% ------------------------------------------------------------
% Plot transition probability heatmap
% ------------------------------------------------------------
fig1 = figure("Color", "w", "Position", [100 100 900 750]);

imagesc(transitionProbability);
axis square;
colorbar;

xticks(1:nElements);
yticks(1:nElements);
xticklabels("E" + string(1:nElements));
yticklabels("E" + string(1:nElements));

xlabel("Following element");
ylabel("Current element");
title("Transition probability: P(next element | current element)");

for i = 1:nElements
    for j = 1:nElements
        if transitionProbability(i,j) > 0
            text(j, i, sprintf("%.2f", transitionProbability(i,j)), ...
                "HorizontalAlignment", "center", ...
                "FontSize", 9);
        end
    end
end

saveas(fig1, fullfile(outDir, "transition_probability_heatmap.png"));

% ------------------------------------------------------------
% Plot overlap probability heatmap
% ------------------------------------------------------------
fig2 = figure("Color", "w", "Position", [100 100 900 750]);

imagesc(overlapProbability);
axis square;
colorbar;

xticks(1:nElements);
yticks(1:nElements);
xticklabels("E" + string(1:nElements));
yticklabels("E" + string(1:nElements));

xlabel("Overlapping element");
ylabel("Reference element");
title("Overlap probability: P(element B overlaps element A | element A occurred)");

for i = 1:nElements
    for j = 1:nElements
        if overlapProbability(i,j) > 0
            text(j, i, sprintf("%.2f", overlapProbability(i,j)), ...
                "HorizontalAlignment", "center", ...
                "FontSize", 9);
        end
    end
end

saveas(fig2, fullfile(outDir, "overlap_probability_heatmap.png"));

% ------------------------------------------------------------
% Plot directed transition graph
% ------------------------------------------------------------
thresholdProb = 0.10;
% Only show transitions with probability >= thresholdProb

A = transitionProbability;
A(A < thresholdProb) = 0;

G = digraph(A, "omitselfloops");

fig3 = figure("Color", "w", "Position", [100 100 1100 800]);

p = plot(G, ...
    "Layout", "layered", ...
    "NodeLabel", "E" + string(1:nElements), ...
    "LineWidth", 1 + 8 * G.Edges.Weight, ...
    "ArrowSize", 12);

title(sprintf("Wave-element transition graph, edges >= %.2f probability", thresholdProb));

saveas(fig3, fullfile(outDir, "transition_graph.png"));

% ------------------------------------------------------------
% Make an event-level transition table
% ------------------------------------------------------------
transitionRows = table;
row = 0;

for i = 1:nElements
    for j = 1:nElements

        if transitionCounts(i,j) == 0
            continue
        end

        row = row + 1;

        transitionRows.fromElement(row, 1) = i;
        transitionRows.toElement(row, 1) = j;
        transitionRows.fromLabel(row, 1) = labels(i);
        transitionRows.toLabel(row, 1) = labels(j);
        transitionRows.count(row, 1) = transitionCounts(i,j);
        transitionRows.probability(row, 1) = transitionProbability(i,j);
        transitionRows.meanLagSec(row, 1) = meanTransitionLag(i,j);
        transitionRows.medianLagSec(row, 1) = medianTransitionLag(i,j);
    end
end

transitionRows = sortrows(transitionRows, "probability", "descend");

writetable(transitionRows, fullfile(outDir, "transition_summary.csv"));

% ------------------------------------------------------------
% Return results
% ------------------------------------------------------------
trans = struct;
trans.events = events;
trans.labels = labels;
trans.transitionCounts = transitionCounts;
trans.transitionProbability = transitionProbability;
trans.transitionLags = transitionLags;
trans.meanTransitionLag = meanTransitionLag;
trans.medianTransitionLag = medianTransitionLag;
trans.overlapCounts = overlapCounts;
trans.overlapProbability = overlapProbability;
trans.overlapDurations = overlapDurations;
trans.meanOverlapDuration = meanOverlapDuration;
trans.transitionSummary = transitionRows;

fprintf("Transition analysis complete.\n");
fprintf("Saved outputs to: %s\n", outDir);

end