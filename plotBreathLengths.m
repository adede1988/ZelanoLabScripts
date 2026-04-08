function [] = plotBreathLengths(outDat, R)

idx = outDat.behDat.sniffOnset(outDat.behDat.goodBreath == 1); 

test = arrayfun(@(x) R.rspDat(x-1000:x+7000), idx, 'uniformoutput', false);
test = cat(1, test{:});

lenVals = outDat.behDat.length(outDat.behDat.goodBreath == 1); 
figure('visible', false, 'position', [0,0,1000,500]);
[~, order] = sort(lenVals);
imagesc(-2:.002:14, [], test(order,:))
caxis([-100,100])
xlabel('time from breath onset (s)')
ylabel('breaths sorted by length')
title(sprintf('Variability in breath lengths (%s)', ...
                  outDat.sessID), ...
          'Interpreter','none');
saveas(gcf,fullfile(outDat.figs, ['breathLengths' '.jpg']));

idx = find(cellfun(@(x) strcmp(x, 'RRint'), outDat.labels));
rrDat = outDat.data(idx, :); 
idx = outDat.behDat.sniffOnset(outDat.behDat.goodBreath == 1); 
test = arrayfun(@(x) rrDat(x-1000:x+5000), idx, 'uniformoutput', false);
test = cat(1, test{:});

% test = (test - mean(test, 2)) ./ std(test, [], 2);


figure('visible', false, 'position', [0,0,1000,500]);
[~, order] = sort(lenVals);
imagesc(-2:.002:14, [], test(order,:))
xlabel('time from breath onset (s)')
ylabel('breaths sorted by length')
colorbar
title(sprintf('Heart Rate Across Breaths (%s)', ...
                  outDat.sessID), ...
          'Interpreter','none');
saveas(gcf,fullfile(outDat.figs, ['HeartByBreathLengths' '.jpg']));


end