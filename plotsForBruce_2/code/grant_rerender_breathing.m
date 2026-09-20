function grant_rerender_breathing(climUp)
% Re-render the 8 breathing group heatmaps (ALL/dupiS1/dupiS23/control x audiobook/focus)
% from their saved maps at a fixed color scale [4 climUp] (default 19.13, = the cueTask
% group scale) so all heatmaps are directly comparable. 60 Hz band zeroed; +200 ms axis.
if nargin<1||isempty(climUp), climUp=19.13; end
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); T=fullfile(proj,'out','tables');
groups={'all','control','dupiS1','dupiS23'}; tasks={'audiobook','focusedBreathing'};
for g=1:numel(groups), for t=1:numel(tasks)
  mp=fullfile(T,sprintf('grant_group_meanZ_%s_%s.mat',groups{g},tasks{t}));
  grant_group_render(mp, tasks{t}, [100 1000], [], [25 60], 57, climUp, 4, groups{g}, 'Time from +200 ms (ms)');
end, end
fprintf('re-rendered 8 breathing heatmaps at clim [4 %.2f]\n', climUp);
end
