% Capture build_behavior_table_{O15,cueTask,threshTask} outputs on deterministic
% synthetic inputs, to use as an oracle when refactoring the shared front matter.
addpath(fileparts(fileparts(mfilename('fullpath'))));

nT = 5; spt = 2; N = nT*spt;
sniffs = zeros(N,6);
sniffs(:,1) = (1:N)'*100;             % onset
sniffs(:,2) = repelem((1:nT)',spt);   % trial num
sniffs(:,3) = repmat((1:spt)',nT,1);  % within-trial idx
sniffs(:,4) = (1:N)'-3;               % TTL offset

sn_O15 = sniffs; sn_O15(:,6) = mod((0:N-1)',3)+1;  % O15 has 3 sniff types
sn_one = sniffs; sn_one(:,6) = ones(N,1);          % cue/thresh have 1

bO = table();
bO.target   = arrayfun(@(x){sprintf('t%d',x)}, (1:nT)');
bO.response = arrayfun(@(x){sprintf('r%d',x)}, (1:nT)');
bO.expScore = (1:nT)' + 0.5;

bC = table();
bC.n            = (1:nT)';
bC.cue          = (1:nT)';
bC.odor         = (nT:-1:1)';
bC.response     = (1:nT)'*2;
bC.response_str = {'Yes'; ''; 'No'; 'Yes'; ''};
bC.type         = ["hit"; "miss"; "cr"; "fa"; "skip"];

bT = table();
bT.trialNum     = (1:nT)';
bT.Odor         = mod((0:nT-1)',3)+1;
bT.pleasantness = (1:nT)'/2;
bT.intensity    = (nT:-1:1)'/2;

o15 = build_behavior_table_O15(sn_O15, bO);
cue = build_behavior_table_cueTask(sn_one, bC);
thr = build_behavior_table_threshTask(sn_one, bT);

save(fullfile(fileparts(mfilename('fullpath')),'builder_oracle.mat'), ...
     'sn_O15','sn_one','bO','bC','bT','o15','cue','thr');
fprintf('oracle saved: o15=%dx%d cue=%dx%d thr=%dx%d\n', ...
        height(o15),width(o15),height(cue),width(cue),height(thr),width(thr));
