function heartBeats = detectBeats(ECGz, beatSep, spec)
% detectBeats  Generic ECG beat detector parameterised by a spec string.
%
%   heartBeats = detectBeats(ECGz, beatSep, spec)
%
%   Replaces the 25 bespoke getBeats_<name> local functions that lived in
%   getSessionParams_breathingTask.m. Each legacy detector had the shape:
%
%       test = find(arrayfun(@(x,y,z) <conds>, ...
%                            ECGz(c1,a1:end-b1), ECGz(c2,a2:end-b2), ...));
%       test = test(diff(test) > beatSep);
%       heartBeats = test;
%
%   The channel slices are a relative-lag alignment. Each conjunct is encoded
%   as 'chan,lag,op,thr' where lag = a-1 (start sample minus 1) and
%   op in {gt,lt}. Conjuncts are joined by '&'. Example (legacy ...JH_1):
%       '1,2,gt,5 & 2,0,gt,4 & 3,0,lt,-0.5'
%
%   This reproduces the legacy slicing by truncating every channel to the
%   common length after applying its lag. (Every legacy detector has at least
%   one conjunct with b==0, so maxLag == k-1 and the truncated length L
%   matches the legacy slice length N-k+1 exactly.)
%
%   Inputs
%     ECGz    : [nChan x nSamples] z-scored (filtered, decimated) ECG.
%     beatSep : minimum sample separation between detected beats.
%     spec    : 'chan,lag,op,thr & chan,lag,op,thr & ...'  (op = gt|lt)
%
%   Output
%     heartBeats : 1xN vector of beat sample indices (within the lagged frame),
%                  same convention as the legacy getBeats_* functions.

    conj = parseSpec(spec);
    if isempty(conj)
        heartBeats = [];
        return;
    end

    N      = size(ECGz, 2);
    maxLag = max([conj.lag]);
    L      = N - maxLag;
    if L <= 1
        heartBeats = [];
        return;
    end

    mask = true(1, L);
    for i = 1:numel(conj)
        v = ECGz(conj(i).chan, 1+conj(i).lag : conj(i).lag+L);   % length L
        if strcmpi(conj(i).op, 'gt')
            mask = mask & (v >  conj(i).thr);
        else
            mask = mask & (v <  conj(i).thr);
        end
    end

    test = find(mask);
    test = test(diff(test) > beatSep);   % NB: drops the last element by construction
    heartBeats = test;
end


function conj = parseSpec(spec)
% Parse 'chan,lag,op,thr & chan,lag,op,thr & ...' into a struct array.
    conj = struct('chan', {}, 'lag', {}, 'op', {}, 'thr', {});
    if isempty(spec)
        return;
    end
    spec = char(spec);
    terms = strsplit(spec, '&');
    for t = 1:numel(terms)
        term = strtrim(terms{t});
        if isempty(term)
            continue;
        end
        parts = strsplit(term, ',');
        if numel(parts) ~= 4
            error('detectBeats:badSpec', ...
                  'Each conjunct must be "chan,lag,op,thr"; got "%s".', term);
        end
        c = struct();
        c.chan = str2double(strtrim(parts{1}));
        c.lag  = str2double(strtrim(parts{2}));
        c.op   = lower(strtrim(parts{3}));
        c.thr  = str2double(strtrim(parts{4}));
        if isnan(c.chan) || isnan(c.lag) || isnan(c.thr) || ...
                ~ismember(c.op, {'gt','lt'})
            error('detectBeats:badSpec', ...
                  'Could not parse conjunct "%s".', term);
        end
        conj(end+1) = c; %#ok<AGROW>
    end
end
