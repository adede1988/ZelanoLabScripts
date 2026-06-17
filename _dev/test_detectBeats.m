function test_detectBeats()
% Validate detectBeats against all 25 legacy getBeats_* bodies (verbatim copies
% from getSessionParams_breathingTask.m). For several beatSep values and a fixed
% synthetic ECGz, detectBeats(ECGz,beatSep,spec) must equal getBeats_<name>(ECGz,beatSep)
% index-for-index.

    addpath('C:\Users\Adam\Documents\GitHub\ZelanoLabScripts');

    % name -> spec (chan,lag,op,thr & ...), derived from the legacy slices.
    M = {
        'getBeats_251205_EEG_NWU_AK',     '1,0,gt,3 & 3,10,gt,3'
        'getBeats_251113_EEG_NWU_GH',     '1,0,gt,2.5 & 3,10,gt,1'
        'getBeats_251009_EEG_NWU_SM',     '1,0,gt,3.5'
        'getBeats_251009_EEG_NWU_JM',     '3,0,lt,0 & 2,2,gt,3 & 1,17,gt,0.5'
        'getBeats_251008_EEG_NWU_GM',     '2,0,gt,2 & 1,3,gt,3 & 3,11,gt,1'
        'getBeats_251008_EEG_NWU_JC',     '3,0,lt,0 & 2,0,gt,2'
        'getBeats_251003_EEG_NWU_TI',     '1,0,lt,-1 & 2,0,gt,4'
        'getBeats_251027_Dupi_NMH_DL_1',  '1,0,gt,4'
        'getBeats_251120_Dupi_NMH_JL_1',  '1,0,gt,2.5'
        'getBeats_250818_Dupi_NMH_JH_1',  '1,2,gt,5 & 2,0,gt,4 & 3,0,lt,-0.5'
        'getBeats_250818_Dupi_NMH_JH_2',  '2,0,gt,3 & 3,1,lt,-3 & 3,13,gt,1'
        'getBeats_250623_DUPI_NMH_KS_2',  '1,0,gt,0.5 & 2,0,gt,1.75 & 3,0,lt,-2 & 2,9,lt,-1'
        'getBeats_250623_Dupi_NMH_KS_1',  '1,0,gt,2 & 2,7,gt,0.75 & 3,1,lt,-2'
        'getBeats_250908_OBE_NWU_AS',     '1,4,gt,2 & 2,0,gt,5 & 3,1,lt,-4'
        'getBeats_250723_EEG_NWU_IN',     '1,4,lt,-1 & 2,2,gt,2 & 3,0,lt,-1'
        'getBeats_250725_EEG_NWU_BN',     '3,11,gt,1 & 3,0,lt,-3 & 2,1,gt,1'
        'getBeats_250815_EEG_NWU_PP',     '1,0,lt,-4 & 2,0,gt,3'
        'getBeats_250819_EEG_NWU_ZL',     '1,0,lt,-2 & 2,0,gt,4 & 3,0,lt,-2'
        'getBeats_250723_EEG_NWU_BK',     '1,4,gt,1.5 & 2,0,gt,2 & 3,2,lt,-3'
        'getBeats_250912_EEG_NWU_JN',     '2,0,gt,4 & 3,0,lt,-4'
        'getBeats_250904_OBE_NWU_TI',     '1,0,lt,-2 & 2,1,gt,3 & 3,2,lt,0'
        'getBeats_250811_Dupi_NMH_TPB_1', '2,2,gt,1 & 3,0,lt,-1 & 2,12,lt,-1'
        'getBeats_250811_Dupi_NMH_TB_2',  '1,11,gt,1 & 2,1,gt,1 & 3,0,lt,-1'
        'getBeats_250929_Dupi_NMH_GH_1',  '3,0,lt,-4'
        'getBeats_251009_OBE_NWU_CP_1',   '3,0,lt,-3'
    };

    rng(7);
    ECGz = randn(3, 200000);
    % Make beats actually occur by injecting structure occasionally (not required
    % for correctness, but exercises the diff>beatSep path on non-empty sets).
    ECGz(:, 1:50:end) = ECGz(:, 1:50:end) * 4;

    beatSeps = [50, 200, 1000];
    nFail = 0; nTest = 0;
    for i = 1:size(M,1)
        fn   = str2func(M{i,1});
        spec = M{i,2};
        for bs = beatSeps
            nTest = nTest + 1;
            try
                legacy = fn(ECGz, bs);
            catch ME
                fprintf('LEGACY ERROR %s (beatSep=%d): %s\n', M{i,1}, bs, ME.message);
                nFail = nFail + 1; continue;
            end
            mine = detectBeats(ECGz, bs, spec);
            if ~isequal(legacy, mine)
                nFail = nFail + 1;
                fprintf('MISMATCH %-32s beatSep=%-5d  legacy=%d mine=%d  firstDiff=%d\n', ...
                    M{i,1}, bs, numel(legacy), numel(mine), firstDiff(legacy, mine));
            end
        end
    end
    fprintf('\ndetectBeats parity: %d/%d comparisons passed.\n', nTest-nFail, nTest);
    if nFail == 0
        fprintf('ALL_DETECTBEATS_PASS\n');
    else
        fprintf('DETECTBEATS_FAIL count=%d\n', nFail);
    end
end

function d = firstDiff(a, b)
    n = min(numel(a), numel(b));
    d = -1;
    for k = 1:n
        if a(k) ~= b(k), d = k; return; end
    end
    if numel(a) ~= numel(b), d = n+1; end
end

% ===================== legacy getBeats_* bodies (verbatim) =====================

function heartBeats = getBeats_251205_EEG_NWU_AK(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x > 3 & y>3, ...
                         ECGz(1,1:end-10),ECGz(3,11:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251113_EEG_NWU_GH(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x > 2.5 & y>1, ...
                         ECGz(1,1:end-10),ECGz(3,11:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251009_EEG_NWU_SM(ECGz, beatSep)
    test = find(arrayfun(@(x) x > 3.5, ...
                         ECGz(1,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251009_EEG_NWU_JM(ECGz, beatSep)
    test = find(arrayfun(@(x, y, z) x < 0 & y>3 & z>.5, ...
                         ECGz(3,1:end-17), ECGz(2,3:end-15), ...
                         ECGz(1,18:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251008_EEG_NWU_GM(ECGz, beatSep)
    test = find(arrayfun(@(x, y, z) x > 2 & y>3 & z>1, ...
                         ECGz(2,1:end-11), ECGz(1,4:end-8), ...
                         ECGz(3,12:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251008_EEG_NWU_JC(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x < 0 & y>2, ...
                         ECGz(3,1:end), ECGz(2,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251003_EEG_NWU_TI(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x < -1 & y>4, ...
                         ECGz(1,1:end), ECGz(2,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251027_Dupi_NMH_DL_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x > 4, ...
                         ECGz(1,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251120_Dupi_NMH_JL_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x > 2.5, ...
                         ECGz(1,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250818_Dupi_NMH_JH_1(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 5 & y > 4 & z < -0.5, ...
                         ECGz(1,3:end), ...
                         ECGz(2,1:end-2), ...
                         ECGz(3,1:end-2)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250818_Dupi_NMH_JH_2(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 3 & y < -3 & z > 1, ...
                         ECGz(2,1:end-13), ...
                         ECGz(3,2:end-12), ...
                         ECGz(3,14:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250623_DUPI_NMH_KS_2(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z,n) x > 0.5 & y > 1.75 & z < -2 & n < -1, ...
                         ECGz(1,1:end-9), ...
                         ECGz(2,1:end-9), ...
                         ECGz(3,1:end-9), ...
                         ECGz(2,10:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250623_Dupi_NMH_KS_1(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 2 & y > 0.75 & z < -2, ...
                         ECGz(1,1:end-7), ...
                         ECGz(2,8:end), ...
                         ECGz(3,2:end-6)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250908_OBE_NWU_AS(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 2 & y > 5 & z < -4, ...
                         ECGz(1,5:end), ...
                         ECGz(2,1:end-4), ...
                         ECGz(3,2:end-3)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250723_EEG_NWU_IN(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x < -1 & y > 2 & z < -1, ...
                         ECGz(1,5:end), ...
                         ECGz(2,3:end-2), ...
                         ECGz(3,1:end-4)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250725_EEG_NWU_BN(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1 & y < -3 & z > 1, ...
                         ECGz(3,12:end), ...
                         ECGz(3,1:end-11), ...
                         ECGz(2,2:end-10)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250815_EEG_NWU_PP(ECGz, beatSep)
    test = find(arrayfun(@(x,y) x < -4 & y > 3, ...
                         ECGz(1,1:end), ...
                         ECGz(2,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250819_EEG_NWU_ZL(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x < -2 & y > 4 & z < -2, ...
                         ECGz(1,1:end), ...
                         ECGz(2,1:end), ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250723_EEG_NWU_BK(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1.5 & y > 2 & z < -3, ...
                         ECGz(1,5:end), ...
                         ECGz(2,1:end-4), ...
                         ECGz(3,3:end-2)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250912_EEG_NWU_JN(ECGz, beatSep)
    test = find(arrayfun(@(x,y) x > 4 & y < -4, ...
                         ECGz(2,1:end), ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250904_OBE_NWU_TI(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x < -2 & y > 3 & z < 0, ...
                         ECGz(1,1:end-2), ...
                         ECGz(2,2:end-1), ...
                         ECGz(3,3:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250811_Dupi_NMH_TPB_1(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1 & y < -1 & z < -1, ...
                         ECGz(2,3:end-10), ...
                         ECGz(3,1:end-12), ...
                         ECGz(2,13:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250811_Dupi_NMH_TB_2(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1 & y > 1 & z < -1, ...
                         ECGz(1,12:end), ...
                         ECGz(2,2:end-10), ...
                         ECGz(3,1:end-11)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250929_Dupi_NMH_GH_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x < -4 , ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251009_OBE_NWU_CP_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x < -3 , ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end
