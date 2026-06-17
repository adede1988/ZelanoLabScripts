function test_beatSpec_map()
% For every CURATED breathing session, confirm the sheet beatSpec (via applyParams)
% reproduces the exact legacy getBeats_<name> that getSessionParams_breathingTask's
% switch assigns to that session. Pure compute, no .mat loads.

    repo = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';
    addpath(repo);
    xlsx = fullfile(repo, 'dataTracking.xlsx');

    % sheet sessionID (lower) -> legacy getBeats_<name> assigned by the switch
    M = {
        '250818_dupi_nmh_jh_1',   'getBeats_250818_Dupi_NMH_JH_1'
        '250623_dupi_nmh_ks_2',   'getBeats_250623_DUPI_NMH_KS_2'
        '250623_dupi_nmh_ks_1',   'getBeats_250623_Dupi_NMH_KS_1'
        '250908_obe_nwu_as',      'getBeats_250908_OBE_NWU_AS'
        '250723_eeg_nwu_in',      'getBeats_250723_EEG_NWU_IN'
        '250725_eeg_nwu_bn',      'getBeats_250725_EEG_NWU_BN'
        '250815_eeg_nwu_pp',      'getBeats_250815_EEG_NWU_PP'
        '250819_eeg_nwu_zl',      'getBeats_250819_EEG_NWU_ZL'
        '250723_eeg_nwu_bk',      'getBeats_250723_EEG_NWU_BK'
        '250912_eeg_nwu_jn',      'getBeats_250912_EEG_NWU_JN'
        '250904_obe_nwu_ti_1',    'getBeats_250904_OBE_NWU_TI'
        '250818_dupi_nmh_jh_2',   'getBeats_250818_Dupi_NMH_JH_2'
        '250811_dupi_nmh_tpb_1',  'getBeats_250811_Dupi_NMH_TPB_1'
        '250811_dupi_nmh_tb_2',   'getBeats_250811_Dupi_NMH_TB_2'
        '250929_dupi_nmh_gh_1',   'getBeats_250929_Dupi_NMH_GH_1'
        '251009_obe_nwu_cp_1',    'getBeats_251009_OBE_NWU_CP_1'
        '251002_dupi_nmh_ab_1',   'getBeats_251009_OBE_NWU_CP_1'
        '251027_dupi_nmh_dl_1',   'getBeats_251027_Dupi_NMH_DL_1'
        '250929_dupi_nmh_gh_2',   'getBeats_251027_Dupi_NMH_DL_1'
        '251002_dupi_nmh_ab_2',   'getBeats_251027_Dupi_NMH_DL_1'
        '251013_dupi_nmh_jn_2',   'getBeats_251027_Dupi_NMH_DL_1'
        '251030_dupi_nmh_db_1',   'getBeats_251027_Dupi_NMH_DL_1'
        '251030_dupi_nmh_db_2',   'getBeats_251027_Dupi_NMH_DL_1'
        '251110_dupi_nmh_pc_1',   'getBeats_251027_Dupi_NMH_DL_1'
        '251120_dupi_nmh_jl_1',   'getBeats_251120_Dupi_NMH_JL_1'
        '250818_dupi_nmh_jh_3',   'getBeats_251027_Dupi_NMH_DL_1'
        '251003_eeg_nwu_ti',      'getBeats_251003_EEG_NWU_TI'
        '251008_eeg_nwu_jc',      'getBeats_251008_EEG_NWU_JC'
        '251008_eeg_nwu_gm',      'getBeats_251008_EEG_NWU_GM'
        '251009_eeg_nwu_jm',      'getBeats_251009_EEG_NWU_JM'
        '251009_eeg_nwu_sm',      'getBeats_251009_EEG_NWU_SM'
        '251027_eeg_nwu_as',      'getBeats_251009_EEG_NWU_SM'
        '251105_eeg_nwu_gl',      'getBeats_251009_EEG_NWU_SM'
        '251110_eeg_nwu_ga',      'getBeats_251009_EEG_NWU_SM'
        '251111_eeg_nwu_vw',      'getBeats_251009_EEG_NWU_SM'
        '251113_eeg_nwu_gh',      'getBeats_251113_EEG_NWU_GH'
        '251118_eeg_nwu_adtest',  'getBeats_251009_EEG_NWU_SM'
        '251202_eeg_nwu_gj',      'getBeats_251009_EEG_NWU_SM'
        '260109_eeg_nwu_aa',      'getBeats_251113_EEG_NWU_GH'
        '251205_eeg_nwu_ak',      'getBeats_251205_EEG_NWU_AK'
        '251208_eeg_nwu_za',      'getBeats_251009_EEG_NWU_SM'
    };

    cfg = applyParams('breathingTask', 'makeOutDat', xlsx);
    sheetIDs = lower(cfg.sessionIDs);

    rng(3); ECGz = randn(3, 150000); ECGz(:,1:37:end) = ECGz(:,1:37:end)*4;
    beatSeps = [60, 300];

    nMis = 0; nChk = 0; nNoOracle = 0;
    for i = 1:numel(sheetIDs)
        id = sheetIDs{i};
        mi = find(strcmp(id, M(:,1)), 1);
        if isempty(mi)
            nNoOracle = nNoOracle + 1;   % guess/future session, no legacy oracle
            continue;
        end
        spec = applyParams('breathingTask', cfg.sessionIDs{i}, xlsx).beatSpec;
        fn   = str2func(M{mi,2});
        for bs = beatSeps
            nChk = nChk + 1;
            legacy = fn(ECGz, bs);
            mine   = detectBeats(ECGz, bs, spec);
            if ~isequal(legacy, mine)
                nMis = nMis + 1;
                fprintf('  MISMATCH %-24s spec="%s" vs %s\n', cfg.sessionIDs{i}, spec, M{mi,2});
            end
        end
    end
    fprintf('beatSpec map: %d checks, %d ok, %d mismatch; %d guess sessions w/o oracle.\n', ...
        nChk, nChk-nMis, nMis, nNoOracle);
    if nMis == 0, fprintf('ALL_BEATSPEC_MAP_PASS\n'); else, fprintf('BEATSPEC_MAP_FAIL\n'); end
end

% ===================== legacy getBeats_* bodies (verbatim) =====================
function heartBeats = getBeats_251205_EEG_NWU_AK(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x > 3 & y>3, ECGz(1,1:end-10),ECGz(3,11:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251113_EEG_NWU_GH(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x > 2.5 & y>1, ECGz(1,1:end-10),ECGz(3,11:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251009_EEG_NWU_SM(ECGz, beatSep)
    test = find(arrayfun(@(x) x > 3.5, ECGz(1,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251009_EEG_NWU_JM(ECGz, beatSep)
    test = find(arrayfun(@(x, y, z) x < 0 & y>3 & z>.5, ECGz(3,1:end-17), ECGz(2,3:end-15), ECGz(1,18:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251008_EEG_NWU_GM(ECGz, beatSep)
    test = find(arrayfun(@(x, y, z) x > 2 & y>3 & z>1, ECGz(2,1:end-11), ECGz(1,4:end-8), ECGz(3,12:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251008_EEG_NWU_JC(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x < 0 & y>2, ECGz(3,1:end), ECGz(2,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251003_EEG_NWU_TI(ECGz, beatSep)
    test = find(arrayfun(@(x, y) x < -1 & y>4, ECGz(1,1:end), ECGz(2,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251027_Dupi_NMH_DL_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x > 4, ECGz(1,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251120_Dupi_NMH_JL_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x > 2.5, ECGz(1,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250818_Dupi_NMH_JH_1(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 5 & y > 4 & z < -0.5, ECGz(1,3:end), ECGz(2,1:end-2), ECGz(3,1:end-2)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250818_Dupi_NMH_JH_2(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 3 & y < -3 & z > 1, ECGz(2,1:end-13), ECGz(3,2:end-12), ECGz(3,14:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250623_DUPI_NMH_KS_2(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z,n) x > 0.5 & y > 1.75 & z < -2 & n < -1, ECGz(1,1:end-9), ECGz(2,1:end-9), ECGz(3,1:end-9), ECGz(2,10:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250623_Dupi_NMH_KS_1(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 2 & y > 0.75 & z < -2, ECGz(1,1:end-7), ECGz(2,8:end), ECGz(3,2:end-6)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250908_OBE_NWU_AS(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 2 & y > 5 & z < -4, ECGz(1,5:end), ECGz(2,1:end-4), ECGz(3,2:end-3)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250723_EEG_NWU_IN(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x < -1 & y > 2 & z < -1, ECGz(1,5:end), ECGz(2,3:end-2), ECGz(3,1:end-4)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250725_EEG_NWU_BN(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1 & y < -3 & z > 1, ECGz(3,12:end), ECGz(3,1:end-11), ECGz(2,2:end-10)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250815_EEG_NWU_PP(ECGz, beatSep)
    test = find(arrayfun(@(x,y) x < -4 & y > 3, ECGz(1,1:end), ECGz(2,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250819_EEG_NWU_ZL(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x < -2 & y > 4 & z < -2, ECGz(1,1:end), ECGz(2,1:end), ECGz(3,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250723_EEG_NWU_BK(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1.5 & y > 2 & z < -3, ECGz(1,5:end), ECGz(2,1:end-4), ECGz(3,3:end-2)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250912_EEG_NWU_JN(ECGz, beatSep)
    test = find(arrayfun(@(x,y) x > 4 & y < -4, ECGz(2,1:end), ECGz(3,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250904_OBE_NWU_TI(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x < -2 & y > 3 & z < 0, ECGz(1,1:end-2), ECGz(2,2:end-1), ECGz(3,3:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250811_Dupi_NMH_TPB_1(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1 & y < -1 & z < -1, ECGz(2,3:end-10), ECGz(3,1:end-12), ECGz(2,13:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250811_Dupi_NMH_TB_2(ECGz, beatSep)
    test = find(arrayfun(@(x,y,z) x > 1 & y > 1 & z < -1, ECGz(1,12:end), ECGz(2,2:end-10), ECGz(3,1:end-11)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_250929_Dupi_NMH_GH_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x < -4 , ECGz(3,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
function heartBeats = getBeats_251009_OBE_NWU_CP_1(ECGz, beatSep)
    test = find(arrayfun(@(x) x < -3 , ECGz(3,1:end)));
    test = test(diff(test) > beatSep); heartBeats = test;
end
