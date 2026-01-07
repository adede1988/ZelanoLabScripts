function [raw, P] = getSessionParams_breathingTask(S)
% getSessionParams_cueTask
%   Loads the prebuilt outDat for the breathing task and aggregates
%   all subject-specific parameters into P (centralizing what
%   was hard-coded in breathingPreProc_scratch.m).
%
% Inputs
%   S.id   : session ID (e.g., '250904_OBE_NWU_TI')
%   S.root : folder containing <S.id>.mat (with outDat inside)
%
% Outputs
%   raw : struct similar to O15 loader outputs
%         .sessID, .fs_raw, .data, .labels, .beh, .TTL (if present)
%   P   : parameters for preprocessing and event handling
%         .task, .type, .fs_target, .debug
%         .computeResp, .rspIDX, .rspFlip
%         .hasEEG, .spikeThresh, .spikeWin, .spikeClean
%         .respThresh, .cuedBackBuff, .adjWin
%         .ttlMap (names to look for if TTL is a table/struct)

    arguments
        S struct
    end

    %% --- Load prebuilt outDat ---
    matPath = fullfile(S.root, S.id, 'preProc', [S.id '_breathingPreProc.mat']);
    if ~exist(matPath, 'file')
        error('getSessionParams_breathingTask:MissingMat', ...
              'Expected MAT at %s (containing outDat).', matPath);
    end
    tmp = load(matPath);
    if ~isfield(tmp, 'outDat')  && ~isfield(tmp, 'out') 
        error('getSessionParams_breathingTask:NoOutDat', ...
              'MAT must contain variable outDat.');
    end
    try
        od = tmp.outDat;
    catch
        od = tmp.out; 
    end


    % raw view to match O15 pipeline expectations
    raw = struct();
    raw.sessID = char(od.sessID);
    raw.fs_raw = od.fs;
    raw.data   = od.data;
    raw.labels = od.labels;
    raw.beh    = od.behDat;

    % create TTLs hard coded based on assumption of 5 min windows: 
    TTL = [0:600000:size(od.data,2)]; 
    TTL(1) = 1; 
    TTL(end) = []; 
    od.TTL = TTL; 

    raw.paths  = struct('root', S.root);

    sessID = char(S.id);
    if contains(sessID, 'OBE', 'IgnoreCase', true) || contains(char(S.root), 'OBE', 'IgnoreCase', true)
        raw.type = 'OBE';
    else
        raw.type = 'Dupi';
    end
    %% --- Defaults (session-agnostic) ---
    P = struct();
    P.task        = 'breathingTask';
    P.type        = raw.type;
    P.fs_target   = 500;
    P.debug       = false;

    % Respiration (optional — leave on so downstream helpers work)
    P.computeResp = true;
    P.rspIDX      = 1;    
    P.rspFlip     = 1;    

    % EEG / spike cleaning defaults
    P.hasMacros = true; 
    P.hasEEG      = true;
    P.spikeThresh = 20;
    P.spikeWin    = 11;
    P.spikeClean  = true;
    P.macroRemove = []; %are there macro channels to be removed? 

    % % Respiration onset metric & alignment defaults
    % P.respThresh   = 500;   % scratch "otherwise" default
    % P.cuedBackBuff = 150;   % scratch "otherwise" default
    % P.adjWin       = 500;   % scratch default, except TI case

  

    %% --- Subject-specific overrides from breathingTaskPreProc_scratch.m ---

    switch raw.sessID

        % ---------------- participants from sessionIDs / rspIDX / rspFlip ----------------
        case '250818_Dupi_NMH_JH_1'
            P.rspIDX      = 3;
            P.rspFlip     = -1;
            P.spikeThresh = 20;
            % P.spikeWin    = 11;
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                getBeats_250818_Dupi_NMH_JH_1(ECGz, beatSep);
    
        case '250623_DUPI_NMH_KS_2'
            P.rspIDX      = 3;
            P.rspFlip     = -1;
            P.spikeThresh = 15;
            % P.spikeWin    = 11;
            P.hasEEG      = true;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250623_DUPI_NMH_KS_2(ECGz, beatSep);
    
        case '250623_Dupi_NMH_KS_1'
            P.rspIDX      = 3;
            P.rspFlip     = -1;
            P.spikeThresh = 15;
            % P.spikeWin    = 11;
            P.hasEEG      = false;
            % P.spikeClean  = false;
            P.macroRemove = 6; 
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250623_Dupi_NMH_KS_1(ECGz, beatSep);
    
        case '250908_OBE_NWU_AS'
            P.rspIDX      = 3;
            P.rspFlip     = -1;
            P.spikeThresh = 20;
            P.spikeWin    = 11;
            P.hasEEG      = true;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250908_OBE_NWU_AS(ECGz, beatSep);
    
        case '250723_EEG_NWU_IN'
            P.rspIDX      = 1;
            P.rspFlip     = -1;
            P.hasEEG      = true;
            P.hasMacros   = false;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250723_EEG_NWU_IN(ECGz, beatSep);
    
        case '250725_EEG_NWU_BN'
            P.rspIDX      = 1;
            P.rspFlip     = -1;
            P.hasEEG      = true;
            P.hasMacros   = false;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250725_EEG_NWU_BN(ECGz, beatSep);
    
        case '250815_EEG_NWU_PP'
            P.rspIDX      = 1;
            P.rspFlip     = -1;
            P.hasEEG      = true;
            P.hasMacros   = false;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250815_EEG_NWU_PP(ECGz, beatSep);
    
        case '250819_EEG_NWU_ZL'
            P.rspIDX      = 1;
            P.rspFlip     =  1;
            P.hasEEG      = true;
            P.hasMacros   = false;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250819_EEG_NWU_ZL(ECGz, beatSep);
    
        case '250723_EEG_NWU_BK'
            P.rspIDX      = 1;
            P.rspFlip     = -1;
            P.hasEEG      = true;
            P.hasMacros   = false;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250723_EEG_NWU_BK(ECGz, beatSep);
    
        case '250912_EEG_NWU_JN'
            P.rspIDX      = 1;
            P.rspFlip     =  1;
            P.hasEEG      = true;
            P.hasMacros   = false;
            % P.spikeClean  = false;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250912_EEG_NWU_JN(ECGz, beatSep);
    
        case '250904_OBE_NWU_TI'
            P.rspIDX      = 3;
            P.rspFlip     = 1;
            P.spikeThresh = 50;
            % P.spikeWin    = 11;
            P.hasEEG      = true;
            % P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250904_OBE_NWU_TI(ECGz, beatSep);
    
        case '250818_Dupi_NMH_JH_2'
            P.rspIDX      = 3;
            P.rspFlip     = 1;
            P.spikeThresh = 20;
            % P.spikeWin    = 11;
            P.hasEEG      = true;
            % P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250818_Dupi_NMH_JH_2(ECGz, beatSep);
    
        case '250811_Dupi_NMH_TPB_1'
            P.rspIDX      = 1;
            P.rspFlip     = 1;
            P.spikeThresh = 10;
            % P.spikeWin    = 9;
            P.hasEEG      = true;
            % P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250811_Dupi_NMH_TPB_1(ECGz, beatSep);
        
        case '250811_Dupi_NMH_TB_2' 
            P.rspIDX      = 1;
            P.rspFlip     = -1;
            P.spikeThresh = 10;
            % P.spikeWin    = 9;
            P.hasEEG      = true;
            % P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250811_Dupi_NMH_TB_2(ECGz, beatSep);

        case '250929_Dupi_NMH_GH_1' 
            P.rspIDX      = 1;
            P.rspFlip     = 1;
            P.spikeThresh = 15;
            % P.spikeWin    = 9;
            P.hasEEG      = true;
            % P.spikeClean  = true;
            P.macroRemove = 6; 
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_250929_Dupi_NMH_GH_1(ECGz, beatSep);


        case '251009_OBE_NWU_CP_1' 
            P.rspIDX      = 1; 
            P.rspFlip     = 1;
            P.macroRemove = 6; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_251009_OBE_NWU_CP_1(ECGz, beatSep);

         case '251002_Dupi_NMH_AB_1' 
            P.rspIDX      = 1; 
            P.rspFlip     = 1;
            P.macroRemove = []; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_251009_OBE_NWU_CP_1(ECGz, beatSep);

        case '251027_Dupi_NMH_DL_1'
            P.rspIDX      = 1; 
            P.rspFlip     = 1;
            P.macroRemove = []; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_251027_Dupi_NMH_DL_1(ECGz, beatSep);
        case '250929_Dupi_NMH_GH_2'
            P.rspIDX      = 3; 
            P.rspFlip     = 1;
            P.macroRemove = [5,6]; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_251027_Dupi_NMH_DL_1(ECGz, beatSep);
        case '251002_Dupi_NMH_AB_2'
            P.rspIDX      = 1; 
            P.rspFlip     = 1;
            P.macroRemove = [5,6]; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_251027_Dupi_NMH_DL_1(ECGz, beatSep);

        case '251013_Dupi_NMH_JN_2'
            P.rspIDX      = 1; 
            P.rspFlip     = 1;
            P.macroRemove = [6]; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.getBeats    = @(ECGz, beatSep)...
                 getBeats_251027_Dupi_NMH_DL_1(ECGz, beatSep);
        otherwise
            % Defaults if you hit an unexpected session:
            error('participant needs parameter specification')
            P.getBeats = @(ECGz, beatSep) error('No heartbeat algorithm defined for session %s', raw.sessID);
    end


end


 % idx = cellfun(@(x) contains(x, 'rsp'), od.labels);
 %    rspDat = od.data(idx,:); 
 %    figure
 %    plot(rspDat(1,:))
 %    hold on 
 %    plot(rspDat(3,:))
   

%% evaluate macros for spike params: 
 % idx = cellfun(@(x) contains(x, 'macro'), od.labels);
 %    figure
 %    macroDat = od.data(idx, :); 
 %    plot(macroDat(1,:))
 %    hold on 
 %    for ii = 2:6
 %        plot(macroDat(ii,:)+(ii-1)*50)
 %    end
 %    legend()
 % 
 %    title([od.sessID ' macros raw'], 'Interpreter','none')





%% ECG function construction: 
% idx = cellfun(@(x) contains(x, 'ECG'), od.labels);
% ECG = od.data(idx, :); 
% 
% d = designfilt('bandpassiir', 'FilterOrder', 4, ...
% 'HalfPowerFrequency1', 5, 'HalfPowerFrequency2', 40, ...
% 'SampleRate', od.fs);
% 
% ECG = filtfilt(d, ECG')'; 
% 
% 
% 
% %plot for custom algorithm design: 
% ECGz = (ECG - mean(ECG, 2)) ./ std(ECG, [], 2); 
% 
% ECGz = ECGz(:, 1:4:end);
% 
% 
% beatSep = od.fs / 20; 
% 
% 
% 
% 
% 
% figure
% plot(ECGz(1,100000:110000), 'color', 'k')
% hold on 
% % plot(ECGz(2,100000:110000), 'color', 'red')
% plot(ECGz(3,100000:110000), 'color', 'green')
% xlim([0 10000])
% xticks([0:1000:10000])
% xticklabels(0:2:20)
% xlabel('Time (s)')
% 
% title(sprintf('ECG beat detection (%s)', ...
%           od.sessID), ...
%             'Interpreter','none');



function heartBeats = getBeats_251027_Dupi_NMH_DL_1(ECGz, beatSep)
    % get within-beat times: JH_1
    test = find(arrayfun(@(x) x > 4, ...
                         ECGz(1,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end



function heartBeats = getBeats_250818_Dupi_NMH_JH_1(ECGz, beatSep)
    % get within-beat times: JH_1
    test = find(arrayfun(@(x,y,z) x > 5 & y > 4 & z < -0.5, ...
                         ECGz(1,3:end), ...
                         ECGz(2,1:end-2), ...
                         ECGz(3,1:end-2)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250818_Dupi_NMH_JH_2(ECGz, beatSep)
    % get within-beat times: JH_2
    test = find(arrayfun(@(x,y,z) x > 3 & y < -3 & z > 1, ...
                         ECGz(2,1:end-13), ...
                         ECGz(3,2:end-12), ...
                         ECGz(3,14:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250623_DUPI_NMH_KS_2(ECGz, beatSep)
    % get within-beat times: KS_2
    test = find(arrayfun(@(x,y,z,n) x > 0.5 & y > 1.75 & z < -2 & n < -1, ...
                         ECGz(1,1:end-9), ...
                         ECGz(2,1:end-9), ...
                         ECGz(3,1:end-9), ...
                         ECGz(2,10:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250623_Dupi_NMH_KS_1(ECGz, beatSep)
    % get within-beat times: KS_1
    test = find(arrayfun(@(x,y,z) x > 2 & y > 0.75 & z < -2, ...
                         ECGz(1,1:end-7), ...
                         ECGz(2,8:end), ...
                         ECGz(3,2:end-6)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250908_OBE_NWU_AS(ECGz, beatSep)
    % get within-beat times: AS
    test = find(arrayfun(@(x,y,z) x > 2 & y > 5 & z < -4, ...
                         ECGz(1,5:end), ...
                         ECGz(2,1:end-4), ...
                         ECGz(3,2:end-3)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250723_EEG_NWU_IN(ECGz, beatSep)
    % get within-beat times: IN
    test = find(arrayfun(@(x,y,z) x < -1 & y > 2 & z < -1, ...
                         ECGz(1,5:end), ...
                         ECGz(2,3:end-2), ...
                         ECGz(3,1:end-4)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250725_EEG_NWU_BN(ECGz, beatSep)
    % get within-beat times: BN
    test = find(arrayfun(@(x,y,z) x > 1 & y < -3 & z > 1, ...
                         ECGz(3,12:end), ...
                         ECGz(3,1:end-11), ...
                         ECGz(2,2:end-10)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250815_EEG_NWU_PP(ECGz, beatSep)
    % get within-beat times: PP
    test = find(arrayfun(@(x,y) x < -4 & y > 3, ...
                         ECGz(1,1:end), ...
                         ECGz(2,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250819_EEG_NWU_ZL(ECGz, beatSep)
    % get within-beat times: ZL
    test = find(arrayfun(@(x,y,z) x < -2 & y > 4 & z < -2, ...
                         ECGz(1,1:end), ...
                         ECGz(2,1:end), ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250723_EEG_NWU_BK(ECGz, beatSep)
    % get within-beat times: BK
    test = find(arrayfun(@(x,y,z) x > 1.5 & y > 2 & z < -3, ...
                         ECGz(1,5:end), ...
                         ECGz(2,1:end-4), ...
                         ECGz(3,3:end-2)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250912_EEG_NWU_JN(ECGz, beatSep)
    % get within-beat times: JN
    test = find(arrayfun(@(x,y) x > 4 & y < -4, ...
                         ECGz(2,1:end), ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250904_OBE_NWU_TI(ECGz, beatSep)
    % get within-beat times: TI
    % (fixed to 2D indexing; original had ECGz(3,3:end, cndi))
    test = find(arrayfun(@(x,y,z) x < -2 & y > 3 & z < 0, ...
                         ECGz(1,1:end-2), ...
                         ECGz(2,2:end-1), ...
                         ECGz(3,3:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250811_Dupi_NMH_TPB_1(ECGz, beatSep)
    % get within-beat times: TPB_1
    % (fixed to 2D indexing; original had ECGz(2,13:end, cndi))
    test = find(arrayfun(@(x,y,z) x > 1 & y < -1 & z < -1, ...
                         ECGz(2,3:end-10), ...
                         ECGz(3,1:end-12), ...
                         ECGz(2,13:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end


function heartBeats = getBeats_250811_Dupi_NMH_TB_2(ECGz, beatSep)
    % get within-beat times: TPB_1
    % (fixed to 2D indexing; original had ECGz(2,13:end, cndi))
    test = find(arrayfun(@(x,y,z) x > 1 & y > 1 & z < -1, ...
                         ECGz(1,12:end), ...
                         ECGz(2,2:end-10), ...
                         ECGz(3,1:end-11)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_250929_Dupi_NMH_GH_1(ECGz, beatSep)
    % get within-beat times: GH 1
    % (fixed to 2D indexing; original had ECGz(2,13:end, cndi))
    test = find(arrayfun(@(x) x < -4 , ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end

function heartBeats = getBeats_251009_OBE_NWU_CP_1(ECGz, beatSep)
    % get within-beat times: CP 1
    % (fixed to 2D indexing; original had ECGz(2,13:end, cndi))
    test = find(arrayfun(@(x) x < -3 , ...
                         ECGz(3,1:end)));
    test = test(diff(test) > beatSep);
    heartBeats = test;
end