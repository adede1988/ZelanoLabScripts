%% Copy all PDF files containing "questionnaire(s)" in the filename
% Crawls participant subfolders under the two study directories and copies
% matching PDFs into the target output directory (organized by study/participant).

roots = {
    'R:\Neurology\Zelano_Lab\Lab_Common\IRB\STU00218720_Recording electrophysiological signals from the human olfactory bulb and olfactory epithelium (OBE)\2025 Participant Documents'
    'R:\Neurology\Zelano_Lab\Lab_Common\IRB\STU00222296_Study of dupilumab mechanism for olfactory restoration (Dupi)\2025 Participant Documents'
    'R:\Neurology\Zelano_Lab\Lab_Common\IRB\STU00222296_Study of dupilumab mechanism for olfactory restoration (Dupi)\2026 Participant Documents'
    };

studyTags = {'OBE','Dupi', 'Dupi'};  % used to prevent naming collisions across studies

outDir = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\questionnaireFiles';

dryRun = false;  % set true to preview without copying

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

results = struct('Source', {}, 'Destination', {}, 'Status', {}, 'Message', {});

for r = 1:numel(roots)
    root = roots{r};
    tag  = studyTags{r};

    if ~exist(root, 'dir')
        warning('Root directory not found: %s', root);
        continue
    end

    % Iterative directory crawl (stack-based) for maximum compatibility
    stack = {root};

    while ~isempty(stack)
        curr = stack{end};
        stack(end) = [];

        listing = dir(curr);

        for i = 1:numel(listing)
            name = listing(i).name;

            % Skip dot dirs
            if listing(i).isdir
                if strcmp(name,'.') || strcmp(name,'..')
                    continue
                end
                stack{end+1} = fullfile(curr, name); %#ok<AGROW>
                continue
            end

            % File handling
            [~, ~, ext] = fileparts(name);

            isPDF = strcmpi(ext, '.pdf');
            hasQuestionnaireWord = contains(lower(name), 'questionnaire'); % matches questionnaire/questionnaires

            if ~(isPDF && hasQuestionnaireWord)
                continue
            end

            src = fullfile(curr, name);

            % Determine participant folder as the first folder below the root
            rel = strrep(src, [root filesep], '');
            parts = strsplit(rel, filesep);
            participant = parts{1};

            % Destination: outDir\StudyTag\Participant\filename.pdf
            destFolder = fullfile(outDir, tag, participant);
            if ~exist(destFolder, 'dir') && ~dryRun
                mkdir(destFolder);
            end

            % Avoid overwriting: if file exists, append _001, _002, ...
            dest = fullfile(destFolder, name);
            if exist(dest, 'file')
                [p, base, e] = fileparts(dest);
                k = 1;
                while true
                    alt = fullfile(p, sprintf('%s_%03d%s', base, k, e));
                    if ~exist(alt, 'file')
                        dest = alt;
                        break
                    end
                    k = k + 1;
                end
            end

            % Copy
            try
                if ~dryRun
                    copyfile(src, dest, 'f');
                end
                results(end+1) = struct( ... %#ok<AGROW>
                    'Source', src, ...
                    'Destination', dest, ...
                    'Status', "OK", ...
                    'Message', "" ...
                );
                fprintf('[%s] %s -> %s\n', tag, src, dest);
            catch ME
                results(end+1) = struct( ... %#ok<AGROW>
                    'Source', src, ...
                    'Destination', dest, ...
                    'Status', "FAIL", ...
                    'Message', string(ME.message) ...
                );
                warning('Failed to copy: %s\n  Reason: %s', src, ME.message);
            end
        end
    end
end

% Summary
okCount   = sum(string({results.Status}) == "OK");
failCount = sum(string({results.Status}) == "FAIL");
fprintf('\nDone. Copied: %d | Failed: %d | Output: %s\n', okCount, failCount, outDir);

% Optional: write a log CSV
% T = struct2table(results);
% writetable(T, fullfile(outDir, 'questionnaire_copy_log.csv'));
