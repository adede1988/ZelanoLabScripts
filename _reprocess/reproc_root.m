function r = reproc_root()
% Backup/log root: E:\reprocBackup on the lab desktop, else <repo>\..\reprocBackup.
    if isfolder('E:\')
        r = 'E:\reprocBackup';
    else
        L = labPaths(); r = fullfile(fileparts(L.repo), 'reprocBackup');
    end
    if ~isfolder(r), mkdir(r); end
end
