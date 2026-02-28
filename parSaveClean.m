function [] = parSaveClean(fn, cleanDat)
    save(fn, 'cleanDat', '-v7.3')
end