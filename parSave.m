function [] = parSave(fn, chanDat)
    save(fn, 'chanDat', '-v7.3')
end