set(0, 'defaultfigurewindowstyle', 'docked')

figure
plot(smoothdata(-chanDat.data(36,548590:933250), 'gaussian', 500))
hold on
plot(smoothdata(-chanDat.data(36,1319861:1703931), 'gaussian', 500))

figure; 
histogram(chanDat.behDat.length(chanDat.behDat.goodBreath==1 & ...
                                strcmp(chanDat.behDat.task,'slowFocus')), ...
                                [0:1:14])
hold on 

histogram(chanDat.behDat.length(chanDat.behDat.goodBreath==1& ...
                                strcmp(chanDat.behDat.task,'fastFocus')), ...
                                [0:1:14])