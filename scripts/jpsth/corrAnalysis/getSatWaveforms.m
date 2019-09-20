jpsthPairs = load('dataProcessed/dataset/JPSTH_PAIRS_CellInfoDB.mat');
jpsthPairs = jpsthPairs.JpsthPairCellInfoDB;
ninfo = load('dataProcessed/dataset/ninfo_nstats_SAT.mat');
ninfo = struct2table(ninfo.ninfo);
%Choice Errors (n=24): find( [ninfo.errGrade] >= 2)
eChoice = ninfo(ninfo.errGrade >= 2,:);
jpsthPairs.tempUnitNum = jpsthPairs.X_unitNum;
A = innerjoin(eChoice,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
jpsthPairs.tempUnitNum = jpsthPairs.Y_unitNum;
B = innerjoin(eChoice,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
eChoicePairs = sortrows([A;B],'unitNum');
clearvars A B

%Timing Errors (n = 25): find( abs([ninfo.rewGrade]) >= 2)
eTiming = ninfo(abs(ninfo.rewGrade) >= 2,:);
jpsthPairs.tempUnitNum = jpsthPairs.X_unitNum;
A = innerjoin(eTiming,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
jpsthPairs.tempUnitNum = jpsthPairs.Y_unitNum;
B = innerjoin(eTiming,jpsthPairs,'LeftKeys','unitNum','RightKeys',{'tempUnitNum'});
eTtmingPairs = sortrows([A;B],'unitNum');

clearvars A B jpsthPairs

