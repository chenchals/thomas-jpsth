ninfo = struct2table(NINFOS.ninfo);
ninfo.rewGrade = abs(ninfo.rewGrade);

% From plotSDF_ErrorCompare_SAT.m
% idxErr = ([ninfo.errGrade] >= 2);
% idxRew = (abs([ninfo.rewGrade]) >= 2 & ~isnan([nstats.A_Reward_tErrStart_Fast]));


dninfo = ninfo(strcmp(ninfo.monkey,'D'),:);
eninfo = ninfo(strcmp(ninfo.monkey,'E'),:);


% choice err only
% errGrade >=2;
dChoice = sum(dninfo.errGrade>=2 & ~(dninfo.rewGrade>=2))
eChoice = sum(eninfo.errGrade>=2 & ~(eninfo.rewGrade>=2))

% both
dBoth = sum(dninfo.errGrade>=2 & dninfo.rewGrade>=2)
eBoth = sum(eninfo.errGrade>=2 & eninfo.rewGrade>=2)

% timing only
dTime = sum(~(dninfo.errGrade>=2) & dninfo.rewGrade>=2)
eTime = sum(~(eninfo.errGrade>=2) & eninfo.rewGrade>=2)



