ninfo = struct2table(NINFOS.ninfo);
ninfo.rewGrade = abs(ninfo.rewGrade);

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



