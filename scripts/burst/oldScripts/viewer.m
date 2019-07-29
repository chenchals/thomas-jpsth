var = oVar25;
% Build augmented jpsth matrix for all pairs
border = 2;
% Pairs are [1,2];[1,3];[1,4];....[n-1,n]
pairs = var.cellPairs;
% fliplr, so that when creating augmented matrix 1st column is cell1 on X
% and cells 2 to n on Y
pairs = fliplr(pairs);
maxPairsPerCell = sum(pairs(:,1)==1);
jpsthSide = size(var.jpsth(1).normJpsth,1);
augmentedJpsth = nan(maxPairsPerCell*(jpsthSide+border));% a square matrix
augmentedByCell = {};
for i = 1: size(pairs,1)
    x = (pairs(i,1)-2)*(jpsthSide+border);%X-Cell
    y = (pairs(i,2)-1)*(jpsthSide+border);%Y-Cell
    augmentedJpsth(x+1:x+jpsthSide,y+1:y+jpsthSide) = flipud(var.jpsth(i).normJpsth);
end
figure()
cLims = [nanmin(augmentedJpsth(:)) nanmax(augmentedJpsth(:))];
alpha = ones(size(augmentedJpsth));
alpha(isnan(augmentedJpsth)) = 0;
im = imagesc(augmentedJpsth,cLims);
colormap(gca);
im.AlphaData = alpha;


f=figure()
p(1,1)=subplot(29,1,1)
p(1,2)=subplot(29,1,2)


% figure()
% for i=1:406
%     v=oVar25.jpsth(i).coinHist(:,2); 
%     v=oVar25.smoothFx(v,19,3);
%     if max(v)>0.25
%         i, 
%         plot(v); 
%         hold on;
%         drawnow;
%     end
% end
% hold off