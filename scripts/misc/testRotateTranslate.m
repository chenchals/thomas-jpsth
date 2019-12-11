% From: https://stackoverflow.com/questions/29780903/2d-body-transformation-and-rotation-in-matlab
% Define A
A = [-2,-2,6,6,-2; -2,2,2,-2,-2; 1 1 1 1 1];
% Define Translation Matrix
trans = @(x,y,z) repmat([x; y; z],[1 5]);
% Define Rotation Matrix
se2 = @(x, y, theta) [
    cosd(theta), -sind(theta), x;
    sind(theta), cosd(theta), y;
    0,        0,           1];
% Calculate Rotated Rect
B = se2(0,0,45) * (A - trans(2,0,0) ) + trans(2,0,0);
% Plot Rectangles
figure; plot(A(1,:),A(2,:),'b')
hold on;
plot(B(1,:),B(2,:),'r')
hold off;
axis equal

xvals = [1:10];
yvals = zeros(1,10);
zvals = ones(1,10);

transOrig = @(x,y,z) repmat([x;y;z;],[1,numel(x)]);
rotOrig =  @(x, y, theta) [
    cosd(theta), -sind(theta), x;
    sind(theta), cosd(theta), y;
    0,            0,           1];

rotMat = rotOrig(0,0,45) * ...
    ([xvals;yvals;zvals] - transOrig(floor(numel(xvals)/2),0,0) ) + transOrig(floor(numel(xvals)/2),0,0);


hold off

figure
plot(rotMat(1,:),rotMat(2,:),'b')



rotOrig(0,0,45)*[xvals;yvals;zvals]



