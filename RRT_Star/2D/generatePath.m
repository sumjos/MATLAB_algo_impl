% RRT* algorithm in 2D with collision avoidance.
%
% Author: Sumit Joshi
%
% nodes:    Contains list of all explored nodes. Each node contains its
%           coordinates, cost to reach and its parent.
%
% Brief description of algorithm:
% 1. Pick a random node q_rand.
% 2. Find the closest node q_near from explored nodes to branch out from, towards
%    q_rand.
% 3. Steer from q_near towards q_rand: interpolate if node is too far away, reach
%    q_new. Check that obstacle is not hit.
% 4. Update cost of reaching q_new from q_near, treat it as Cmin. For now,
%    q_near acts as the parent node of q_new.
% 5. From the list of 'visited' nodes, check for nearest neighbors with a
%    given radius, insert in a list q_nearest.
% 6. In all members of q_nearest, check if q_new can be reached from a
%    different parent node with cost lower than Cmin, and without colliding
%    with the obstacle. Select the node that results in the least cost and
%    update the parent of q_new.
% 7. Add q_new to node list.
% 8. Continue until maximum number of nodes is reached or goal is hit.

function pathNodes =  generatePath()
clearvars
close all
%initialization
global minXIndex;
global trackInterp;
global TestTrack;
global Xo;
global Yo;
global Xf;
global Yf;
global safeDistance;
global r
global lookupDistance;
global forwardBiasPercent;
global intervalForLineIntersectionCheck;
global maxPermittedAngle;

%Tuning paramerters
lookupDistance = 200;
forwardBiasPercent = 0.6;
safeDistance = 0.1;
intervalForLineIntersectionCheck = 0.5;
r = 60;
maxPermittedAngle = 87;
EPS = 500;

numNodes = 3000;
TestTrack = load('TestTrack.mat');
%replace later by inputs
TestTrack = TestTrack.TestTrack;
% obs = {[1122.76376860723,488.631892695626;
%     1123.67026629858,484.488818920667;
%     1124.75187949990,484.725474100659;
%     1123.84538180854,488.868547875618],...
%     [1377.72095614637,701.071605966588;
%     1382.74280526035,696.843352410747;
%     1384.03227680587,698.374843238537;
%     1379.01042769188,702.603096794378]};
obs = generateRandomObstacles(4, TestTrack);
minXIndex = 12;
%load intial and final coocrdinates
Xo = 287;
Yo = -176;
Xf = TestTrack.cline(1, end);
Yf = TestTrack.cline(2, end);

%create and save x,y interpolation
interpxy(TestTrack.bl(1,:), TestTrack.bl(2,:), TestTrack.br(1,:), TestTrack.br(2,:));
trackInterp = load('interpXY');

%obstacle = [500,150,200,200];
obstacle = obs;
q_start.coord = [Xo Yo];
q_start.cost = 0;
q_start.parent = 0;
q_goal.coord = [Xf Yf];
q_goal.cost = 0;

nodes(1) = q_start;


figure(1) 
axis([TestTrack.bl(1, minXIndex) - 10 TestTrack.br(1, end) + 10 ...
    TestTrack.bl(2, 1) - 10 TestTrack.bl(2, end) + 10])
hold on
plotTrack();
plotObstacles(obstacle);

i = 0;
counter = 0;
pathCompleted = false;
finishLine = @(x) TestTrack.bl(2, end) + ((TestTrack.bl(2, end) - TestTrack.br(2, end))...
    *(x - TestTrack.bl(1, end)) / (TestTrack.bl(1, end) - TestTrack.br(1, end)));
Xfinish = TestTrack.bl(1, end) : 0.1 : TestTrack.br(1, end);
Yfinish = cell2mat(arrayfun(finishLine, Xfinish,  'UniformOutput', false));
plot(Xfinish, Yfinish, '-k', 'LineWidth', 1);
while i <= numNodes && ~pathCompleted
    i = i + 1;
    %biasing function
    counter = counter + 1;
    q_rand = getBiasedQrand(nodes);
    if ~checkIfPointInsideTrack(q_rand(1), q_rand(2))
        i = i - 1;
        continue;
    end
    plot(q_rand(1), q_rand(2), 'x', 'Color',  [0 0.4470 0.7410]);
    % Break if node is already reached or has crossed line
    for j = 1:1:length(nodes)
        y_calc = finishLine(nodes(j).coord(1));
        y_act = nodes(j).coord(2);
        if nodes(j).coord(2) > finishLine(nodes(j).coord(1))
            pathCompleted = true;
            break;
        end
    end
    
    % Pick the closest node from existing list to branch out from
    ndist_temp = [];
    

%% Single minimum point check
    for j = 1:1:length(nodes)
        n = nodes(j);
        tmp = dist(n.coord, q_rand);
        ndist_temp = [ndist_temp tmp];
    end
    
    %sorting ndist
    [val, minId] = min(ndist_temp);
     q_near = nodes(minId);
        q_new.coord = steer(q_rand, q_near.coord, val, EPS);
%         q_new.parent = NaN;
    tatti = noCollision(q_rand, q_near.coord, obstacle);
    if tatti == 1

        
        line([q_near.coord(1), q_new.coord(1)], [q_near.coord(2), q_new.coord(2)], 'Color', 'k', 'LineWidth', 1);
        drawnow
        hold on
        q_new.cost = dist(q_new.coord, q_near.coord) + q_near.cost;
        
        % Within a radius of r, find all existing nodes
        q_nearest = [];
        neighbor_count = 1;
        for j = 1:1:length(nodes)
            if noCollision(q_new.coord, nodes(j).coord, obstacle)==1 && dist(nodes(j).coord, q_new.coord) <= r
                q_nearest(neighbor_count).coord = nodes(j).coord;
                q_nearest(neighbor_count).cost = nodes(j).cost;
                neighbor_count = neighbor_count+1;
            end
        end
        
        % Initialize cost to currently known value
        q_min = q_near;
        C_min = q_new.cost;
        
        % Iterate through all nearest neighbors to find alternate lower
        % cost paths
        
        for k = 1:1:length(q_nearest)
            if ((q_nearest(k).cost + dist(q_nearest(k).coord, q_new.coord)) < C_min)
                q_min = q_nearest(k);
                C_min = q_nearest(k).cost + dist(q_nearest(k).coord, q_new.coord);
                line([q_min.coord(1), q_new.coord(1)], [q_min.coord(2), q_new.coord(2)], 'Color', 'g');
                hold on
            end
        end
        q_new.cost = C_min;
        
        % Update parent to least cost-from node
        for j = 1:1:length(nodes)
            if nodes(j).coord == q_min.coord
                q_new.parent = j;
            end
        end
        
        % Append to nodes
        nodes = [nodes q_new];
        
        % Check if nodes in vicinty can be traversed at lower cost through
        % q_new
        for k = 1:1:length(q_nearest)
            c1 = q_new.cost + dist(q_nearest(k).coord, q_new.coord);
            c2 = q_nearest(k).cost;
            if q_new.cost + dist(q_nearest(k).coord, q_new.coord) < q_nearest(k).cost
                % Update parent to least cost-from node
                for j = 1:1:length(nodes)
                    if nodes(j).coord == q_nearest(k).coord
                        % since last new node was appended to end
                        nodes(j).parent = length(nodes);
                        break;
                    end
                end
            end
        end
    end

end

D = [];
for j = 1:1:length(nodes)
    tmpdist = dist(nodes(j).coord, q_goal.coord);
    D = [D tmpdist];
end

% Search backwards from goal to start to find the optimal least cost path
[val, idx] = min(D);
q_final = nodes(idx);
q_goal.parent = idx;
q_end = q_goal;
nodes = [nodes q_goal];
figure(2)
axis([TestTrack.bl(1, minXIndex) - 10 TestTrack.br(1, end) + 10 ...
    TestTrack.bl(2, 1) - 10 TestTrack.bl(2, end) + 10])
hold on
plotTrack();
plot(Xfinish, Yfinish, '-k', 'LineWidth', 1);
plotObstacles(obstacle);
pathNodes{1} = q_end.coord;
pathIdx = 2;
while q_end.parent ~= 0
    start = q_end.parent;
    line([q_end.coord(1), nodes(start).coord(1)], [q_end.coord(2), nodes(start).coord(2)], 'Color', 'r', 'LineWidth', 1);
    q_end = nodes(start);
    pathNodes{pathIdx} = q_end.coord;
    pathIdx = pathIdx + 1;
end

end


%%%%%%%%%%%%%
%% New Interpolation function with end extension
function interpxy(Xl, Yl, Xr, Yr)
global minXIndex;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Extension at the end
ext_dist=100;

%Left line
slope_l_end=(Yl(end)-Yl(end-1))/(Xl(end)-Xl(end-1));
Xl_new=Xl(end)+ext_dist*cos(atan(slope_l_end));
Yl_new=Yl(end)+ext_dist*sin(atan(slope_l_end));

%Right line
slope_r_end=(Yr(end)-Yr(end-1))/(Xr(end)-Xr(end-1));
Xr_new=Xr(end)+ext_dist*cos(atan(slope_r_end));
Yr_new=Yr(end)+ext_dist*sin(atan(slope_r_end));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%linear interpolation strategy
%left track boundary interp
Xl1 = Xl(1:minXIndex);
Yl1 = Yl(1:minXIndex);
Xl2 = [Xl(minXIndex:end) Xl_new];
Yl2 = [Yl(minXIndex:end) Yl_new];
interpolLeft1 = @(x) interp1(Xl1, Yl1,x);
interpolLeft2 = @(x) interp1(Xl2, Yl2, x);

%left track boundary interp
Xr1 = [Xr(1:minXIndex) Xl(minXIndex)];
Yr1 = [Yr(1:minXIndex) Yr(minXIndex)];
Xr2 = [Xr(minXIndex:end) Xr_new];
Yr2 = [Yr(minXIndex:end) Yr_new];
interpolRight1 = @(x) interp1(Xr1, Yr1,x);
interpolRight2 = @(x) interp1(Xr2, Yr2, x);

save('interpXY', 'interpolLeft1', 'interpolLeft2', ...
    'interpolRight1', 'interpolRight2');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function isInsideTrack = checkIfPointInsideTrack(x, y)
global minXIndex;
global trackInterp;
global TestTrack;
global safeDistance;
isInsideTrack = true;
%check in lower and upper regions of track
if y <= TestTrack.br(2, minXIndex)
    y_calc_bl = trackInterp.interpolLeft1(x);
    y_calc_br = trackInterp.interpolRight1(x);
    if isnan(y_calc_br) || y_calc_br < (y + safeDistance)...
            || (~isnan(y_calc_bl) && (y - safeDistance) < y_calc_bl)
        isInsideTrack = false;
        return;
    end
else
    y_calc_bl = trackInterp.interpolLeft2(x);
    y_calc_br = trackInterp.interpolRight2(x);
    if (isnan(y_calc_bl) || ((y + safeDistance) > y_calc_bl)...
            || ((~isnan(y_calc_br) && ((y - safeDistance) < y_calc_br))))
        isInsideTrack = false;
        return;
    end
end
end

function plotTrack()
global TestTrack;
global minXIndex;
global trackInterp;
%track lower section
Xl1 = linspace(TestTrack.bl(1, 1), TestTrack.bl(1, minXIndex), 2000);
Yl1 = cell2mat(arrayfun(trackInterp.interpolLeft1, Xl1, 'UniformOutput', false));
Xr1 = linspace(TestTrack.br(1, 1), TestTrack.br(1, minXIndex), 2000);
Yr1 =  cell2mat(arrayfun(trackInterp.interpolRight1, Xr1, 'UniformOutput', false));
%track upper section
Xl2 = linspace(TestTrack.bl(1, minXIndex), TestTrack.bl(1, end), 2000);
Yl2 = cell2mat(arrayfun(trackInterp.interpolLeft2, Xl2, 'UniformOutput', false));
Xr2 = linspace(TestTrack.br(1, minXIndex), TestTrack.br(1, end), 2000);
Yr2 = cell2mat(arrayfun(trackInterp.interpolRight2, Xr2, 'UniformOutput', false));
plot(Xl1, Yl1, 'k', Xr1, Yr1, 'k', Xl2, Yl2, 'k', Xr2, Yr2, 'k');

end

function plotObstacles(obstacles, removePrev)
global obsFigHandle;
for idx = 1:length(obstacles)
    obs = obstacles{idx};
    E1x = obs(1:2,1);
    E1y = obs(1:2,2);
    E2x = obs(2:3,1);
    E2y = obs(2:3,2);
    E3x = obs(3:4,1);
    E3y = obs(3:4,2);
    E4x = [obs(1,1) obs(4,1)];
    E4y = [obs(1,2) obs(4,2)];
    if nargin == 2 && removePrev
        set(obsFigHandle, 'Visible', 'off');
    end
    obsFigHandle = plot(E1x, E1y, 'k', E2x, E2y, 'k', E3x, E3y, 'k', E4x, E4y, 'k');
    
end
end

function nc = noCollision(n2, n1, obstacles)
A = [n1(1) n1(2)];
B = [n2(1) n2(2)];
nc = true;
if isLineIntersectingTrack(A, B)
    nc = false;
    return;
end
for idx = 1:length(obstacles)
    obs = obstacles{idx};
    
    C1 = [obs(1,1),obs(1,2)];
    D1 = [obs(4,1),obs(4,2)];
    C2 = [obs(1,1),obs(1,2)];
    D2 = [obs(2,1),obs(2,2)];
    C3 = [obs(3,1),obs(3,2)];
    D3 = [obs(2,1),obs(2,2)];
    C4 = [obs(3,1),obs(3,2)];
    D4 = [obs(4,1),obs(4,2)];
    
    % Check if path from n1 to n2 intersects any of the four edges of the
    % obstacle
    
    ints1 = ccw(A,C1,D1) ~= ccw(B,C1,D1) && ccw(A,B,C1) ~= ccw(A,B,D1);
    ints2 = ccw(A,C2,D2) ~= ccw(B,C2,D2) && ccw(A,B,C2) ~= ccw(A,B,D2);
    ints3 = ccw(A,C3,D3) ~= ccw(B,C3,D3) && ccw(A,B,C3) ~= ccw(A,B,D3);
    ints4 = ccw(A,C4,D4) ~= ccw(B,C4,D4) && ccw(A,B,C4) ~= ccw(A,B,D4);
    if ints1==0 && ints2==0 && ints3==0 && ints4==0
        nc = true;
    else
        nc = false;
        break;
    end
end
end

function isIntersecting = isLineIntersectingTrack(A, B)
global intervalForLineIntersectionCheck;
global maxPermittedAngle;
isIntersecting = false;
step = intervalForLineIntersectionCheck;
if (A(1) > B(1))
    step = -step;
end
if atand(abs(A(2) - B(2))/abs((A(1) - B(1)))) > maxPermittedAngle
    isIntersecting = true;
    return;
end

%disp(atand(abs(A(2) - B(2))/(A(1) - B(1))));
Yfun = @(x) A(2) + (A(2) - B(2))*(x - A(1))/(A(1) - B(1));
X = A(1) : step : B(1);
for idx = 1 : length(X)
    y = Yfun(X(idx));
    if ~checkIfPointInsideTrack(X(idx), y)
        isIntersecting = true;
        return;
    end
end
end

function val = ccw(A,B,C)
val = (C(2)-A(2)) * (B(1)-A(1)) > (B(2)-A(2)) * (C(1)-A(1));
end

function Qrand = getBiasedQrand(nodes)
global Yo;
global TestTrack;
global minXIndex;
global lookupDistance;
global forwardBiasPercent;

X_absMin =  TestTrack.bl(1, minXIndex);
%Xmax =  TestTrack.br(1,end) + 2;
Y_absMin =  Yo;
% Ymax =  TestTrack.bl(2,end) + 2;
%q_rand = [Xmin + floor(rand(1)*(Xmax - Xmin)) Ymin + floor(rand(1)*(Ymax - Ymin))];

%baising logic: is used to generate bluster of random points in front of vehicle
%Can modify this strategy to make it better
q = nodes(1);
for idx = 1:length(nodes)
    
    %Max cost strategy
%         if (nodes(idx).cost > q.cost)
%             q = nodes(1,idx);
%         end
    
    % ax+by sum strategy, where a and b are weights on coordinates
    a = 0.7; b = 1-a;
    if ((a*nodes(idx).coord(1) + b*nodes(idx).coord(2)) > ...
            (a*q.coord(1) + b*q.coord(2)))
        q = nodes(idx);
    end
end

Xmin = q.coord(1) - (1- forwardBiasPercent)*lookupDistance;
Ymin = q.coord(2) - (1- forwardBiasPercent)*lookupDistance;
if Xmin < X_absMin
    Xmin = X_absMin;
end
if Ymin < Y_absMin
    Ymin = Y_absMin;
end
%roi of rand
roi = {[Xmin, Ymin; Xmin + lookupDistance, Ymin; ...
    Xmin + lookupDistance, Ymin + lookupDistance;
    Xmin, Ymin + lookupDistance]};
%plotObstacles(roi, true);
Qrand = [Xmin + floor(rand(1)*lookupDistance) ...
    Ymin + floor(rand(1)*lookupDistance)];

end

function A = steer(qr, qn, val, eps)
qnew = [0 0];

% Steer towards qn with maximum step size of eps
if val >= eps
    qnew(1) = qn(1) + ((qr(1)-qn(1))*eps)/dist(qr,qn);
    qnew(2) = qn(2) + ((qr(2)-qn(2))*eps)/dist(qr,qn);
else
    qnew(1) = qr(1);
    qnew(2) = qr(2);
end
A = [qnew(1), qnew(2)];
end

function d = dist(q1,q2)
d = sqrt((q1(1)-q2(1))^2 + (q1(2)-q2(2))^2);
end
