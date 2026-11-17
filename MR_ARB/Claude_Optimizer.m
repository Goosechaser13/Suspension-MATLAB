%% Suspension MR optimization: tune bellcrank spring + pushrod points to hit a TARGET motion ratio
clear; clc; close all;

%% ---------------- TARGET ----------------
targetMR = 2;    % desired motion ratio (damper travel / wheel travel). Change me.

%% Fixed geometry (in.) Can paste param. from SuspensionAnalysis.m
MSpring    = [10.1;15.11;0];
MBell      = [11.5;7.92;0];
MPush      = [21.72;4.21;0];
MBellSpr  = [16.24;11.76;0];
MBellPush = [16.33;9.56;0];
MArmL      = [10.35;2.75;0];
MArmU      = [12;7.67;0];
MArmL_End  = MPush;
MArmU_End  = [21.72;6.62;0];

x_rest = 7; x_min = 5.94; x_max = 8.19; % Damper ranges
n = 40; % test resolution
xDamper = linspace(x_min,x_max,n);

% Package fixed data for the solver
G.MSpring = MSpring; G.MBell = MBell; G.MPush = MPush;
G.MArmL = MArmL; G.MArmU = MArmU;
G.MArmL_End = MArmL_End; G.MArmU_End = MArmU_End;
G.x = xDamper; G.x_rest = x_rest;
G.target = targetMR;

%% Optimization setup
% design variables: [BellSpr_x BellSpr_y BellPush_x BellPush_y]
p0 = [MBellSpr(1:2); MBellPush(1:2)];
searchRadius = 3;          % in. each point may move +/- this much
G.minClear = 0.5;          % min perpendicular distance of pushrod point below the line (in.)
lb = p0 - searchRadius;
ub = p0 + searchRadius;

w_mean = 1;      % weight on (mean(MR) - target)^2   -> hit the target on average
w_var  = 1;      % weight on var(MR)                 -> keep MR flat across travel
obj   = @(p) objective(p,G,w_mean,w_var);
nlcon = @(p) constraints(p,G);

opts = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'MaxIterations',500,'MaxFunctionEvaluations',5000, ...
    'StepTolerance',1e-10,'OptimalityTolerance',1e-8);

%% Sanity check: local solver must match DamperX_to_Points
if exist('DamperX_to_Points','file') == 2
    try
        pts0 = [G.MSpring, G.MBell, G.MPush, [p0(3:4);0], [p0(1:2);0], ...
                G.MArmL, G.MArmU, G.MArmL_End, G.MArmU_End];
        maxErr = 0;
        for xi = G.x(1:3:end)
            a = damperToPoints(pts0,xi);
            b = DamperX_to_Points(pts0,xi);
            maxErr = max(maxErr, max(abs(a(:)-b(:))));
        end
        fprintf('Solver check vs DamperX_to_Points: max point difference = %.2e in.\n',maxErr);
        if maxErr > 1e-6
            warning('Local solver does not match DamperX_to_Points!');
        end
    catch ME
        fprintf('Solver check skipped (%s)\n',ME.message);
    end
else
    fprintf('Solver check skipped (DamperX_to_Points.m not on path).\n');
end

%% Multi-start (first start = your current design)
nStarts = 25;
rng(1);
starts = [p0, p0 + (rand(4,nStarts-1)-0.5)*2*searchRadius*0.5];
% put every start on the rest-length constraint (spring point x_rest from the chassis mount)
for k = 1:nStarts
    v = starts(1:2,k) - MSpring(1:2);
    starts(1:2,k) = MSpring(1:2) + x_rest*v/norm(v);
end
starts = min(max(starts,lb),ub);

bestJ = inf; bestP = p0;
for k = 1:nStarts
    try
        [pk,Jk,flag] = fmincon(obj,starts(:,k),[],[],[],[],lb,ub,nlcon,opts);
        [~,~,okk] = solveMR(pk,G);
        if flag > 0 && okk && Jk < bestJ
            bestJ = Jk; bestP = pk;
        end
    catch
        % infeasible/unsolvable start, skip
    end
end

%% Results
fprintf('\nTarget MR           : %.3f\n',targetMR);
fprintf('Original BellSpring : [%.3f, %.3f]\n',p0(1),p0(2));
fprintf('Original BellPush   : [%.3f, %.3f]\n',p0(3),p0(4));
fprintf('Optimized BellSpring: [%.3f, %.3f]\n',bestP(1),bestP(2));
fprintf('Optimized BellPush  : [%.3f, %.3f]\n',bestP(3),bestP(4));

[MR0,wh0,ok0] = solveMR(p0,G);
[MR1,wh1,ok1] = solveMR(bestP,G);
if ok0, fprintf('Original : mean MR = %.4f, std = %.4f, RMS error vs target = %.4f\n', ...
        mean(MR0),std(MR0),sqrt(mean((MR0-targetMR).^2))); end
if ok1, fprintf('Optimized: mean MR = %.4f, std = %.4f, RMS error vs target = %.4f\n', ...
        mean(MR1),std(MR1),sqrt(mean((MR1-targetMR).^2))); end
fprintf('Pushrod length: %.3f -> %.3f in\n', ...
    norm(MPush(1:2)-p0(3:4)), norm(MPush(1:2)-bestP(3:4)));
fprintf('Spring length at rest (optimized): %.4f in\n', norm(MSpring(1:2)-bestP(1:2)));

figure;
if ok0, plot(wh0,MR0,'b-','LineWidth',1.5); hold on; end
if ok1, plot(wh1,MR1,'r-','LineWidth',1.5); end
yline(targetMR,'k--');
xlabel('Wheel travel from rest (in.)'); ylabel('Motion ratio (damper/wheel)');
legend('Original','Optimized',sprintf('Target = %.2f',targetMR)); grid on;
title('Instantaneous motion ratio');

%% ---------------- Local functions ----------------
function J = objective(p,G,w_mean,w_var)
    [MR,~,ok,nBad] = solveMR(p,G);
    if ~ok
        % graded penalty (more unreachable/non-monotonic samples = worse) so SQP/multistart can still rank bad points
        J = 1e3*(1 + nBad/numel(G.x)); return;
    end
    J = w_mean*(mean(MR)-G.target)^2 + w_var*var(MR);
end

function [c,ceq] = constraints(p,G)
spr  = p(1:2);          % BellSpring point
push = p(3:4);          % BellPush point
piv  = G.MBell(1:2);    % bellcrank pivot

d = spr - piv;          % line from pivot to spring point
r = push - piv;
crossZ = d(1)*r(2) - d(2)*r(1);     % negative => pushrod point is below/right of the line
signedDist = crossZ / norm(d);      % perpendicular distance (negative = below)

% Pushrod point must be at least minClear below the line:  signedDist <= -minClear
c = signedDist + G.minClear;

% Spring length at rest must equal x_rest
ceq = norm(G.MSpring(1:2) - spr) - G.x_rest;
end

function [MR,wheelH,ok,nBad] = solveMR(p,G)
    % Same procedure as SuspensionAnalysis.m: sweep damper length, solve linkage with
    % DamperX_to_Points logic, take wheel height (row 2, point 10), MR = d(damper)/d(wheel).
    ok = false; MR = []; wheelH = []; nBad = 0;
    pts = [G.MSpring, G.MBell, G.MPush, [p(3:4);0], [p(1:2);0], ...
           G.MArmL, G.MArmU, G.MArmL_End, G.MArmU_End];
    n = numel(G.x);
    wh = nan(n,1);
    for i = 1:n
        try
            np = damperToPoints(pts,G.x(i));
            wh(i) = np(2,10);
        catch
            % linkage can't reach this damper length
        end
    end
    nBad = sum(isnan(wh));
    if nBad > 0, return; end
    try
        whRest = damperToPoints(pts,G.x_rest);
    catch
        nBad = n; return;
    end
    wh = wh - whRest(2,10);          % wheel height relative to rest (constant offset; doesn't change MR)
    dw = diff(wh);
    if ~(all(dw>1e-6) || all(dw<-1e-6))   % wheel height must be monotonic
        nBad = 1 + sum(sign(dw) ~= sign(median(dw)));
        return;
    end
    % Same formula as SuspensionAnalysis (gradient(x,wheelHeight)); abs() only removes the sign
    % convention (wheel up = damper compress gives a negative value there).
    MR = abs(gradient(G.x(:),wh));
    wheelH = wh;
    ok = true;
end

function newPoints = damperToPoints(points,x)
    % Mirrors DamperX_to_Points.m without needing circcirc (Mapping TB) or rotz (Robotics TB)
    newPoints = points;
    % spring-side bellcrank point
    L_bs = norm(points(:,5)-points(:,2));
    [xo,yo] = circInt(points(1,1),points(2,1),x,points(1,2),points(2,2),L_bs);
    [newPoints(1,5),newPoints(2,5)] = pickNearest(xo,yo,points(1,5),points(2,5));
    % bellcrank rotation
    v_old = points(1:2,5)-points(1:2,2);
    v_new = newPoints(1:2,5)-points(1:2,2);
    dA = atan2d(v_new(2),v_new(1)) - atan2d(v_old(2),v_old(1));
    R = [cosd(dA) -sind(dA) 0; sind(dA) cosd(dA) 0; 0 0 1];
    newPoints(:,4) = points(:,2) + R*(points(:,4)-points(:,2));
    % pushrod / lower arm end
    L_arm = norm(points(:,3)-points(:,6));
    L_pr  = norm(points(:,3)-points(:,4));
    [xo,yo] = circInt(newPoints(1,4),newPoints(2,4),L_pr,points(1,6),points(2,6),L_arm);
    [newPoints(1,3),newPoints(2,3)] = pickNearest(xo,yo,points(1,3),points(2,3));
    newPoints(:,8) = newPoints(:,3);
    % upper arm end
    L_kn = norm(points(:,9)-points(:,8));
    L_au = norm(points(:,9)-points(:,7));
    [xo,yo] = circInt(newPoints(1,8),newPoints(2,8),L_kn,points(1,7),points(2,7),L_au);
    [newPoints(1,9),newPoints(2,9)] = pickNearest(xo,yo,points(1,9),points(2,9));
    v_knuckle = (newPoints(:,9) - newPoints(:,8)) / norm(newPoints(:,9) - newPoints(:,8));
    newPoints(:,10) = v_knuckle*0.97 + newPoints(:,8);
end

function [xo,yo] = circInt(x1,y1,r1,x2,y2,r2)
    d = hypot(x2-x1,y2-y1);
    if d > r1+r2 || d < abs(r1-r2) || d == 0
        xo = []; yo = []; return;
    end
    a = (r1^2 - r2^2 + d^2)/(2*d);
    h = sqrt(max(r1^2 - a^2,0));
    xm = x1 + a*(x2-x1)/d;  ym = y1 + a*(y2-y1)/d;
    xo = [xm + h*(y2-y1)/d, xm - h*(y2-y1)/d];
    yo = [ym - h*(x2-x1)/d, ym + h*(x2-x1)/d];
end

function [px,py] = pickNearest(xo,yo,x0,y0)
    if isempty(xo), error('No intersection'); end
    [~,k] = min(hypot(xo-x0,yo-y0));
    px = xo(k); py = yo(k);
end