%% Claude_ARB_Optimizer.m  (v2 - flattens MR_ARB over the whole travel)
% Tunes the ARB layout so the ARB motion ratio
%     MR_ARB = d(ARB end height) / d(wheel height)
% stays as close as possible to targetMR over the ENTIRE damper sweep.
% MR_ARB is computed exactly like SuspensionAnalysis.m (DamperX_to_Points +
% gradient of point 13's Y vs wheel height).
%
% Needs: DamperX_to_Points.m on the path (with the branch-selection fix for
% newPoints(:,13)), Mapping Toolbox (circcirc, wrapToPi).
% Optimization Toolbox is used if present (lsqnonlin / fminimax), else fminsearch.

clear; clc;

%% ===== USER SETTINGS =====
targetMR  = 0.75;        % desired ARB MR
objective = 'minmax';    % 'minmax' = minimize the WORST |MR - target| (shrinks peak-to-peak)
                         % 'rms'    = minimize RMS error

% What the optimizer may change (true = free)
vary.MARB_z    = true;   % ARB pivot position along Z
vary.L_ARB     = true;   % ARB arm length
vary.L_Drop    = true;   % droplink length (arm REST ANGLE floats to suit)
vary.BellARB_r = true;   % droplink attach point on bellcrank: radius from bellcrank pivot
vary.BellARB_a = true;   % droplink attach point on bellcrank: angle about bellcrank pivot
vary.MARB_y    = true;   % ARB pivot height
vary.MARB_x    = true;  % ARB pivot position along X

% true -> arm is horizontal at rest, droplink length is derived (L_Drop ignored)
lockArmHorizontal = false;

% Force-transfer (twist efficiency) options.
% eta = fraction of the droplink force that becomes ARB twisting torque
%     = |torque about the bar axis| / (link force * arm length)
% eta = 1 only when the droplink is perpendicular to the arm AND lies in the
% arm's plane (no X-direction component). The rest loads the bar/bearings
% sideways or along the arm instead of twisting it.
twistWeight = 1;     % 0 = off. >0 adds twistWeight*(1 - eta) at every sweep point
minTwistEff = 0.6;     % 0 = off. e.g. 0.95 penalizes any sweep point with eta < 0.95

names = {'MARB_z','L_ARB','L_Drop','BellARB_r','BellARB_a','MARB_y','MARB_x'};
% Bounds (inches; BellARB_a in degrees) -- PLACEHOLDERS, set to your packaging limits
%       PivotPos, L_ARB, L_Drop, BellARB_r, BellARB angle, pivot h,pivot pos
lb = [  4.2, 3, 3.5, 1.5, 110, 3.85, -11.5];
ub = [  6.0, 7.0, 5, 4, 185, 4,  -8.5];

nStarts = 20;            % random restarts (start #1 = current design). Lower for a quick test.
n       = 25;            % sweep points (finer = more accurate MR; slower)

%% ===== GEOMETRY (keep identical to SuspensionAnalysis.m) =====
MSpring = [-9.359109; 15.858556 ;0]; % Spring mount to chassis X Y Z in.
MBell = [-9.3;7;0]; %Bellcrank pivot point
MPush = [-21.72;4.21;0]; % Where pusharm attaches to lower arm

MBellSpr = [-12.403644; 10.67471 ;0]; %Point where bellcrank connects to spring
MBellPush = [-13.763625; 9.208269 ;0]; % bellcrank connects to rod
MBellARB = [ -12.45 ; 8.56 ; 0];

MArmL = [-7.98395 ; 3 ; 6.35]; % Where lower arm connects to chassis
MArmU = [-9.633238 ; 7.305803 ; 6.35]; % upper arm to chassis

MArmL_End = [-20.567206; 3.142493; 0]; % Where knuckle connects to lower arm
MArmU_End = [-21.776179 ; 10.6 ; -0]; %Where knuckle connects to upper arm

minWheelTravel = 1.9685; %minimum travel in inches set by rules
v_knuckle = (MArmU_End-MArmL_End)/norm(MArmU_End-MArmL_End);
MWheel = v_knuckle*(3.12) + MArmL_End;

MARB = [-10.6 ; 3.8 ; 4.32];
L_ARB = 4.53;
MARB_End = [-10.6 ; 5.02 ; .03];

points = [MSpring,MBell,MPush,MBellPush,MBellSpr, ...
    MArmL, MArmU, MArmL_End , MArmU_End,MWheel,MBellARB,MARB,MARB_End];

x_min = 5.94;
x_max = 8.19;

%% ===== SWEEP SETUP =====
% The ARB doesn't affect the linkage, so wheel height is computed once.
x = linspace(x_min, x_max, n);
tp = zeros(3, size(points,2), n);
for i = 1:n
    tp(:,:,i) = DamperX_to_Points(points, x(i));
end
wheelH = squeeze(tp(2,10,:)) - MWheel(2);
wheelH = wheelH(:);

%% ===== DESIGN VARIABLES =====
% p = [Zp, L, Ldl, rB, aB(deg), Yp, Xp]
L_Drop0 = norm(MARB_End - MBellARB);
bv      = MBellARB - MBell;
p0   = [MARB(3), L_ARB, L_Drop0, norm(bv(1:2)), atan2d(bv(2),bv(1)), MARB(2), MARB(1)];
free = [vary.MARB_z, vary.L_ARB, vary.L_Drop && ~lockArmHorizontal, ...
        vary.BellARB_r, vary.BellARB_a, vary.MARB_y, vary.MARB_x];
idx  = find(free);
if isempty(idx), error('Nothing to optimize: enable at least one variable.'); end

outside = idx(p0(idx) < lb(idx) | p0(idx) > ub(idx));
for k = outside
    warning('Current %s = %.4f is outside your bounds [%.3f, %.3f].', names{k}, p0(k), lb(k), ub(k));
end

lbF = lb(idx);  ubF = ub(idx);
q0  = min(max(p0(idx), lbF), ubF);

opt = struct('target',targetMR, 'horizontal',lockArmHorizontal, ...
             'twistWeight',twistWeight, 'minTwistEff',minTwistEff);
res = @(q) arbResiduals(q, idx, p0, points, x, wheelH, opt);

%% ===== OPTIMIZE (multi-start) =====
rng(1);
bestScore = inf;  qBest = q0;
for s = 1:nStarts
    if s == 1
        qs = q0;
    else
        qs = lbF + rand(1,numel(idx)).*(ubF - lbF);
    end
    [q, sc] = runSolver(res, qs, lbF, ubF, objective);
    fprintf('start %2d: %s error = %.4f\n', s, objective, sc);
    if sc < bestScore
        bestScore = sc;  qBest = q;
    end
end

%% ===== RESULTS =====
p = p0;  p(idx) = qBest;
ptsOpt = buildARB(points, p, lockArmHorizontal);
[MR0, Y0, eta0] = sweepARB(points, x, wheelH);
[MR1, Y1, eta1] = sweepARB(ptsOpt, x, wheelH);

restAng1 = atan2d(ptsOpt(2,13)-ptsOpt(2,12), -(ptsOpt(3,13)-ptsOpt(3,12)));
L_Drop1  = norm(ptsOpt(:,13) - ptsOpt(:,11));

fprintf('\n=== RESULT (target MR = %.3f, objective = %s) ===\n', targetMR, objective);
for k = idx
    span = ub(k) - lb(k);
    flag = '';
    if abs(p(k)-lb(k)) < 0.01*span || abs(p(k)-ub(k)) < 0.01*span
        flag = '   <-- AT BOUND (limiting the result)';
    end
    fprintf('  %-10s %9.4f   (was %9.4f)%s\n', names{k}, p(k), p0(k), flag);
end
fprintf('Droplink length : %.4f in  (was %.4f)\n', L_Drop1, L_Drop0);
fprintf('Arm rest angle  : %.2f deg above horizontal (was 0)\n', restAng1);
fprintf('MR over travel  : %.3f to %.3f   (was %.3f to %.3f)\n', min(MR1), max(MR1), min(MR0), max(MR0));
fprintf('Peak-to-peak MR : %.3f   (was %.3f)\n', max(MR1)-min(MR1), max(MR0)-min(MR0));
fprintf('Std dev of MR   : %.4f   (was %.4f)\n', std(MR1), std(MR0));
fprintf('Mean MR         : %.3f   (was %.3f)\n', mean(MR1), mean(MR0));
fprintf('Twist efficiency: min %.3f, mean %.3f   (was min %.3f, mean %.3f)\n', ...
        min(eta1), mean(eta1), min(eta0), mean(eta0));
if bestScore >= 50
    warning('Best result still hit the infeasible penalty - widen bounds or change the target.');
end

fprintf('\nPaste into SuspensionAnalysis.m (replace the existing lines):\n');
fprintf('MBellARB = [%.6f ; %.6f ; %.6f];\n', ptsOpt(:,11));
fprintf('MARB     = [%.6f ; %.6f ; %.6f];\n', ptsOpt(:,12));
fprintf('MARB_End = [%.6f ; %.6f ; %.6f];\n', ptsOpt(:,13));

figure(10); clf;
subplot(1,3,1);
plot(wheelH, MR0, '--', wheelH, MR1, '-', 'LineWidth', 1.5); hold on;
yline(targetMR, 'k:');
grid on; legend('Original','Optimized','Target','Location','best');
xlabel('Wheel height from rest'); ylabel('ARB MR'); title('ARB Instantaneous MR');
subplot(1,3,2);
plot(wheelH, Y0 - Y0(1), '--', wheelH, Y1 - Y1(1), '-', 'LineWidth', 1.5);
grid on; legend('Original','Optimized','Location','best');
xlabel('Wheel height from rest'); ylabel('ARB end height (rel. to first point)');
title('ARB end travel');
subplot(1,3,3);
plot(wheelH, eta0, '--', wheelH, eta1, '-', 'LineWidth', 1.5);
grid on; ylim([0 1.05]); legend('Original','Optimized','Location','best');
xlabel('Wheel height from rest'); ylabel('Twist efficiency \eta');
title('Droplink force into ARB twist');

%% ===== LOCAL FUNCTIONS =====
function r = arbResiduals(q, idx, p0, pts0, x, wheelH, opt)
    nTerms = 1 + (opt.twistWeight > 0) + (opt.minTwistEff > 0);
    p = p0;  p(idx) = q;
    try
        pts = buildARB(pts0, p, opt.horizontal);
        [MR, ~, eta] = sweepARB(pts, x, wheelH);
        r = MR(:) - opt.target;
        if opt.twistWeight > 0
            r = [r; opt.twistWeight*(1 - eta(:))];
        end
        if opt.minTwistEff > 0
            r = [r; 10*max(0, opt.minTwistEff - eta(:))];
        end
    catch
        r = 50*ones(numel(x)*nTerms,1);    % infeasible geometry -> big penalty
    end
end

function [MR, Y, eta] = sweepARB(pts, x, wheelH)
    Y   = zeros(numel(x),1);
    eta = zeros(numel(x),1);
    for i = 1:numel(x)
        np   = DamperX_to_Points(pts, x(i));
        Y(i) = np(2,13);
        a = np(:,13) - np(:,12);               % ARB arm (pivot -> end)
        d = np(:,11) - np(:,13);               % droplink (ARB end -> bellcrank)
        d = d / norm(d);
        eta(i) = abs(a(2)*d(3) - a(3)*d(2)) / norm(a);   % torque about X per (force * arm length)
    end
    MR = gradient(Y, wheelH);          % same definition as SuspensionAnalysis.m
end

function pts = buildARB(pts0, p, horizontal)
    % Rebuilds droplink bellcrank point (col 11), ARB pivot (col 12), ARB end (col 13).
    Zp = p(1); L = p(2); Ldl = p(3); rB = p(4); aB = p(5); Yp = p(6); Xp = p(7);
    pts = pts0;
    pts(:,11) = pts0(:,2) + rB*[cosd(aB); sind(aB); 0];   % rides on the bellcrank
    pts(:,12) = [Xp; Yp; Zp];
    if horizontal
        pts(:,13) = [Xp; Yp; Zp - L];
        return
    end
    % Droplink length fixed -> intersect (arm circle about pivot) with
    % (droplink circle about the bellcrank ARB point) in the Y-Z plane.
    dx = Xp - pts(1,11);
    r2 = Ldl^2 - dx^2;
    if r2 <= 0, error('Droplink shorter than its X offset.'); end
    [yo,zo] = circcirc(Yp, Zp, L, pts(2,11), pts(3,11), sqrt(r2));
    if isempty(yo) || any(isnan(yo)), error('ARB arm / droplink circles do not intersect.'); end
    [~,k] = min(hypot(yo - Yp, zo - (Zp - L)));           % branch closest to horizontal
    pts(:,13) = [Xp; yo(k); zo(k)];
end

function [q, score] = runSolver(res, q0, lb, ub, objective)
    if license('test','Optimization_Toolbox')
        o1 = optimoptions('lsqnonlin','Display','off', ...
            'FiniteDifferenceStepSize',1e-4,'FunctionTolerance',1e-10,'StepTolerance',1e-8);
        q = lsqnonlin(res, q0, lb, ub, o1);
        if strcmp(objective,'minmax')
            o2 = optimoptions('fminimax','Display','off', ...
                'MaxIterations',300,'MaxFunctionEvaluations',6000);
            try
                qb = fminimax(@(v) twoSided(res(v)), q, [],[],[],[], lb, ub, [], o2);
                if scoreOf(res(qb),objective) < scoreOf(res(q),objective), q = qb; end
            catch
            end
        end
    else
        clip = @(v) min(max(v,lb),ub);
        sc   = @(v) scoreOf(res(clip(v)),objective) + 1e3*sum((v-clip(v)).^2);
        q = fminsearch(sc, q0, optimset('TolX',1e-8,'TolFun',1e-10,'MaxFunEvals',6000));
        q = clip(q);
    end
    score = scoreOf(res(q), objective);
end

function F = twoSided(r)
    F = [r; -r];                       % fminimax minimizes max(F) = max|r|
end

function s = scoreOf(r, objective)
    if strcmp(objective,'minmax'), s = max(abs(r)); else, s = sqrt(mean(r.^2)); end
end