%% 
%% ARB_Geometry_Optimizer.m
% Tunes the ARB geometry so that
%   1) ARB motion ratio is as constant as possible over the damper sweep
%   2) as much droplink force as possible goes into ARB torsion
%   3) the bar sees as little bending (cantilever) as possible
%
% Design variables (7):
%   MARBx, MARBy, MARBz   ARB axis point / lever-arm root   [in]
%   L_ARB                 lever arm length                    [in]
%   L_Drop                droplink length                     [in]
%   BellARB_r, BellARB_ang  droplink pickup on the bellcrank, polar
%                           about the bellcrank pivot (XY plane) [in, deg]
%
% Needs: DamperX_to_Points.m (with the branch patch), Mapping Toolbox (circcirc,
%        wrapToPi), Optimization Toolbox (fmincon)
%
% Modelling assumptions:
%   - ARB axis runs along X. The lever arm lies in the plane x = MARBx.
%   - At rest the arm end is found from L_ARB and L_Drop (circle intersection
%     in the YZ plane), so L_Drop and L_ARB can be set independently.
clear; clc; close all;

%% ---------------- Linkage (copy of SuspensionAnalysis.m, keep in sync) ---------MSpring = [-9.359109; 15.858556 ;0]; % Spring mount to chassis X Y Z in.
MSpring = [-8.74; 15.858556 ;0]; % Spring mount to chassis X Y Z in.
MBell = [-9.3;7.29;0]; %Bellcrank pivot point
MPush = [-18.67;3.84;0]; % Where pusharm attaches to lower arm

MBellSpr = [-12.82; 10.57 ;0]; %Point where bellcrank connects to spring
MBellPush = [-13.99; 8.96 ;0]; % bellcrank connects to rod
MBellARB = [ -13.29 ; 7.82 ; 0];

MArmL = [-7.96;3.01;5.6]; % Where lower arm connects to chassis
MArmU = [-9.61;7.3;5.6]; % upper arm to chassis

MArmL_End = [-20.73;4.04;0]; % Where knuckle connects to lower arm
MArmU_End = [-19.73;11.58;0]; %Where knuckle connects to upper arm

minWheelTravel = 1.9685; %minimum travel in inches set by rules
v_knuckle = (MArmU_End-MArmL_End)/norm(MArmU_End-MArmL_End);
MWheel = v_knuckle*(3.26) + MArmL_End;

MARB = [-11.85; 2.23 ; 2.9];
L_ARB = 2.8;
MARB_End = [-11.85 ; 2.23 ; MARB(3)-L_ARB];

% Starting ARB design (the one at the bottom of SuspensionAnalysis.m)
MBellARB0 = MBellARB;
MARB0     = MARB;
MARB_End0 = MARB_End;

%% ---------------- Settings -----------------------------------------------------
S.x_min = 5.94;  S.x_max = 8.19;  S.n = 15;   % damper sweep
S.mrMode   = 'Y';   % 'Y'     = ARB end Y travel / wheel travel (same as SuspensionAnalysis)
                    % 'angle' = ARB arm rotation [deg] / wheel travel [in]
S.targetMR = [];    % []     = just make MR constant (level chosen by optimiser, kept in MRband)
                    % number = hold that MR (units follow mrMode). Sign matters: the wheel
                    %          height falls as the damper extends, so MR is likely negative -
                    %          check the baseline MR printed below before choosing a target.
S.MRband   = [0.5 1.5];   % allowed mean|MR| as a fraction of the start design's (targetMR = [] only)
S.xBearing = -7.15;  % <<< PLACEHOLDER. X of the chassis bearing/support nearest the lever arm.
                    %     Overhang = MARBx - xBearing. Set this to your real value.
S.w = struct('mr',1, 'eta',0, 'bend',0);   % cost weights
S.nStarts = 30;     % random multi-start count (first start = current design)

% Bounds are the start design +/- span. Edit freely.
%               MARBx MARBy MARBz L_ARB L_Drop BellARB_r BellARB_ang(deg)
span       = [   .5     0.4    0.2    .5    2     0      0 ];
varNames   = {'MARBx','MARBy','MARBz','L_ARB','L_Drop','BellARB_r','BellARB_ang'};

%% ---------------- Setup --------------------------------------------------------
B.base  = [MSpring,MBell,MPush,MBellPush,MBellSpr,MArmL,MArmU,MArmL_End,MArmU_End,MWheel];
B.MBell = MBell;
B.End0  = MARB_End0;
B.MWheel = MWheel;
B.x = linspace(S.x_min, S.x_max, S.n);

dBell = MBellARB0(1:2) - MBell(1:2);
v0 = [MARB0(1), MARB0(2), MARB0(3), norm(MARB_End0-MARB0), ...
      norm(MARB_End0-MBellARB0), norm(dBell), atan2d(dBell(2),dBell(1))];
lb = v0 - span;  ub = v0 + span;
lb(4) = max(lb(4), 0.5);  lb(5) = max(lb(5), 1.0);  lb(6) = max(lb(6), 0.5);

% Baseline
S.MR0 = NaN;
S.debug = true;                       % baseline: show the real error if something fails
[~, i0] = evalDesign(v0, B, S);
S.debug = false;                      % optimisation: failures just get a big penalty
S.MR0 = mean(abs(i0.MR));
fprintf('Rebuilt rest ARB end vs. your MARB_End: %.2e in off (should be ~0)\n', ...
        norm(i0.pts(:,13) - MARB_End0));

%% ---------------- Optimise (normalised variables, multi-start) -----------------
toV  = @(z) lb + z.*(ub-lb);
fun  = @(z) evalDesign(toV(z), B, S);
opts = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
        'FiniteDifferenceStepSize',1e-4,'StepTolerance',1e-9,'MaxFunctionEvaluations',3000);
rng(1);
z0s = [(v0-lb)./(ub-lb); rand(S.nStarts-1, numel(v0))];
bestJ = inf;  bestZ = [];
for k = 1:S.nStarts
    z0 = z0s(k,:);
    if fun(z0) >= 1e3, continue; end          % skip infeasible starts
    try
        [z,J] = fmincon(fun, z0, [],[],[],[], zeros(1,numel(v0)), ones(1,numel(v0)), [], opts);
    catch
        continue
    end
    if J < bestJ, bestJ = J; bestZ = z; end
    fprintf('start %2d: J = %.4f   (best %.4f)\n', k, J, bestJ);
end
if isempty(bestZ), warning('No feasible optimum found - showing start design.'); vBest = v0;
else,              vBest = toV(bestZ); end
[~, ib] = evalDesign(vBest, B, S);

%% ---------------- Report -------------------------------------------------------
disp(table(varNames', v0', vBest', lb', ub', ...
    'VariableNames', {'Variable','Start','Optimised','LowerBound','UpperBound'}));
report = @(nm,i) fprintf('%-10s MR mean %.3f | MR spread (max-min) %.3f | torsion eff min/mean %.3f/%.3f | bend/torsion max %.3f\n', ...
    nm, mean(i.MR), max(i.MR)-min(i.MR), min(i.eta), mean(i.eta), max(i.bend));
fprintf('\n'); report('Start', i0); report('Optimised', ib);

pb = ib.pts;
fprintf('\n%% Paste into SuspensionAnalysis.m:\n');
fprintf('MBellARB = [%.6f ; %.6f ; %.6f];\n', pb(:,11));
fprintf('MARB     = [%.6f ; %.6f ; %.6f];\n', pb(:,12));
fprintf('MARB_End = [%.6f ; %.6f ; %.6f];\n', pb(:,13));

figure(10); clf;
subplot(1,3,1); plot(i0.wheel,i0.MR,'--',ib.wheel,ib.MR,'-','LineWidth',1.5);
legend('Start','Optimised','Location','best'); title('Instantaneous ARB MR');
xlabel('Wheel travel from rest [in]'); ylabel(['MR (' S.mrMode ')']);
subplot(1,3,2); plot(i0.wheel,i0.eta,'--',ib.wheel,ib.eta,'-','LineWidth',1.5);
title('Droplink force into torsion (1 = all)'); xlabel('Wheel travel from rest [in]'); ylim([0 1.05]);
subplot(1,3,3); plot(i0.wheel,i0.bend,'--',ib.wheel,ib.bend,'-','LineWidth',1.5);
title('Bar bending / torsion moment'); xlabel('Wheel travel from rest [in]');

%% ---------------- Local functions ----------------------------------------------
function [J, info] = evalDesign(v, B, S)
% v = [MARBx MARBy MARBz L_ARB L_Drop BellARB_r BellARB_ang]
info = [];
try
    pts = buildPoints(v, B);
    n = numel(B.x);
    P = zeros(3, size(pts,2), n);
    for i = 1:n
        P(:,:,i) = DamperX_to_Points(pts, B.x(i));
    end
catch ME
    if S.debug, rethrow(ME); end
    J = 1e3;  return
end
wheel = squeeze(P(2,10,:)) - B.MWheel(2);               % n x 1
dW = diff(wheel);
if ~(all(dW > 0) || all(dW < 0)) || ~isreal(P)   % must be monotonic in either direction
    if S.debug, error('Wheel height is not monotonic over the damper sweep (or P is complex).'); end
    J = 1e3;  return
end
p11 = squeeze(P(:,11,:));                                % droplink bellcrank end, 3 x n
p13 = squeeze(P(:,13,:));                                % droplink ARB end,       3 x n
p12 = pts(:,12);                                         % ARB axis point
r   = p13 - p12;

% --- motion ratio
switch S.mrMode
    case 'Y',     q = p13(2,:)';
    case 'angle', q = rad2deg(unwrap(atan2(r(3,:), r(2,:))))';
    otherwise,    error('mrMode must be ''Y'' or ''angle''');
end
MR = gradient(q, wheel);

% --- droplink force (unit) acting on the ARB end, moments about the bearing
F  = p11 - p13;  F = F ./ vecnorm(F);
rb = p13 - [S.xBearing; p12(2); p12(3)];                 % bearing point -> arm end
M  = cross(rb, F, 1);                                    % moment of the force about the bearing
Mt = abs(M(1,:));                                        % torsion (about X)
Mb = vecnorm(M(2:3,:), 2, 1);                            % bending (about Y and Z)
eta  = Mt ./ vecnorm(r(2:3,:), 2, 1);                    % fraction of force that twists the bar
bend = Mb ./ max(Mt, 1e-3);                              % bending moment per unit torsion

% --- cost
mAbs = max(mean(abs(MR)), 1e-6);
if isempty(S.targetMR)
    mrCost = std(MR) / mAbs;
    band = 0;
    if ~isnan(S.MR0)
        lo = S.MRband(1)*S.MR0;  hi = S.MRband(2)*S.MR0;
        band = max(0, lo-mAbs)/lo + max(0, mAbs-hi)/hi;
    end
else
    mrCost = sqrt(mean((MR - S.targetMR).^2)) / abs(S.targetMR);
    band = 0;
end
J = S.w.mr*mrCost + 10*band + S.w.eta*((1-mean(eta)) + (1-min(eta))) + S.w.bend*max(bend);

info = struct('pts',pts,'wheel',wheel,'MR',MR,'eta',eta','bend',bend','J',J);
end

function pts = buildPoints(v, B)
MARB = v(1:3)';  L_ARB = v(4);  L_Drop = v(5);  r = v(6);  ang = v(7);
MBellARB = B.MBell + [r*cosd(ang); r*sind(ang); 0];
dx = MBellARB(1) - MARB(1);
if L_Drop <= abs(dx), error('droplink shorter than its X offset'); end
Lyz = sqrt(L_Drop^2 - dx^2);                              % droplink length seen in the YZ plane
[yo, zo] = circcirc(MARB(2), MARB(3), L_ARB, MBellARB(2), MBellARB(3), Lyz);
if isempty(yo) || any(isnan(yo)), error('arm and droplink circles do not meet'); end
[~, k] = min(hypot(yo - B.End0(2), zo - B.End0(3)));      % same side as the current design
MARB_End = [MARB(1); yo(k); zo(k)];
pts = [B.base, MBellARB, MARB, MARB_End];
end

