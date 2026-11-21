function newPoints = DamperX_to_Points(points, x)
newPoints = points;

% Spring/bellcrank point: circle around spring mount (length x) and bellcrank pivot
L_bell_spring = norm(points(:,5) - points(:,2));
[xo,yo] = circcirc(points(1,1),points(2,1),x, ...
    points(1,2),points(2,2),L_bell_spring);
[newPoints(1,5), newPoints(2,5)] = pickNearest(xo,yo,points(1,5),points(2,5));

% Bellcrank rotation (CCW positive)
v_old = points(1:2,5) - points(1:2,2);
v_new = newPoints(1:2,5) - points(1:2,2);
dAngle = atan2d(v_new(2),v_new(1)) - atan2d(v_old(2),v_old(1));
newPoints(:,4) = points(:,2) + rotz(dAngle)*(points(:,4) - points(:,2));
newPoints(:,11) = points(:,2) + rotz(dAngle)*(points(:,11) - points(:,2));

% Pushrod / lower arm end
L_ArmL    = norm(points(1:2,3) - points(1:2,6));
L_pushrod = norm(points(:,3) - points(:,4));
[xo,yo] = circcirc(newPoints(1,4),newPoints(2,4),L_pushrod, ...
    points(1,6),points(2,6),L_ArmL);
[newPoints(1,3), newPoints(2,3)] = pickNearest(xo,yo,points(1,3),points(2,3));

% Lower arm end
v_old = points(1:2,3) - points(1:2,6);
v_new = newPoints(1:2,3) - points(1:2,6);
dAngle = atan2d(v_new(2),v_new(1)) - atan2d(v_old(2),v_old(1));
dAngle = atan2d(sind(dAngle), cosd(dAngle));
newPoints(:,8) = points(:,6) + rotz(dAngle)*(points(:,8) - points(:,6));


% Upper arm end
L_knuckle = norm(points(:,9) - points(:,8));
L_ArmU    = norm(points(:,9) - points(:,7));
[xo,yo] = circcirc(newPoints(1,8),newPoints(2,8),L_knuckle, ...
    points(1,7),points(2,7),L_ArmU);
[newPoints(1,9), newPoints(2,9)] = pickNearest(xo,yo,points(1,9),points(2,9));

% wheel location
v_knuckle = (newPoints(:,9)-newPoints(:,8) ) / norm(newPoints(:,9)-newPoints(:,8) );
newPoints(:,10) = v_knuckle*(3.26) + newPoints(:,8);

% Find ARB rotation
L_DropLink = norm( points(:,13)-points(:,11) );
v_ARB = points(:,13)-points(:,12);
v_Drop = newPoints(:,11)-points(:,12);

a = v_ARB(2)*v_Drop(2) + v_ARB(3)*v_Drop(3);
b = v_ARB(2)*v_Drop(3) - v_ARB(3)*v_Drop(2);
c = (norm(v_ARB)^2 + norm(v_Drop)^2 - L_DropLink^2) /2 - v_ARB(1)*v_Drop(1);

%check if solution is reachable
R = hypot(a,b);
if abs(c) > R
    error('Target distance is unreachable by rotating around the x-axis')
end

%calculate the two exact angles
phi = atan2(b,a);
alpha = acos(c/R);
theta_1 = wrapToPi(phi + alpha);
theta_2 = wrapToPi(phi - alpha);
Rx = @(th) [1, 0, 0; 0, cos(th), -sin(th); 0, sin(th), cos(th)];
E1 = points(:,12) + (Rx(theta_1) * v_ARB);
E2 = points(:,12) + (Rx(theta_2) * v_ARB);
newPoints(:,13) = E2;

end

function [px,py] = pickNearest(xo,yo,x0,y0)
if isempty(xo) || any(isnan(xo))
    error('Linkage cannot reach this damper length (circles do not intersect).');
end
[~,k] = min(hypot(xo - x0, yo - y0));
px = xo(k); py = yo(k);
end