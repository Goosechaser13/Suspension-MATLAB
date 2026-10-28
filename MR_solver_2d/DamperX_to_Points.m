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

% Pushrod / lower arm end
L_ArmL    = norm(points(:,3) - points(:,6));
L_pushrod = norm(points(:,3) - points(:,4));
[xo,yo] = circcirc(newPoints(1,4),newPoints(2,4),L_pushrod, ...
    points(1,6),points(2,6),L_ArmL);
[newPoints(1,3), newPoints(2,3)] = pickNearest(xo,yo,points(1,3),points(2,3));

% Lower arm end follows the pushrod mount (same point for now)
newPoints(:,8) = newPoints(:,3);

% Upper arm end
L_knuckle = norm(points(:,9) - points(:,8));
L_ArmU    = norm(points(:,9) - points(:,7));
[xo,yo] = circcirc(newPoints(1,8),newPoints(2,8),L_knuckle, ...
    points(1,7),points(2,7),L_ArmU);
[newPoints(1,9), newPoints(2,9)] = pickNearest(xo,yo,points(1,9),points(2,9));

% wheel location
v_knuckle = (newPoints(:,9)-newPoints(:,8) ) / norm(newPoints(:,9)-newPoints(:,8) );
newPoints(:,10) = v_knuckle*(0.97) + newPoints(:,8);
end

function [px,py] = pickNearest(xo,yo,x0,y0)
if isempty(xo) || any(isnan(xo))
    error('Linkage cannot reach this damper length (circles do not intersect).');
end
[~,k] = min(hypot(xo - x0, yo - y0));
px = xo(k); py = yo(k);
end