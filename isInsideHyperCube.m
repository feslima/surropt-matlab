% isInside = 1 (inside specified tolerance)
% isInside = 2 (at limit of hypercube)
% isInside = 0 (outside hypercube)

function isInside = isInsideHyperCube(point,lb,ub)
% Determine if a point is inside a hypercube following a specified
% tolerance

% tolerance set to % of design range
tol = 1e-6; % 0.00001

% Normalize the bounds and point
% absolUB = norm(ub-lb);
absolUB = abs(ub);
ubNorm = ub ./ absolUB;
lbNorm = lb ./ absolUB;
pointNorm = point ./ absolUB;

% Check if point is outside the bounds
if any(bsxfun(@gt,pointNorm,ubNorm)) || any(bsxfun(@lt,pointNorm, lbNorm))
    isInside = 0;
    return
end

% Check if point is inside or at the bounds
if any(abs(pointNorm - ubNorm) <= tol) || any(abs(pointNorm - lbNorm) <= tol) % point is at hypercube limit
    isInside = 2;
else
    isInside = 1; % inside hypercube
end