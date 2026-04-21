function [phi, gradPhi] = radiusPenaltyCost(r, rSafe, rOn, epsBar)

rmag = norm(r);

% zero penalty outside activation region
if rmag >= rOn
    phi = 0;
    gradPhi = zeros(3,1);
    return
end

% avoid divide by zero
if rmag < epsBar
    phi = 1;
    gradPhi = zeros(3,1);
    return
end

Delta = rOn - rSafe;
eta = (rOn - rmag)/Delta;

% clip to [0,1]
eta = min(max(eta,0),1);

% smooth bounded penalty
phi = eta^4;

% derivative dphi/drmag
dphi_drmag = -4*eta^3 / Delta;

% chain rule
gradPhi = dphi_drmag * (r/rmag);

end