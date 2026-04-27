function aJ2 = accelJ2(r_eci, mu, J2, Re)
% J2 acceleration in inertial coordinates
% Inputs must all be in consistent nondimensional units

    x = r_eci(1);
    y = r_eci(2);
    z = r_eci(3);

    r2 = x^2 + y^2 + z^2;
    r1 = sqrt(r2);
    z2 = z^2;

    factor = -(3/2) * J2 * mu * Re^2 / r1^5;

    common = 1 - 5*z2/r2;

    ax = factor * x * common;
    ay = factor * y * common;
    az = factor * z * (3 - 5*z2/r2);

    aJ2 = [ax; ay; az];
end