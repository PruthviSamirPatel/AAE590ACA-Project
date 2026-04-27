function aJ2 = accelJ2_cart(r, mu, J2, Re)
    x = r(1);
    y = r(2);
    z = r(3);

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