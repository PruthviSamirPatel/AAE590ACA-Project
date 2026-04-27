function xdot = starlink_reduced_elem_dyn(~, x, u, slot, cfg)
%STARLINK_REDUCED_ELEM_DYN Reduced circular low-thrust + J2 dynamics.
%
% State:
%   x(1) = a [km]
%   x(2) = relative RAAN = Omega_sat - Omega_slot(t) [rad]
%   x(3) = relative M    = M_sat     - M_slot(t)     [rad]
%
% Control:
%   u = tangential acceleration [km/s^2]

    a = x(1);

    if a <= cfg.aMin && u < 0
        u = 0;
    end
    if a >= cfg.aMax && u > 0
        u = 0;
    end

    n = circ_mean_rate(a, cfg.Earth);
    nF = circ_mean_rate(slot.aF, cfg.Earth);

    OmDot = j2_raan_rate(a, cfg.inc_rad, cfg.Earth);
    OmDotF = j2_raan_rate(slot.aF, cfg.inc_rad, cfg.Earth);

    adot = 2*u/n;
    relOmdot = OmDot - OmDotF;
    relMdot = n - nF;

    xdot = [adot; relOmdot; relMdot];
end
