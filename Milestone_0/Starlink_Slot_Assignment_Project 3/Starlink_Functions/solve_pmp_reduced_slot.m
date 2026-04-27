function sol = solve_pmp_reduced_slot(sat, slot, cfg)
%SOLVE_PMP_REDUCED_SLOT Reduced-order path-constrained PMP solution.
%
% This implements the PMP structure for the simplified Starlink deployment
% model where only SMA is actuated and RAAN is obtained from J2:
%
%   minimize t_f
%   subject to:
%       a_dot = 2 u_T / n(a), |u_T| <= uMax
%       Omega_dot = J2 secular nodal rate
%       a >= aParking
%
% For lower-RAAN slots, the minimum-time strategy is:
%   1) coast on the lowest permitted SMA, a = aParking, so RAAN regresses
%      as fast as possible,
%   2) then apply maximum positive tangential thrust until a = aFinal.
%
% This is a bang/path-arc structure:
%   lower-SMA boundary arc -> max-thrust bang arc.

    est = estimate_pmp_pair(sat, slot, cfg);

    tCoast = est.tCoast;
    tRaise = est.tRaise;
    tf = est.tf;

    x0 = [sat.a0; ...
          wrapToPiLocal(deg2rad(sat.raan0_deg - slot.raan0_deg)); ...
          wrapToPiLocal(deg2rad(sat.M0_deg    - slot.M0_deg))];

    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

    if tCoast > 0
        t1 = linspace(0, tCoast, max(3, ceil(tCoast/900)));
        x1 = zeros(numel(t1),3);
        for k = 1:numel(t1)
            tau = t1(k);
            x1(k,:) = [sat.a0, ...
                       x0(2) + (j2_raan_rate(sat.a0,cfg.inc_rad,cfg.Earth) - j2_raan_rate(slot.aF,cfg.inc_rad,cfg.Earth))*tau, ...
                       x0(3) + (circ_mean_rate(sat.a0,cfg.Earth) - circ_mean_rate(slot.aF,cfg.Earth))*tau];
        end
        u1 = zeros(numel(t1),1);
        xStartRaise = x1(end,:)';
    else
        t1 = 0;
        x1 = x0';
        u1 = 0;
        xStartRaise = x0;
    end

    if tRaise > 0
        [t2local, x2] = ode45(@(t,x) starlink_elem_dyn(t, x, cfg.uMax_km_s2, slot, cfg), ...
                              [0 tRaise], xStartRaise, opts);
        t2 = tCoast + t2local;
        u2 = cfg.uMax_km_s2*ones(numel(t2),1);

        % Avoid duplicate sample at switching time.
        if numel(t1) > 1
            tHist = [t1(:); t2(2:end)];
            xHist = [x1; x2(2:end,:)];
            uHist = [u1(:); u2(2:end)];
        else
            tHist = t2(:);
            xHist = x2;
            uHist = u2(:);
        end
    else
        tHist = t1(:);
        xHist = x1;
        uHist = u1(:);
    end

    final = xHist(end,:)';
    sol.method = "Reduced PMP";
    sol.satId = sat.id;
    sol.slotId = slot.id;
    sol.tCoast = tCoast;
    sol.tRaise = tRaise;
    sol.tf = tf;
    sol.t = tHist;
    sol.a = xHist(:,1);
    sol.relRAAN = arrayfun(@wrapToPiLocal, xHist(:,2));
    sol.relM = arrayfun(@wrapToPiLocal, xHist(:,3));
    sol.u = uHist;
    sol.deltaV_mps = 1000 * trapz(tHist, abs(uHist));
    sol.finalAError_km = final(1) - slot.aF;
    sol.finalRAANError = wrapToPiLocal(final(2));
    sol.finalMError = wrapToPiLocal(final(3));
    sol.arrivalEpoch = cfg.launchEpoch + seconds(tf);
    sol.status = est.status;

    if abs(sol.finalAError_km) < cfg.aTol_km && abs(rad2deg(sol.finalRAANError)) < cfg.raanTol_deg
        sol.success = true;
    else
        sol.success = false;
    end
end
