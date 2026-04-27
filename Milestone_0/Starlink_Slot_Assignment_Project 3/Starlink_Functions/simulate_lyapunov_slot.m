function sol = simulate_lyapunov_slot(sat, slot, pmp, cfg)
%SIMULATE_LYAPUNOV_SLOT Lyapunov SMA tracking for assigned slot.
%
% Reference:
%   a_ref(t) = aParking during the RAAN-drift coast phase,
%              aFinal after the PMP switching time.
%
% Lyapunov function:
%   V = 1/2 (a - a_ref)^2
%
% Since a_dot = 2u/n(a), choosing
%   u = -sat(k_a (a - a_ref), uMax)
% gives V_dot <= 0 whenever a_ref is constant.
%
% This is not intended to be a global optimizer. It is a feedback realization
% of the same Starlink-style RAAN-drift-plus-raise strategy.

    x0 = [sat.a0; ...
          wrapToPiLocal(deg2rad(sat.raan0_deg - slot.raan0_deg)); ...
          wrapToPiLocal(deg2rad(sat.M0_deg    - slot.M0_deg))];

    tEnd = pmp.tf*(1 + cfg.lyap.extraTimeFraction) + cfg.lyap.maxExtraDays*86400;
    opts = odeset('RelTol',1e-9, 'AbsTol',1e-11, ...
                  'Events', @(t,x) lyap_terminal_event(t,x,slot,pmp,cfg));

    [t, x] = ode45(@(t,x) lyap_dyn(t,x,slot,pmp,cfg), [0 tEnd], x0, opts);

    u = zeros(numel(t),1);
    for k = 1:numel(t)
        u(k) = lyap_control(t(k), x(k,:)', pmp, cfg);
    end

    final = x(end,:)';
    sol.method = "Lyapunov";
    sol.satId = sat.id;
    sol.slotId = slot.id;
    sol.tf = t(end);
    sol.t = t;
    sol.a = x(:,1);
    sol.relRAAN = arrayfun(@wrapToPiLocal, x(:,2));
    sol.relM = arrayfun(@wrapToPiLocal, x(:,3));
    sol.u = u;
    sol.deltaV_mps = 1000 * trapz(t, abs(u));
    sol.finalAError_km = final(1) - slot.aF;
    sol.finalRAANError = wrapToPiLocal(final(2));
    sol.finalMError = wrapToPiLocal(final(3));
    sol.arrivalEpoch = cfg.launchEpoch + seconds(sol.tf);

    if abs(sol.finalAError_km) < cfg.lyap.tolA_km && ...
       abs(rad2deg(sol.finalRAANError)) < cfg.lyap.tolRAAN_deg
        sol.status = "Converged";
        sol.success = true;
    else
        sol.status = "Ended at max time or tolerance not reached";
        sol.success = false;
    end
end

function xdot = lyap_dyn(t, x, slot, pmp, cfg)
    u = lyap_control(t, x, pmp, cfg);
    xdot = starlink_elem_dyn(t, x, u, slot, cfg);
end

function u = lyap_control(t, x, pmp, cfg)
    if t < pmp.tCoast
        aRef = pmp.a(1);
    else
        aRef = pmp.a(end);
    end

    eA = x(1) - aRef;

    % Saturated Lyapunov feedback.
    uUnsat = -cfg.lyap.kA * eA;
    u = max(min(uUnsat, cfg.uMax_km_s2), -cfg.uMax_km_s2);

    % Do not lower below the drift shell in this scenario.
    if x(1) <= cfg.aMin && u < 0
        u = 0;
    end

    % During the planned coast, hold the lower SMA exactly. This preserves
    % the intended J2 differential precession strategy.
    if t < pmp.tCoast
        u = 0;
    end
end

function [value, isterminal, direction] = lyap_terminal_event(t, x, slot, pmp, cfg)
    aErr = abs(x(1) - slot.aF) - cfg.lyap.tolA_km;
    raanErr = abs(wrapToPiLocal(x(2))) - deg2rad(cfg.lyap.tolRAAN_deg);

    % Do not allow termination before the planned PMP time scale.
    timeGate = t - 0.99*pmp.tf;

    value = max([aErr, raanErr, -timeGate]);
    isterminal = 1;
    direction = -1;
end
