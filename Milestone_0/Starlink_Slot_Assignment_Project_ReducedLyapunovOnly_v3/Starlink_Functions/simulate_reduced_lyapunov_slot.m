function sol = simulate_reduced_lyapunov_slot(sat, slot, cfg)
%SIMULATE_REDUCED_LYAPUNOV_SLOT Propagate one reduced Lyapunov transfer.
%
% This function deliberately avoids the previous "instant jump" behavior by
% using the Project-3 low-thrust acceleration and by plotting on a fixed
% time grid. If uMax is set to 1 m/s^2, the raise really only lasts about
% 112 seconds, so it will appear vertical on a days-scale plot.

    plan = estimate_reduced_lyapunov_pair(sat, slot, cfg);
    plan.aF = slot.aF;

    x0 = [sat.a0; ...
          wrapToPiLocal(deg2rad(sat.raan0_deg - slot.raan0_deg)); ...
          wrapToPiLocal(deg2rad(sat.M0_deg    - slot.M0_deg))];

    % End time based on the planned coast and low-thrust raise duration.
    tEnd = plan.tCoast + ...
           (1 + cfg.lyap.extraRaiseTimeFraction)*plan.tRaise + ...
           cfg.lyap.maxExtraDays*86400;

    opts = odeset('RelTol',1e-10, 'AbsTol',1e-12, ...
                  'MaxStep', cfg.lyap.dtPlot, ...
                  'Events', @(t,x) reduced_lyap_terminal_event(t,x,plan,cfg));

    [t, x] = ode45(@(t,x) reduced_lyap_ode(t,x,slot,plan,cfg), [0 tEnd], x0, opts);

    % Uniform output grid for clean project plots. Include the coast-to-
    % control switch and final point exactly so the mode transition is clear.
    if numel(t) > 1
        tPlot = unique([0; plan.tCoast; (0:cfg.lyap.dtPlot:t(end))'; t(end)]);
        tPlot = tPlot(tPlot >= 0 & tPlot <= t(end));
        xPlot = interp1(t, x, tPlot, 'pchip');
    else
        tPlot = t;
        xPlot = x;
    end

    u = zeros(numel(tPlot),1);
    mode = zeros(numel(tPlot),1); % 0 coast/hold, 1 Lyapunov thrust
    for k = 1:numel(tPlot)
        u(k) = reduced_lyapunov_control(tPlot(k), xPlot(k,:)', plan, cfg);
        if abs(u(k)) > 1e-14
            mode(k) = 1;
        end
    end

    final = xPlot(end,:)';

    rawDV = 1000 * trapz(tPlot, abs(u));
    cleanupDV = cleanup_delta_v_equivalent(final(2), final(3), cfg);

    sol = struct();
    sol.method = "Reduced Lyapunov";
    sol.satId = sat.id;
    sol.slotId = slot.id;
    sol.t = tPlot(:);
    sol.a = xPlot(:,1);
    sol.relRAAN = arrayfun(@wrapToPiLocal, xPlot(:,2));
    sol.relM = arrayfun(@wrapToPiLocal, xPlot(:,3));
    sol.u = u(:);
    sol.mode = mode(:);
    sol.tCoast = plan.tCoast;
    sol.tControlled = max(0, tPlot(end) - plan.tCoast);
    sol.tf = tPlot(end);
    sol.arrivalEpoch = cfg.launchEpoch + seconds(sol.tf);
    sol.rawDeltaV_mps = rawDV;
    sol.cleanupDeltaV_mps = cleanupDV;
    sol.deltaV_mps = rawDV + cleanupDV;
    sol.finalAError_km = final(1) - slot.aF;
    sol.finalRAANError = wrapToPiLocal(final(2));
    sol.finalMError = wrapToPiLocal(final(3));
    sol.minAltitude_km = min(sol.a - cfg.Earth.radius);
    sol.maxAltitude_km = max(sol.a - cfg.Earth.radius);
    sol.plan = plan;

    if abs(sol.finalAError_km) <= cfg.lyap.tolA_km && ...
       abs(rad2deg(sol.finalRAANError)) <= cfg.raanTol_deg
        sol.success = true;
        sol.status = "Converged";
    else
        sol.success = false;
        sol.status = "SMA reached; residual slot error reported as cleanup DV";
    end
end

function xdot = reduced_lyap_ode(t, x, slot, plan, cfg)
    u = reduced_lyapunov_control(t, x, plan, cfg);
    xdot = starlink_reduced_elem_dyn(t, x, u, slot, cfg);
end

function [value, isterminal, direction] = reduced_lyap_terminal_event(t, x, plan, cfg)
    if t < plan.tCoast
        value = 1;              % keep integrating during coast
    else
        value = abs(x(1) - plan.aF) - cfg.lyap.tolA_km;
    end
    isterminal = 1;
    direction = -1;
end
