function u = reduced_lyapunov_control(t, x, plan, cfg)
%REDUCED_LYAPUNOV_CONTROL Scalar tangential Lyapunov feedback.
%
% During coast:
%   u = 0
%
% During controlled phase:
%   V = 1/2 (a - aF)^2
%   u = sat[-kA (a-aF), +/- uMax]
%
% Since a_dot = 2u/n(a), this gives
%   V_dot = (a-aF) 2u/n(a) <= 0
% whenever the control is not blocked by safety logic.

    if t < plan.tCoast
        u = 0;
        return
    end

    aRef = plan.aF;
    eA = x(1) - aRef;

    u = -cfg.lyap.kA * eA;
    u = max(min(u, cfg.uMax_km_s2), -cfg.uMax_km_s2);

    % Do not lower below the drift shell.
    if x(1) <= cfg.aMin && u < 0
        u = 0;
    end

    % Avoid overshoot above the final shell in this scenario.
    if x(1) >= aRef && u > 0
        u = 0;
    end
end
