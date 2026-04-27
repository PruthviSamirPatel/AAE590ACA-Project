function tRaise = compute_raise_time(a0, aF, uTangential, cfg)
%COMPUTE_RAISE_TIME Analytic circular low-thrust raise time.
%
% For circular tangential thrust:
%   da/dt = 2*u/n(a), n(a)=sqrt(mu/a^3)
%
% Integrating from a0 to aF with constant positive u:
%   t = sqrt(mu)/u * (1/sqrt(a0) - 1/sqrt(aF))

    if aF < a0
        sgn = -1;
    else
        sgn = 1;
    end

    u = abs(uTangential);
    if u <= 0
        error('compute_raise_time: uTangential must be nonzero.');
    end

    tRaise = sqrt(cfg.Earth.mu)/u * abs(1/sqrt(a0) - 1/sqrt(aF));

    if sgn < 0
        % same duration for lowering; sign is handled by the controller.
        tRaise = abs(tRaise);
    end
end
