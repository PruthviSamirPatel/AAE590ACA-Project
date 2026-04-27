function n = circ_mean_rate(a, Earth)
%CIRC_MEAN_RATE Circular Keplerian mean motion [rad/s].
%
% Vectorized: a may be scalar or array. MATLAB's integral() passes
% vector-valued quadrature nodes, so use elementwise operations.
    n = sqrt(Earth.mu ./ (a.^3));
end
