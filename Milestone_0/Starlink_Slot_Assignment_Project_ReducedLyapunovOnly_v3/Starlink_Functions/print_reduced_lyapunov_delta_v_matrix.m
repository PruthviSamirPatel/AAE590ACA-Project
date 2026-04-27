function print_reduced_lyapunov_delta_v_matrix(DV, cfg)
%PRINT_REDUCED_LYAPUNOV_DELTA_V_MATRIX Print requested matrix.

    fprintf('\n================ TOTAL DELTA-V MATRIX [m/s] ================\n');
    fprintf('Rows = satellites, columns = slots.\n');
    fprintf('Values = raw Lyapunov raise Delta-V + residual cleanup equivalent Delta-V.\n');
    fprintf('Assignment objective applies %.1f %% margin to each selected transfer.\n\n', ...
        100*cfg.assignment.dvMarginFraction);

    print_matrix_with_labels(DV.total_mps, cfg, 'Reduced Lyapunov total Delta-V matrix [m/s]');

    fprintf('\n--- Raw Lyapunov Delta-V only [m/s] ---\n');
    print_matrix_body(DV.raw_mps, cfg);

    fprintf('\n--- Cleanup equivalent Delta-V only [m/s] ---\n');
    print_matrix_body(DV.cleanup_mps, cfg);

    fprintf('\n--- Final RAAN error matrix [deg] ---\n');
    print_matrix_body(DV.raanErr_deg, cfg);

    fprintf('\n--- Final mean-anomaly residual matrix [deg] ---\n');
    print_matrix_body(DV.mErr_deg, cfg);
end

function print_matrix_with_labels(A, cfg, titleStr)
    fprintf('--- %s ---\n', titleStr);
    print_matrix_body(A, cfg);
end

function print_matrix_body(A, cfg)
    fprintf('%14s', '');
    for j = 1:numel(cfg.slots)
        fprintf('%12s', cfg.slots(j).id);
    end
    fprintf('\n');

    for i = 1:numel(cfg.sats)
        fprintf('%14s', cfg.sats(i).id);
        for j = 1:numel(cfg.slots)
            fprintf('%12.4f', A(i,j));
        end
        fprintf('\n');
    end
end
