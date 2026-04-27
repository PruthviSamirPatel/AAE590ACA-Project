function writematrix_with_labels(A, cfg, filename)
%WRITEMATRIX_WITH_LABELS Save matrix as CSV with row/column labels.

    fid = fopen(filename, 'w');
    if fid < 0
        error('Could not open %s for writing.', filename);
    end

    fprintf(fid, 'Satellite/Slot');
    for j = 1:numel(cfg.slots)
        fprintf(fid, ',%s', cfg.slots(j).id);
    end
    fprintf(fid, '\n');

    for i = 1:numel(cfg.sats)
        fprintf(fid, '%s', cfg.sats(i).id);
        for j = 1:numel(cfg.slots)
            fprintf(fid, ',%.10f', A(i,j));
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
end
