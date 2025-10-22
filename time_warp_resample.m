function y = time_warp_resample(t, v, tempo_scale, new_len)
    % time-warp by resampling to new_len using timestamps
    % new duration = old_duration / tempo_scale
    t = t(:); v = v(:);
    % ensure unique, monotonic sample points
    [t, ia] = unique(t, 'stable');
    v = v(ia);
    % build new time grid compressed/expanded by tempo_scale
    t0 = t(1);
    t1 = t(end);
    dur = (t1 - t0) / tempo_scale;
    t_new = linspace(t0, t0 + dur, new_len).';
    y = interp1(t, v, t_new, 'pchip', 'extrap');
end