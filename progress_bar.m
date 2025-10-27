% Text-based progress bar used to showcase the progress of the current iteration
% to the user. 

function progress_bar(operation_name, i, N, start_time)

    persistent last_length last_update
    bar_length = 50;
    update_interval = 0.2;

    % First call initialization.
    if isempty(last_update)
        last_update = 0;
    end

    % Clamp iteration index.
    i = max(0, min(i, N));
    frac = i / max(1, N);

    % Check time since last update.
    elapsed = toc(start_time);
    if (elapsed - last_update < update_interval) && (i < N)
        return  % Skip update if called too soon.
    end
    last_update = elapsed;

    % Build bar and times.
    filled = round(frac * bar_length);
    bar = [repmat('#', 1, filled) repmat('-', 1, bar_length-filled)];
    pct = frac * 100;

    if i > 0
        estTotal = elapsed / frac;
        eta = estTotal - elapsed;
    else
        eta = NaN;
    end

    % Simple elapsed and ETA formatting.
    elapsed_str = format_time(elapsed);
    eta_str     = format_time(eta);

    % Display string.
    msg = sprintf('%s: [%s] %5.1f%% | Iterations: %3d/%d | Elapsed: %s | ETA: %s', ...
                  operation_name, bar, pct, i, N, elapsed_str, eta_str);

    % Backspace erase.
    if ~isempty(last_length)
        fprintf(repmat('\b', 1, last_length));
    end

    fprintf('%s', msg);
    drawnow limitrate
    last_length = length(msg);

    % Cleanup.
    if i == N
        fprintf('\n');
        clear last_length last_update
    end
end

function s = format_time(t)
    if isnan(t)
        s = 'Estimating...';
    elseif t < 60
        s = sprintf('%4.2fs', t);
    elseif t < 3600
        s = sprintf('%4.2fmin', t/60);
    else
        s = sprintf('%4.2fh', t/3600);
    end
end