function [logs, log_n] = append_log(logs, log_n, S)
    log_n = log_n + 1;
    if log_n > numel(logs), logs = [logs, logs]; end
    fns = fieldnames(logs);
    for f = 1:numel(fns)
        fn = fns{f};
        if isfield(S,fn), logs(log_n).(fn) = S.(fn); end
    end
end