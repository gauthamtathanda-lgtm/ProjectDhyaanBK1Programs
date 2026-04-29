function H = computeDFA(x)
% computeDFA  Hurst Exponent via Detrended Fluctuation Analysis
%
% Usage:  H = computeDFA(x)
%   x — 1D signal vector (>= 200 samples recommended)
%
% Returns: H in (0,1)
%   H > 0.5 → persistent long-range correlations (meditative coherence)
%   H ~ 0.5 → random / uncorrelated
%   H < 0.5 → anti-persistent
%
% Relationship: HFD + H ≈ 2 (cross-validation check)

    x    = double(x(:));
    N    = length(x);
    y    = cumsum(x - mean(x));        % integrate mean-centred signal

    nMin = max(4, floor(N * 0.02));
    nMax = floor(N / 4);
    ns   = unique(round(logspace(log10(nMin), log10(nMax), 20)));
    ns   = ns(ns >= 4 & ns <= nMax);

    if numel(ns) < 4
        H = NaN; return;
    end

    Fn = zeros(size(ns));
    for i = 1:numel(ns)
        n       = ns(i);
        nBlocks = floor(N / n);
        rms2    = zeros(1, nBlocks);
        for b   = 1:nBlocks
            seg   = y((b-1)*n+1 : b*n);
            t     = (1:n)';
            A     = [t, ones(n,1)];
            c     = A \ seg;
            res   = seg - A*c;
            rms2(b) = mean(res.^2);
        end
        Fn(i) = sqrt(mean(rms2));
    end

    valid = Fn > 0 & ~isnan(Fn);
    if sum(valid) < 4
        H = NaN; return;
    end
    p = polyfit(log(ns(valid)), log(Fn(valid)), 1);
    H = p(1);
end