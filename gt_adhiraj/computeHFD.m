function hfd = computeHFD(x, kmax)
% computeHFD  Higuchi (1988) Fractal Dimension
% VALIDATED: line=1.000, sinusoid~1.000, Gaussian noise~2.000, Brownian~1.500
%
% Usage:  hfd = computeHFD(x, kmax)
%   x    — 1D signal vector (double)
%   kmax — max lag: use 10 for broadband, 5 for narrow-band filtered signals
%
% Returns: hfd in range [1, 2]

    if nargin < 2, kmax = 10; end
    x  = double(x(:));
    N  = length(x);
    Lk = zeros(1, kmax);

    for k = 1:kmax
        Lm_sum = 0;
        for m = 1:k
            idx  = m:k:N;
            Xm   = x(idx);
            Mf   = length(Xm) - 1;
            if Mf < 1, continue; end
            nrm      = (N - 1) / (floor((N - m) / k) * k);
            Lm_sum   = Lm_sum + (sum(abs(diff(Xm))) / k) * nrm;
        end
        Lk(k) = Lm_sum / k;
    end

    valid = Lk > 0;
    if sum(valid) < 2
        hfd = NaN; return;
    end
    logk  = log(find(valid));          % log(k)
    logL  = log(Lk(valid));
    n     = sum(valid);
    slope = (n*sum(logk.*logL) - sum(logk)*sum(logL)) / ...
            (n*sum(logk.^2) - sum(logk)^2 + 1e-15);
    hfd   = -slope;                    % HFD is the negative slope
end