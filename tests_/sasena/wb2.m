function wb = wb2(x, y_best, surr_model)
% y_best - best sampled value (for constrained cases, this is the best
% feasible sample)
% s - mean squared error estimate at x
% surr_model - dacefit surrogate model of objective function

% calculate the Expected Improvement (EI) avoiding mathematical underflow
[y_pred, ~, s2] = predictor(x, surr_model);
s = abs(s2) ^ 0.5;

a = (1/sqrt(2))*((y_best-y_pred)/s);

if a < -5  % if a << -1
    % MacLaurin series expansion to avoid floating point underflow
    term = zeros(1,20);
    for i = 1:20
        term(i) = ((-1)^(i-1)*(factorial(2*(i-1)) / ...
            (2^(i-1)*factorial(i-1)))/(2^(i-1)))*a^(-(2*(i-1)+1));
    end
    
    B = (y_best-y_best)*(1/(2*sqrt(pi)))*sum(term) + (s/sqrt(2*pi));
    C = exp(-(1/2)*a^2);
    EI = B * C;
    
else
    B = (y_best-y_best)*(0.5+0.5*erf(a));
    C = s*(1/sqrt(2*pi))*exp(-(1/2)*a^2);
    
    EI = B + C + realmin;
end

% wb is to be MINIMIZED
if s > eps
    % equation 4.10 - Sasena
    wb = y_pred - EI;
else
    wb = y_pred;
end

end