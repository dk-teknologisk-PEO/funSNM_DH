function [sigma_points, Wm, Wc] = generate_sigma_points(x, P, alpha, beta, kappa)
%GENERATESIGMAPOINTS Generates scaled sigma points for the UKF.
    n = length(x);
    lambda = alpha^2 * (n + kappa) - n;
    
    sigma_points = zeros(n, 2*n + 1);
    
    try
        % Use Cholesky decomposition. P must be positive semi-definite.
        U = chol((n + lambda) * P);
    catch
        % If P is not pos-def, use SVD as a robust fallback.
        [u, s, v] = svd((n + lambda) * P);
        U = u * sqrt(s);
    end

    sigma_points(:, 1) = x;
    for i = 1:n
        sigma_points(:, i + 1)   = x + U(i, :)';
        sigma_points(:, i + n + 1) = x - U(i, :)';
    end
    
    % Calculate weights for mean and covariance
    Wm = zeros(1, 2*n + 1);
    Wc = zeros(1, 2*n + 1);
    
    Wm(1) = lambda / (n + lambda);
    Wc(1) = lambda / (n + lambda) + (1 - alpha^2 + beta);
    
    for i = 2:(2*n + 1)
        Wm(i) = 1 / (2 * (n + lambda));
        Wc(i) = 1 / (2 * (n + lambda));
    end
end