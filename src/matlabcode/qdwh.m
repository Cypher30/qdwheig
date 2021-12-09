function U = qdwh(A)
% QDWH algorithm to compute polar decomposition of A

% Initialization
U = A / normest(A, 3e-1);
L = 0.9 / condest(A);
[nrow, ncol] = size(A);
D = zeros(nrow + ncol, ncol);
D(nrow + 1:end, :) = eye(ncol);
% count = 1;
tol1 = 10 * eps / 2;
tol2 = tol1^(1 / 3);
if nrow == ncol && norm(A - A', 'fro') / norm(A, 'fro') < 100 * eps / 2
    symm = 1;
else
    symm = 0;
end

% Iteration
while 1
    LL = L^2;
    gamma = (4 * (1 - LL) / LL^2)^(1/3);
    sqg = sqrt(1 + gamma);
    a = sqg + 1 / 2 * sqrt(8 - 4 * gamma + 8 * (2 - LL) / (LL * sqg));
    b = (a - 1)^2 / 4;
    c = a + b - 1;
    sqc = sqrt(c);
    D(1:nrow, :) = sqc * U;
    [Q, ~] = qr(D, 0);
    quo = b / c;
    Utemp = quo * U + 1 / sqc * (a - quo) * Q(1:nrow, :) * Q(nrow + 1:end, :)';
    if norm(Utemp - U, 'fro') <= tol2 && abs(1 - L) <= tol1
        break;
    else
        if symm
            Utemp = (Utemp + Utemp') / 2;
        end
        U = Utemp;
        L = L * (a + b * LL) / (1 + c * LL);
%         count = count + 1;
    end
end
% count
if symm
    U = (U + U') / 2;
end
end