function [U1, U2] = subspaceit(U, U0)
% Subspace iteration for computing the invarian subspace of U.

% Initialization
n = length(U);
k = round(norm(U, 'fro')^2);
if nargin < 2
    UU = U * randn(n, min(k + 3, n));
else
    UU = U * U0;
end

% Subspace iteration
[UU, ~] = qr(UU, 0);
UU = U * UU;
[UU, ~] = qr(UU);

U1 = UU(:, 1:k);
U2 = UU(:, k + 1:end);
end