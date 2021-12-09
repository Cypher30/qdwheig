function [V, D] = sdc(H)
%
% Initialization
M = length(H);
shift = median(diag(H));
tol = 10 * eps / 2;
block = 32;
if M <= block
    [V, D] = eig(H);
    [eigs, Ix] = sort(diag(D), 'descend');
    V = V(:, Ix);
    D = diag(eigs);
    return
end

% Polar decomposition using QDWH
U = qdwh(H - shift * eye(M));

% Subspace Iteration
normH = norm(H, 'fro');
U = (U + eye(M)) / 2;
[U1, U2] = subspaceit(U);
re = norm(U2' * H * U1, 'fro') / normH;
if re > tol
    % Second subspace iteration
    [U1, U2] = subspaceit(U, U1);
    re = norm(U2' * H * U1, 'fro') / normH;
end
% If after two iteration, the accuracy cannot reach our demand, 
% we simply redo it
if re > tol
    for i = 1:2
        [U1t, U2t] = subspaceit(U);
        if norm(U2t' * H * U1t, 'fro') / normH < re
            U1 = U1t;
            U2 = U2t;
        end
    end
end

% After spliting the spectrum, we work on each subspace
% First, deal with the case that U1 or U2 has only single coloumn
% We could simply find the correponding eigenvalue
eigvals = zeros(M, 1);
if length(U1(1, :)) == 1
    eigvals(1) = U1' * H * U1;
    key = 1;
end
if length(U2(1, :)) == 1
    eigvals(M) = U2' * H * U2;
end

% Second, if U1 or U2 not in the case above, we recursively work on the
% corresponding subspace
if length(U1(1, :)) > 1
    [Ua, eigval1] = sdc(U1' * (H * U1));
    U1 = U1 * Ua;
    key = length(U1(1, :));
    eigvals(1:key) = diag(eigval1);
end
if length(U2(1, :)) > 1
    [Ua, eigval2] = sdc(U2' * (H * U2));
    U2 = U2 * Ua;
    eigvals(key + 1:end) = diag(eigval2);
end
V = [U1, U2];
V = 3/2*V-V*(V'*V)/2;
D = diag(eigvals);
end