function qdwheig_test(N)
% Test function for qdwh

% Initialization
H = rand(N, N);
H = H' + H;

% Eigvalue problem solved by MATLAB function eig and qdwh-eig
[V0, D0] = eig(H);
[V1, D1] = sdc(H);

% Outcome comparing
t = tiledlayout(2, 1);
nexttile;
imagesc(log10(abs(V0' * V0 - eye(N))));
colorbar;
axis square;
nexttile;
imagesc(log10(abs(V1' * V1 - eye(N))));
title(t, "visualized orthogonality of computed eigvector");
colorbar;
axis square;

normH = norm(H, 'fro');
er0 = norm(H - V0 * D0 * V0', 'fro') / normH;
er1 = norm(H - V1 * D1 * V1', 'fro') / normH;
or0 = norm(V0' * V0 - eye(N), 'fro') / sqrt(N);
or1 = norm(V1' * V1 - eye(N), 'fro') / sqrt(N);
fprintf("eig relative error: %d\n sdc relative eroor: %d\n", er0, er1);
fprintf("eig orthogonality: %d\n sdc orthogonality: %d\n", or0, or1);