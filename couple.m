L = chebop(0, 10*pi);
L.op = @(x, u, v) [diff(u) - v; diff(v) + u];
L.lbc = @(u, v) [u-1; v];
rhs = [0; 0];
U = L\rhs;

clf, plot(U), grid on

xxx