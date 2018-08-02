syms A_i A_s L P E;
A_i = 100;
A_s = 25;
L = 4000;
E = 210e3;
P = -1e3;
f = @(x) P/E * L * (log(A_i*A_i)- 2*log((A_i*(L-x)+x*A_s)/L)) / 2* (A_i -A_s);
fplot(f,[0,L]);
