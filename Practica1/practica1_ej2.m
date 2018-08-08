syms A E L;
hold on
L = 60;
A = 2;
E = 30e6;
t_x = @(x) -10 * x;
maxelem = 50;
for N=3:maxelem
    nodes = 0:L/N:L;
    elem = zeros(N,2);
    for i=1:N
        elem(i,:) = [i i+1];
    end
    k_g = zeros(N+1);
    fzas = zeros(1,N+1);
    for i=1:N
        k_e = A * E / L * [1 -1 ; -1 1];
        k_g(elem(i,:), elem(i,:)) = k_e + k_g(elem(i,:),elem(i,:));
        fzapernode = integral(t_x,nodes(i), nodes(i+1));
        fzas(i) = 1/3 * fzapernode;
        fzas(i+1) = 2/3 * fzapernode;
    end
    bc = zeros(1,N+1);
    bc(end) = 1;
    k_red = k_g(~bc,~bc);
    despl = k_red ^ (-1) * fzas(1:end-1)';
    figure(1)
    plot(nodes,[despl' , 0]);
end
 
    