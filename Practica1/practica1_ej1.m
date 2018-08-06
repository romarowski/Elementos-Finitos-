syms A_i A_s L P E delta_L;
hold on
A_i = 100;
delta_L = -0.3251;
A_s = 25;
L = 4000;
E = 210e3;
P = -1e3;
f = @(x) P/E * L * (log(A_i*A_i)- 2*log((A_i*(L-x)+x*A_s)/L)) / 2* (A_i -A_s);
%fplot(f,[0,L]);
maxelem = 50;
erramax = zeros(1,maxelem -2);
for N=3:maxelem % N es la cantidad de elementos
    nodes= 0:L/N:L;  % me dice posicion de los nodos, tiene una casillero mas que la cant de elementos
    elem =zeros(N,2);
    for i=1:N
        elem(i,:) = [i i+1];
    end
    areaele = @(x) A_i * (L-x)/L + A_s * x/L;
    k_g = zeros(N+1);
    for i=1:N
        k_e = areaele(nodes(i+1)) * E / (L/N) * [1 -1 ; -1 1];
        k_g(elem(i,:), elem(i,:)) =  k_e + k_g(elem(i,:), elem(i,:));   
    end
    bc = zeros(1,N+1);
    bc(1) = 1;
    fzas = zeros(1,N);
    fzas(end) = P;
    k_red = k_g(~bc,~bc);
    despl = k_red^(-1) * fzas'; % finalmente la mattriz de despl tiene la misma cantidad de elem que nodos porque
                                % el primer nodo esta fijo y despl=0
    erramax(N-2) = (delta_L - despl(N)) / delta_L;
    figure(1)
    plot([0,despl'],nodes);
    
end