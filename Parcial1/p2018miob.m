a=25; E=79e3; Iz=a^4/12; A=a^2;
nodeDofs=[1 2;3 4;5 6;7 8];
elem=[1 2;2 3;3 4;4 1;1 3];
nodes=[0 0; 1e3 0; 1e3 .5e3;0 .5e3];
dofs=8;

k_g=zeros(dofs);
for i=1:5
    ubic= [nodeDofs(elem(i,1),:) nodeDofs(elem(i,2),:)];
    Le=norm( nodes(elem(i,2),:) -nodes(elem(i,1),:));
    X = E*A/Le;
    Kel = X*[1 -1;-1 1];
    dts=nodes(elem(i,2),:) -nodes(elem(i,1),:)/Le;
    T=[dts 0 0;0 0 dts];
    k_g(ubic,ubic) = T'*Kel*T + k_g(ubic,ubic);

end

    x=[3 7 5]; c=[1 2 4 6 8];
    k_cc=k_g(c,c); k_cx= k_g(c,x); k_xx=k_g(x,x);
    d_c=[0  0  0  -1.5 -1.5]';
    d_x=k_xx^(-1)*(-k_cx'*d_c);
    r_x=k_cc*d_c+k_cx*d_x;

    despl=zeros(8,1);
    despl(x)=d_x; despl(c)=d_c;

    fzas=zeros(8,1);
    fzas(c)=r_x; 

fza_aprox=fzas(8)+fzas(6);
