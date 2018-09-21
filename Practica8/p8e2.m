h=20; Tinf=50;

syms x y func;
func=[1 x y x*y];
nodes=[-1 -1;1 -1;1 1;-1 1];
A=sym('A',[4,4]);

for i=1:4
   A(i,:) = subs(func,[x y],nodes(i,:)); 
end

N=func/A;

B(1,:)=diff(N,x);
B(2,:)=diff(N,y);

k_h=[25 0;0 25];

k_loc = int(int(B'*k_h*B ,y,-1,1),x,-1,1) +  int(subs(h*(N'*N),x,1),y,-1,1);


q=int(subs(h*50*(N'),x,1),y,-1,1);

c=[1 4]; x=[2 3];
Tc=[100 100]';
k_xx=k_loc(x,x); k_xc=k_loc(x,c); k_cc=k_loc(c,c); k_cx=k_loc(c,x);

qc=q(x);
T=zeros(4,1);
T(x)=k_xx^-1*(qc-k_xc*Tc);
T(c)=Tc;

