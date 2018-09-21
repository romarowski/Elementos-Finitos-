clear 

E=30e6; nu=.25; t=.5;

syms x y func;
func=[1 x y x*y];
nodes=[-1.5 -1;1.5 -1;1.5 1;-1.5 1];
A=sym('A',[4,4]);

for i=1:4
   A(i,:) = subs(func,[x y],nodes(i,:)); 
end

N=func/A;


B=sym(zeros(3,8));
B(1,1:2:8) = diff(N,x);
B(2,2:2:8) = diff(N,y);
B(3,1:2:8) = diff(N,y); B(3,2:2:8) = diff(N,x);

N_c=sym(zeros(2,8));
N_c(1,1:2:8)=N; N_c(2,2:2:8)=N;

C=E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];

K=int(int(B'*C*B*.5,y,-1,1),x,-1.5,1.5);

q=zeros(8,1); q(6)=-300*.5;q(8)=-300*.5;

N_f=subs(N_c,y,1);
fzas=int(N_f'*N_f*q,x,-1.5,1.5); 

fzas(6) = fzas(6)-1000;

bc=zeros(8,1); bc(1:2)=[1 1]; bc(7:8)=[1 1]; bc(4)=1;

k_red=double(K(~bc,~bc));
fzas_red=double(fzas(~bc));

despl=zeros(8,1);
despl(~bc)= fzas_red'/k_red;


