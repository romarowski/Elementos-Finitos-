syms b a x y;
b=0.5; a=1;
%p= @(x,y) [1  x  y  x^2  x*y y^2]; LST
p= @(x,y) [1 x y x*y];
func=[1 x y x*y];
A=[p(0,0);p(a,0);p(a,b);p(0,b)];

N=func/A;


m=1;
for j=1:4
             B(1,j)=diff(N(m),x);
             m=m+1;
end
       
    
m=1;
for j=1:4
             B(2,j)=diff(N(m),y);
             m=m+1;
end
k_hat=[25 0;0 25];

K=int(int(B'*k_hat*B,x,-a,a),y,-b,b);

nodeDofs=[1 2 5 4;2 3 6 5;4 5 8 7;5 6 9 8];

k_g=zeros(9);
for i=1:4
    k_g(nodeDofs(i,:),nodeDofs(i,:))= double(K) +k_g(nodeDofs(i,:),nodeDofs(i,:));
    
end
c=[7 8 9]; x =[1 2 3 4 5 6];

k_xc=k_g(x,c); k_xx=k_g(x,x); k_cc=k_g(c,c); k_cx=k_xc';

flujo_c = [0 0 0 0 1000 0];
temp_c = [100 100 100]';
temp_x= inv(k_xx)*(flujo_c-k_xc*temp_c);

temp_tot= [temp_x; temp_c]

for i=1:4
T(i,:)=temp_tot(nodeDofs(i,:))
end


