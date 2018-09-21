clear 

clear

syms x y func;
func=[1 x y x*y];
nodes=[-1 -.5;1 -.5;1 .5;-1 .5];
A=sym('A',[4,4]);

for i=1:4
   A(i,:) = subs(func,[x y],nodes(i,:)); 
end

N=func/A;

B(1,:)=diff(N,x);
B(2,:)=diff(N,y);

k_h=[25 0;0 25];

k_loc = int(int(B'*k_h*B,y,-.5,.5),x,-1,1);

nodesElem= [1 2 5 4;2 3 6 5;5 6 9 8;4 5 8 7];
k_g=zeros(9);
for i=1:4
   k_g(nodesElem(i,:),nodesElem(i,:)) =  k_g(nodesElem(i,:),nodesElem(i,:)) + double(k_loc);
end

% c=[7 8 9]; x=[1 2 3 4 5 6];
% 
% qc=[0 0 0 0 1000 0]'; Tc=[100 100 100]';
% k_xx=k_g(x,x); k_xc=k_g(x,c); k_cc=k_g(c,c); k_cx=k_g(c,x);
% 
% T=zeros(9,1);
% T(x)=k_xx^-1*(qc-k_xc*Tc);
% T(c)=Tc;
% 
% nodes=[-2 -1;0 -1;2 -1;-2 0;0 0;2 0;-2 1;0 1;2 1];
% variable=zeros(4);
% for i=1:4
%     variable(i,:) = T(nodesElem(i,:));
% end
% 
% 
% q=zeros(9,1);
% q(x)=qc;
% q(c)=k_cc*Tc+k_cx*T(x);
% 
% variable2=zeros(4);
% for i=1:4
%     variable2(i,:) = q(nodesElem(i,:));
% end

%% parte B
fq=1000*int(int(N,y,-.5,.5),x,-1,1);

c=[7 8 9]; x=[1 2 3 4 5 6];

qc=zeros(9,1);
for i=1:4
    qc(nodesElem(i,:))= qc(nodesElem(i,:)) +[500; 500; 500; 500];
end
qc=qc(1:6); Tc=[100 100 100]';
k_xx=k_g(x,x); k_xc=k_g(x,c); k_cc=k_g(c,c); k_cx=k_g(c,x);

T=zeros(9,1);
T(x)=k_xx^-1*(qc-k_xc*Tc);
T(c)=Tc;

nodes=[-2 -1;0 -1;2 -1;-2 0;0 0;2 0;-2 1;0 1;2 1];
variable=zeros(4);
for i=1:4
    variable(i,:) = T(nodesElem(i,:));
end


q=zeros(9,1);
q(x)=qc;
q(c)=k_cc*Tc+k_cx*T(x);

variable2=zeros(4);
for i=1:4
    variable2(i,:) = q(nodesElem(i,:));
end
