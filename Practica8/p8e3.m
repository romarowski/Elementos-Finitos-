% clear
% nodes_cart=[2 2;1 -5;4 3;-1.5 2.5];
% 
% syms x y c e;
% 
% func=[1 c e c*e];
% nodes=[-1 -1;1 -1; 1 1;-1 1];
% 
% for i=1:4
%    A(i,:) = subs(func,[c e],nodes(i,:)); 
% end
% 
% N=func/A;
% 
% x_loc = N*nodes_cart(:,1);
% y_loc = N*nodes_cart(:,2);
% 
% i=1; j=1;
% x_cart=zeros(232,1); y_cart=zeros(232,1);
% for n=1-1:0.1:1
%     for m=-1:0.1:1
%     x_cart(i,j)=double(subs(x_loc,[c e],[m n]));
% 
%     y_cart(i,j)=double(subs(y_loc,[c e],[m n]));
%     i=i+1;
%     end
%     j=j+1;
% end

%% Q8

clear
nodes_cart=[;1 -5;4 3;-1.5 2.5];

syms x y c e;

func=[1 c e c^2 c*e e^2 c^2*e c*e^2];
nodes=[-1 -1;1 -1; 1 1;-1 1;0 -1;1 0;0 1;-1 0];

for i=1:4
   A(i,:) = subs(func,[c e],nodes(i,:)); 
end

N=func/A;

x_loc = N*nodes_cart(:,1);
y_loc = N*nodes_cart(:,2);

i=1; j=1;
x_cart=zeros(232,1); y_cart=zeros(232,1);
for n=1-1:0.1:1
    for m=-1:0.1:1
    x_cart(i,j)=double(subs(x_loc,[c e],[m n]));

    y_cart(i,j)=double(subs(y_loc,[c e],[m n]));
    i=i+1;
    end
    j=j+1;
end