clear;clc;

elementos=load('ElementosT1.txt');
nodos=load('NodosT1.txt');
nodos=nodos*1000;

%% Propiedades del Material
E=1;
nu=0.3;
%Plane Strain
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu      0;
                           nu  1-nu      0;
                            0    0  0.5-nu];
%% Definiciones
nDofNod = 2;                    % grados de libertad por nodo
nel = size(elementos,1);         % elementos
nNod = size(nodos,1);           % nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
%% Condiciones de borde
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(19:24,:) = true;
%% Cargas
R = zeros(nNod,nDofNod);        % Vector de cargas

% Puntos de Gauss
rsInt = 3*ones(1,1);
[wpg, upg, npg] = gauss(rsInt);


%Carga de Compresion superrior
for iele=6:10
    r=0;
    nodesEle=nodos(elementos(iele,:),:);
    for ipg=1:npg
    eta=1;
    ksi=upg(ipg);
    
    N=shapefuns([ksi eta],'Q4');
    
    dN=shapefunsder([ksi eta],'Q4');
    
    jac=dN*nodesEle;
    
    ene(1,1:2:7)=N;
    ene(2,2:2:8)=N; 
    
    q=zeros(4,2);
    q(3:4,2)=-7200;
    q=reshape(q',[8 1]);
    
    I=ene'*ene*q*jac(1,1);
    
    r=r+wpg(ipg)*I;
    
    
    
    end
    R(elementos(iele,:),:)= R(elementos(iele,:),:) + reshape(r,[2 4])';
end


