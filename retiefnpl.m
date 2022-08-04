% Estructuras Reticuladas 3D
% Elasticidad finita /Material de Saint Venant-Kirchhoff\
% Deformación de Green-Lagrange /Hipótesis de pequeñas deformaciones\
% Plasticidad /Algoritmo de mapeo de retorno\
% Ing. Sergio A. Merlino Chiozza
%
clc; clear; clf; close all;
%
% Matriz A Nodos, Coordenadas, Restricciones y Fuerzas:
% Restringido(0), Libre(1)
%
NNodos=3; NElem=2;
%
A=[1 0 0 0 0 0 0 0 0 0;
   2 1 1 0 1 1 0 0 0 0;
   3 2 0 0 0 0 0 0 0 0];
%
% Matriz B de Conectividad:
%
B=[1 1 2;
   2 2 3];
%
x = [0 1 2];
y = [0 1 0];
z = [0 0 0];
figure(1);
plot(x,y,'o-','markersize',8,'markerfacecolor','r','linewidth',2);
text(0.05,0.03,'A','fontsize',12);
text(0.88,0.97,'B','fontsize',12);
text(1.90,0.03,'C','fontsize',12);
title('Plasticidad / Reticulado de von Mises');
hold on;
%
% Intervalos de Carga/Descarga
%
Rango=[0,28000,28000,-28300,-28300,29000,29000,-29000];
%
% Valores iniciales
%
for u=1:NElem
ep(u)=0;
alfa(u)=0;
ftrial(u)=0;
end
%
%
deltaF=10;
deltad=zeros(NElem,NNodos);
deltadn=zeros(NNodos,3);
q=zeros(NNodos,3);
epsilonn=zeros(NElem,1);
%
% Propiedades de la Estructura
% Material y Geometría
%
% Largo L inicial de los elementos
%
for u=1:NElem
            n1=B(u,2);
            n2=B(u,3);
            xn1=A(n1,2);
            yn1=A(n1,3);
            zn1=A(n1,4);
            xn2=A(n2,2);
            yn2=A(n2,3);
            zn2=A(n2,4);
            L(u)=sqrt((xn1-xn2)^2+(yn1-yn2)^2+(zn1-zn2)^2);
end
%
% Área de las barras
%
Ae=[pi*0.010^2/4;
    pi*0.012^2/4];
%
% Módulo de Elasticidad de Young
% Tensión de fluencia / Módulo de plasticidad
% Densidad de fuerza de masa
%
Ee=210e9; sy=250e6; mK=21e8; rho=7850;
%
%
for rango=1:2:length(Rango)
% Inicio del proceso de carga, hasta el valor final FF
Fj=Rango(rango);
jj=1;
%
%
% Matriz A de nodos, coordenadas, restricciones y fuerzas:
%
A=[1 0 0 0 0 0 0 0 0  0;
   2 1 1 0 1 1 0 0 Fj 0;
   3 2 0 0 0 0 0 0 0  0];
%
% Carga final
%
FF=Rango(rango+1);
SS=cell(floor(Rango(rango+1)-Rango(rango))/deltaF+1,1);
ee=cell(floor(Rango(rango+1)-Rango(rango))/deltaF+1,1);
%
while Fj~=FF
%
clc; deltad=zeros(NElem,NNodos);
%
% Aún no se alcanzó el equilibrio (equi=0)
%
equi=1;
%
% Se carga con el valor de Fj actualizado
%
A(:,9)=[0 Fj 0]';
%
% Aquí se inicia la iteración (método del gradiente conjugado)
%
cont=0;
stop=0;
%
while stop==0
%
% Matriz Kr de la estructura (matriz reducida por condiciones de borde)
%
[Kr,VEDr,T,BB]=matrixfn(ftrial,Ae,L,Ee,mK,A,B,NElem,NNodos);
%
% Matriz de Fuerzas Nodales
%
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
            w=w+1;
            FN(w)=A(q,7+h);
        end
    end
end
%
if cont==0
%
% Posición inicial
%
g=[1 2];
%
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
           w=w+1;
           A(q,1+h)=g(w);
        end
    end
end
end
%
% Método del Gradiente Conjugado
%
[A,cont,stop]=gradconj(A,Kr,FN,VEDr,NNodos,cont);
%
end
%
% Matriz de Fuerzas Nodales
%
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
            w=w+1;
            FN(w)=A(q,7+h);
        end
    end
end
%
% Desplazamiento nodales
%
Q=Kr\FN';
w=0;
ww=0;
for n=1:NNodos
    for h=1:3
        ww=ww+1;
        if A(n,4+h)==1
            w=w+1;
            q(ww)=Q(w);
        else
            q(ww)=0;
        end
    end
end
%
while equi~=0
%
clc;
%
% /\/\/\ Deformaciones en cada elemento /\/\/\
% /\/\/\  Algoritmo de mapeo de retorno /\/\/\
%
for u=1:NElem
    n1=B(u,2);
    n2=B(u,3);
    epsilonn(u)=epsilonn(u)+BB{u}*[deltad(u,n1) deltad(u,n2)]';
    [s(u),ep(u),alfa(u),ftrial(u)]=amr(epsilonn(u),ep(u),alfa(u),Ee,sy,mK);
end
%
% Fuerzas internas para cada elemento
%
for u=1:NElem
    vfi(u,:)=[-1 1]*s(u)*Ae(u);
end
%
fiit=zeros(NNodos,3);
for u=1:NElem
    n1=B(u,2);
    n2=B(u,3);
    xn1=A(n1,2);
    yn1=A(n1,3);
    zn1=A(n1,4);
    xn2=A(n2,2);
    yn2=A(n2,3);
    zn2=A(n2,4);
    le=sqrt((xn1-xn2)^2+(yn1-yn2)^2+(zn1-zn2)^2);
    l=(xn2-xn1)/le;
    m=(yn2-yn1)/le;
    n=(zn2-zn1)/le;
    fiit(n1,:)=fiit(n1,:)+[vfi(u,1)*l vfi(u,1)*m vfi(u,1)*n];
    fiit(n2,:)=fiit(n2,:)+[vfi(u,2)*l vfi(u,2)*m vfi(u,2)*n];
end
%
fprintf(' |Análisis Plástico-Elasticidad finita\n');
fprintf(' |Material de Saint Venant-Kirchhoff\n');
fprintf(' |Deformación de Green-Lagrange\n');
fprintf(' |Algoritmo de mapeo de retorno\n');
fprintf('\n');
fprintf(' fuerzas nodales - internas\n');
disp(fiit);
fprintf(' alfa\n');
disp(alfa);
fprintf(' ep\n');
disp(ep);
fprintf(' epsilon\n');
disp(epsilonn');
fprintf(' fluencia\n');
disp(ftrial);
%
% Reacción
%
R=[T(1,:)*q' T(2,:)*q' T(3,:)*q'; 0 0 0 ; T(7,:)*q' T(8,:)*q' T(9,:)*q']-[A(1,8) A(1,9) A(1,10); 0 0 0; A(3,8) A(3,9) A(3,10)];
%
% Fuerzas externas
%
AR=[A(:,8) A(:,9) A(:,10)];
AR=[A(:,8) A(:,9) A(:,10)]+R;
if  all(abs(AR(2,:)-fiit(2,:))<=1e-10)
    disp('/\equilibrio/\');
    fprintf('\n');
    equi=0;
    if Fj<FF
        Fj=Fj+deltaF;
      else
        Fj=Fj-deltaF;
    end
else
        equi=1;
        w=0;
        for n=1:NNodos
            for h=1:3
                if A(n,4+h)==1
                    w=w+1;
                    FND(:,h)=fiit(n,h)-AR(n,h);
                end
            end
        end
%
deltadd=-Kr\FND';
%
        w=0;
        for n=1:NNodos
            for h=1:3
                if A(n,4+h)==1
                    w=w+1;
                        if abs(deltadd(w))<1e-18
                            deltadn(n,h)=0;
                        else
                            deltadn(n,h)=deltadd(w);
                        end
                else
                    deltadn(n,h)=0;
                end
            end
        end
% Tabla de cosenos directores
% le=longitud del elemento
% l, m, n cosenos directores de cada elemento
%
for t=1:NElem
    n1=B(t,2);
    n2=B(t,3);
    xn1=A(n1,2);
    yn1=A(n1,3);
    zn1=A(n1,4);
    xn2=A(n2,2);
    yn2=A(n2,3);
    zn2=A(n2,4);
    le=sqrt((xn1-xn2)^2+(yn1-yn2)^2+(zn1-zn2)^2);
    l=(xn2-xn1)/le;
    m=(yn2-yn1)/le;
    n=(zn2-zn1)/le;
    deltad(t,n1)=deltadn(n1,1)*l+deltadn(n1,2)*m+deltadn(n1,3)*n;
    deltad(t,n2)=deltadn(n2,1)*l+deltadn(n2,2)*m+deltadn(n2,3)*n;
end
end
end
%
    SS{jj}=s';
    ee{jj}=epsilonn';
    jj=jj+1;
%
end
%
    for t = 1:(jj-1)
        x1(t)=ee{t}(1);
        y1(t)=SS{t}(1);
    end
    for t = 1:(jj-1)
        x2(t)=ee{t}(2);
        y2(t)=SS{t}(2);
    end
        figure(2);
        plot(x1,y1,'LineWidth',1.5);
        title('Plasticidad / Barra 1 --algoritmo de retorno-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
        figure(3);
        plot(x2,y2,'LineWidth',1.5);
        title('Plasticidad / Barra 2 --algoritmo de retorno-- \sigma(\epsilon)');
        ylabel('Tensión \sigma');
        xlabel('Deformación \epsilon');
        hold on;
end
figure(2);
  legend({'I - Carga','II - Descarga/Carga','III - Descarga/Carga','IV - Descarga/Carga'},'location','southeast');
figure(3);
  legend({'I - Carga','II - Descarga/Carga','III - Descarga/Carga','IV - Descarga/Carga'},'location','southeast');
%
% /\Resultados/\
%
fprintf(' |Análisis Plástico-Elasticidad finita\n');
fprintf(' |Material de Saint Venant-Kirchhoff\n');
fprintf(' |Deformación de Green-Lagrange\n');
fprintf(' |Algoritmo de mapeo de retorno\n');
fprintf('\n');
% 
fprintf('Matriz K (reducida por condiciones de frontera)\n');
fprintf('\n');
    disp(Kr);
fprintf('Posición final:\n');
fprintf('\n'); 
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
           w=w+1;
           disp((A(q,1+h))');
        end
    end
end
fprintf('Desplazamiento vertical:\n');
fprintf('\n');
disp(1-A(2,3));
x = [0 A(2,2)/1.05 2]; % se amplifica el desplazamiento
y = [0 A(2,3)/1.05 0]; % para su visualización
z = [0 0 0];        
figure(1);
plot(x,y,'ro--','markersize',8,'markerfacecolor','b','linewidth',2);
title('Plasticidad / Reticulado de von Mises');
hold on;