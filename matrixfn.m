%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% Matriz Kr de la estructura (matriz reducida por condiciones de borde)
% Ing. Sergio A. Merlino Chiozza
%
function [Kr,VEDr,T,BB]=matrixfn(ftrial,Ae,L,Ee,mK,A,B,NElem,NNodos)
%
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
    C(t,1) = t;
    C(t,2) = le;
    C(t,3) = l;
    C(t,4) = m;
    C(t,5) = n;
end
%
% Matriz de deformación unitaria-desplazamiento del elemento
%
BB=cell(NElem,1);
for u = 1:NElem
     BB{u}=1./C(u,2).*[-1 1];
end
%
% Construcción de las matrices K de cada elemento
%
K = cell(NElem,1);
VEDe = cell(NElem,1);
for k = 1:NElem
         lee=C(k,2);
         le=C(k,3);
         me=C(k,4);
         ne=C(k,5);
         %
         % Derivada Segunda de la Energía de Deformación respecto de E
         % Componente S11 del Tensor de Cosserat para cada elemento
         %
         if ftrial(u)<=0
                        d2phie=Ee;
                        Se = Ee*0.5*(lee^2-L(k)^2)/L(k)^2;
                  else
                        d2phie=Ee*mK/(Ee+mK);
                        Se = Ee*mK/(Ee+mK)*0.5*(lee^2-L(k)^2)/L(k)^2;
         end
         %
         %
                ss=[(le^2)/2 le*me le*ne -le^2 -le*me -le*ne;
                    0 (me^2)/2 me*ne -le*me -me^2 -me*ne;
                    0 0 (ne^2)/2 -le*ne -me*ne -ne^2;
                    0 0 0 (le^2)/2 le*me le*ne;
                    0 0 0 0 (me^2)/2 me*ne;
                    0 0 0 0 0 (ne^2)/2];
         K{k}=Ae(k)/L(k)*d2phie*(lee^2/L(k)^2)*(ss+ss') + Ae(k)/L(k)*Se*[1 0 0 -1 0 0;
                                                                       0 1 0 0 -1 0;
                                                                       0 0 1 0 0 -1;
                                                                      -1 0 0 1 0 0;
                                                                       0 -1 0 0 1 0;
                                                                       0 0 -1 0 0 1];
%
%   Vector energía de deformación de cada elemento
%
VEDe{k}=Ae(k)*lee/L(k)*Se*[-le -me -ne le me ne];
%    
% 
end  
%
% Matriz K de la estructura, ensamble de las matrices elementales
%
T=zeros(NNodos*3);
for n=1:NElem
   for x=0:2
    for y=0:2
        T(B(n,2)*3-x,B(n,2)*3-y)=T(B(n,2)*3-x,B(n,2)*3-y)+K{B(n,1)}(3-x,3-y);
        T(B(n,2)*3-x,B(n,3)*3-y)=T(B(n,2)*3-x,B(n,3)*3-y)+K{B(n,1)}(3-x,6-y);
        T(B(n,3)*3-x,B(n,3)*3-y)=T(B(n,3)*3-x,B(n,3)*3-y)+K{B(n,1)}(6-x,6-y);
        T(B(n,3)*3-x,B(n,2)*3-y)=T(B(n,3)*3-x,B(n,2)*3-y)+K{B(n,1)}(6-x,3-y);
    end        
   end
end
%
% Vector energía VED de la estructura
% Ensamble de los vectores VEDe elementales
%
VED=zeros(NNodos*3);
for n=1:NElem
   for j=0:2
        VED(B(n,2)*3-j,1)=VED(B(n,2)*3-j,1) + VEDe{B(n,1)}(3-j);
        VED(B(n,3)*3-j,1)=VED(B(n,3)*3-j,1) + VEDe{B(n,1)}(6-j);
   end
end
%
% Vector Energía de Deformación VED Reducido
%
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
            w=w+1;
            VEDr(w,1)=VED(q*3-3+h,1);
        end
    end
end
clc;
%
% Matriz K reducida
%
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
            w=w+1;
            Kr1(w,:)=T(q*3-3+h,:);
        end
    end
end
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
            w=w+1;
            Kr(:,w)=Kr1(:,q*3-3+h);
        end
    end
end