% Método del Gradiente Conjugado
% Ing. Sergio A. Merlino Chiozza
%
function [A,cont,stop]=gradconj(A,Kr,FN,VEDr,NNodos,cont)
  %
% Método del Gradiente Conjugado
%
xk=(zeros(size(FN)))';
r=-(FN'-VEDr);
p=-r;
stopgc=0;
while norm(r)>1e-10 && stopgc==0
    alfa=(r'*r)/(p'*Kr*p);
    xk=xk+alfa*p;
    r=r+alfa*Kr*p;
    beta=(r'*r)/((r-alfa*Kr*p)'*(r-alfa*Kr*p));
    p=-r+beta*p;
    if p'*Kr*p < 0
        stopgc=1;
    end
end
pk=xk;
alfak=1;
w=0;
for q=1:NNodos
    for h=1:3
        if A(q,4+h)==1
           w=w+1;
           A(q,1+h)=A(q,1+h) + alfak*pk(w);
        end
    end
end
%
% Condición de parada
%
R=FN-VEDr';
DX=Kr\R';
    rk=Kr*alfak*pk+VEDr-FN';
    if norm(DX)<1e-12
        stop=1;
        else stop=0;
    end
    cont=cont+1;