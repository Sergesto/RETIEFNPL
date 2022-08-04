%
% Algoritmo de mapeo de retorno
% Ing. Sergio A. Merlino Chiozza
%
function [sn,epn,alfan,ftrial]=amr(epsilonn,ep,alfa,Ee,sy,mK)
strial=Ee*(epsilonn-ep);
ftrial=abs(strial)-sy-mK*alfa;
if ftrial<=0 % paso el�stico
    sn=strial;
    epn=ep;
    alfan=alfa;
else % paso pl�stico (mapeo de retorno)
    dgamma=ftrial/(Ee+mK);
    sn=strial-dgamma*Ee*sign(strial);
    epn=ep+dgamma*sign(strial);
    alfan=alfa+dgamma;
end
end