%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% Algoritmo de mapeo de retorno
% Ing. Sergio A. Merlino Chiozza
%
function [sn,epn,alfan,ftrial]=amr(epsilonn,ep,alfa,Ee,sy,mK)
strial=Ee*(epsilonn-ep);
ftrial=abs(strial)-sy-mK*alfa;
if ftrial<=0 % paso elástico
    sn=strial;
    epn=ep;
    alfan=alfa;
else % paso plástico (mapeo de retorno)
    dgamma=ftrial/(Ee+mK);
    sn=strial-dgamma*Ee*sign(strial);
    epn=ep+dgamma*sign(strial);
    alfan=alfa+dgamma;
end
end