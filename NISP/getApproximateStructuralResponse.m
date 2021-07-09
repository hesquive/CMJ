function g=getApproximateStructuralResponse(coeff,Psi)
% Hugo Esquivel, 2021
% -

numRandomBasisVectors=length(Psi);

g=sym(0);

for k=1:numRandomBasisVectors
    g=g+coeff(k)*Psi{k};
end

g=expand(vpa(g));
end
