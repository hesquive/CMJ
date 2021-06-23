function fapprox=getApproximateStructuralResponse(coeff,Psi)
% Hugo Esquivel, 2021
% -

numRandomBasisVectors=length(Psi);

fapprox=sym(0);

for k=1:numRandomBasisVectors
    fapprox=fapprox+coeff(k)*Psi{k};
end

fapprox=expand(vpa(fapprox));
end
