function coeff=getCoefficientsOfApproximateStructuralResponse(x,f,Psinum)
% Hugo Esquivel, 2021
% -

numRandomBasisVectors=length(Psinum);

coeff=zeros(numRandomBasisVectors,1);

for k=1:numRandomBasisVectors
    coeff(k)=sum(f.*Psinum{k}(x(:,1),x(:,2)))/sum(Psinum{k}(x(:,1),x(:,2)).^2);
end
end
