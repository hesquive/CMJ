function Psi=getRandomBasisVectors(numRandomBasisVectorsPerDimension)
% Hugo Esquivel, 2021
% -
% Remarks:
% 1. Because the random variables are uniformly distributed, the orthogonal system to use over each random dimension is 
%    the "Legendre-Chaos system".
% 2. For simplicity, a full tensor product is used below to get the required random basis on the probability space.

syms x1 x2

numRandomVariables=2;
numRandomBasisVectors=numRandomBasisVectorsPerDimension^numRandomVariables;

Psi=cell(numRandomBasisVectors,1);

k=0;

for j=1:numRandomBasisVectorsPerDimension
    for i=1:numRandomBasisVectorsPerDimension
        k=k+1;
        
        Psi{k}(x1,x2)=legendreP(i-1,x1)*legendreP(j-1,x2);
    end
end
end
