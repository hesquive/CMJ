function Psi=getRandomBasisVectors(numRandomBasisVectorsPerDimension,probabilityInfo)
% Hugo Esquivel, 2021
% -
% Remark:
% - For simplicity, a full tensor product is used below to get the required random basis over the probability space.

numRandomVariables=length(probabilityInfo.name);
numRandomBasisVectors=numRandomBasisVectorsPerDimension^numRandomVariables;

pol=cell(numRandomBasisVectorsPerDimension,numRandomVariables);

syms x

for j=1:numRandomVariables
    for i=1:numRandomBasisVectorsPerDimension
        if strcmpi(probabilityInfo.name{j},'Uniform')
            pol{i,j}(x)=legendreP(i-1,x);
            
        elseif strcmpi(probabilityInfo.name{j},'Beta')
            alpha=probabilityInfo.pars{j}(1);
            beta=probabilityInfo.pars{j}(2);
            
            alphaJ=beta-1;
            betaJ=alpha-1;
            
            pol{i,j}(x)=jacobiP(i-1,alphaJ,betaJ,x);
            
        elseif strcmpi(probabilityInfo.name{j},'Normal')
            pol{i,j}(x)=hermiteH(i-1,x/sqrt(2))/2^((i-1)/2);
        end
    end
end

if numRandomVariables==2
    syms x1 x2
end

Psi=cell(numRandomBasisVectors,1);

k=0;

for j=1:numRandomBasisVectorsPerDimension
    for i=1:numRandomBasisVectorsPerDimension
        k=k+1;
        
        Psi{k}(x1,x2)=pol{i,1}(x1)*pol{j,2}(x2);
    end
end
end
