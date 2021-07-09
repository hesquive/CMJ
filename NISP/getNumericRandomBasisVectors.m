function Psinum=getNumericRandomBasisVectors(Psi,probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

if numRandomVariables==2
    syms x1 x2
end

numRandomBasisVectors=length(Psi);

Psinum=cell(numRandomBasisVectors,1);

for i=1:numRandomBasisVectors
    if i==1
        Psinum{i}=@(x1,x2) 1.0  +0*x1+0*x2;
    else
        Psinum{i}=matlabFunction(Psi{i},'Vars',[x1,x2]);
    end
end
end
