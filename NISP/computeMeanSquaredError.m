function eps=computeMeanSquaredError(x,f,g,probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

if numRandomVariables==2
    syms x1 x2
end

numRealizations=size(x,1);

g=matlabFunction(g,'Vars',[x1,x2]);

eps=sum((f-g(x(:,1),x(:,2))).^2)/numRealizations;
end
