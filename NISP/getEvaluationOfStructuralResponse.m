function f=getEvaluationOfStructuralResponse(fact,x,probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

if numRandomVariables==2
    syms x1 x2
end

f=matlabFunction(fact,'Vars',[x1,x2]);
f=f(x(:,1),x(:,2)); % numerical evaluation of fact using x as input
end
