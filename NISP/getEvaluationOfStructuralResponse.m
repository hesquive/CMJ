function f=getEvaluationOfStructuralResponse(fact,x)
% Hugo Esquivel, 2021
% -

syms x1 x2

f=matlabFunction(fact,'Vars',[x1,x2]);
f=f(x(:,1),x(:,2)); % numerical evaluation of fnum using x as input
end
