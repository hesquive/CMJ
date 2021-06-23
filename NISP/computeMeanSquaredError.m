function eps=computeMeanSquaredError(x,f,fapprox)
% Hugo Esquivel, 2021
% -

syms x1 x2

numRealizations=size(x,1);

fapprox=matlabFunction(fapprox,'Vars',[x1,x2]);

eps=sum((f-fapprox(x(:,1),x(:,2))).^2)/numRealizations;
end
