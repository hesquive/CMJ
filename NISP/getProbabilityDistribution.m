function probabilityInfo=getProbabilityDistribution(probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

syms x [1,numRandomVariables]

pdf=cell(numRandomVariables,1);
supp=cell(numRandomVariables,1);

for i=1:numRandomVariables
    if strcmpi(probabilityInfo.name{i},'Uniform') % uniform distribution (standard version)
        pdf{i}(x(i))=1/2  +0*x(i);
        supp{i}=[-1,1];
        
    elseif strcmpi(probabilityInfo.name{i},'Beta') % beta distribution (standard version)
        alpha=probabilityInfo.pars{i}(1);
        beta=probabilityInfo.pars{i}(2);
        
        pdf{i}(x(i))=(x(i)+1)^(alpha-1)*(1-x(i))^(beta-1)/(2^(alpha+beta-1)*betafun(alpha,beta));
        supp{i}=[-1,1];
        
    elseif strcmpi(probabilityInfo.name{i},'Normal') % normal distribution (standard version)
        pdf{i}(x(i))=exp(-x(i)^2/2)/sqrt(2*sym(pi));
        supp{i}=[-Inf,Inf];
    end
end

probabilityInfo.pdf=pdf;
probabilityInfo.supp=supp;
end

function z=betafun(x,y)
z=beta(x,y);
end
