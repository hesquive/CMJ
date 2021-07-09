function x=getRealizationsFromRandomVariables(numRealizations,probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

x=zeros(numRealizations,numRandomVariables);

for i=1:numRandomVariables
    if strcmpi(probabilityInfo.name{i},'Uniform') % uniform distribution (defined in [-1,1])
        x(:,i)=-1+2*rand(numRealizations,1);
        
    elseif strcmpi(probabilityInfo.name{i},'Beta') % beta distribution (defined in [-1,1])
        alpha=probabilityInfo.pars{i}(1);
        beta=probabilityInfo.pars{i}(2);
        
        x(:,i)=-1+2*betarnd(alpha,beta,numRealizations,1);
        
    elseif strcmpi(probabilityInfo.name{i},'Normal') % normal distribution (standard version)
        x(:,i)=normrnd(0,1,numRealizations,1);
    end
end
end
