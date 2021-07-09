function s=getProbabilitySupports_Plot(probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

s=zeros(numRandomVariables,2);

for i=1:numRandomVariables
    s(i,:)=[-1,1];
    
    if strcmpi(probabilityInfo.name{i},'Normal')
        s(i,:)=3*s(i,:);
    end
end
end
