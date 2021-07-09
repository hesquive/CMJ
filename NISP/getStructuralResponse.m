function f=getStructuralResponse(probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

if numRandomVariables==2
    syms x1 x2
end

p1=legendreP(2,x1)+jacobiP(3,3,5,x1);
p2=legendreP(5,x2)+jacobiP(2,2,7,x2);

% Structural response of interest (notice that this f is a function of x1 and x2):
f=p1*p2/1e2;

for i=1:numRandomVariables
    if strcmpi(probabilityInfo.name{i},'Normal')
        f=f/50;
    end
end
end
