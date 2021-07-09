function plotActualStructuralResponse(fact,probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

if numRandomVariables==2
    fact=matlabFunction(fact);
    
    numPoints=50;
    
    s=getProbabilitySupports_Plot(probabilityInfo);
    
    x1=linspace(s(1,1),s(1,2),numPoints);
    x2=linspace(s(2,1),s(2,2),numPoints);
    
    [x1grid,x2grid]=meshgrid(x1,x2);
    factgrid=fact(x1grid,x2grid);
    
    figure
    surf(x1grid,x2grid,factgrid)
    colormap(turbo)
    title('\bfseries{Actual Structural Response, $f$}','Interpreter','LaTeX')
    xlabel('$x_1$','Interpreter','LaTeX')
    ylabel('$x_2$','Interpreter','LaTeX')
    zlabel('$f(x_1,x_2)$','Interpreter','LaTeX')
end
end
