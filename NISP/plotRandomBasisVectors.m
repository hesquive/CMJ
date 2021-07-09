function plotRandomBasisVectors(Psinum,probabilityInfo)
% Hugo Esquivel, 2021
% -

numRandomVariables=length(probabilityInfo.name);

if numRandomVariables==2
    numRandomBasisVectors=length(Psinum);
    
    numPoints=50;
    
    s=getProbabilitySupports_Plot(probabilityInfo);
    
    x1=linspace(s(1,1),s(1,2),numPoints);
    x2=linspace(s(2,1),s(2,2),numPoints);
    
    [x1grid,x2grid]=meshgrid(x1,x2);
    
    for k=1:numRandomBasisVectors
        Psikgrid=Psinum{k}(x1grid,x2grid)  +0*x1grid+0*x2grid;
        
        figure
        surf(x1grid,x2grid,Psikgrid)
        colormap(turbo)
        title(sprintf('\\bfseries{Random Basis Vector $\\Psi_{%d}$}',k),'Interpreter','LaTeX')
        xlabel('$x_1$','Interpreter','LaTeX')
        ylabel('$x_2$','Interpreter','LaTeX')
        zlabel(sprintf('$\\Psi_{%d}(x_1,x_2)$',k),'Interpreter','LaTeX')
    end
end
end
