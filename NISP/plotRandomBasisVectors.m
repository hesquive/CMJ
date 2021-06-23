function plotRandomBasisVectors(Psinum)
% Hugo Esquivel, 2021
% -

numRandomBasisVectors=length(Psinum);
numPoints=50;

x1=linspace(-1,1,numPoints); % x1 = [-1,1]
x2=linspace(-1,1,numPoints); % x2 = [-1,1]

[x1grid,x2grid]=meshgrid(x1,x2);

for k=1:numRandomBasisVectors
    Psikgrid=Psinum{k}(x1grid,x2grid)  +0*x1grid+0*x2grid;
    
    figure
    surf(x1grid,x2grid,Psikgrid)
    colormap(jet)
    title(sprintf('\\bfseries{Random Basis Vector $\\Psi_{%d}$}',k),'Interpreter','LaTeX')
    xlabel('$x_1$','Interpreter','LaTeX')
    ylabel('$x_2$','Interpreter','LaTeX')
    zlabel(sprintf('$\\Psi_{%d}(x_1,x_2)$',k),'Interpreter','LaTeX')
end
end
