function bathy = bathySpatialRegularization(bathy)
%BATHYSPATIALREGULARIZATION compute fCombined depth estimates
%
% bathy = BATHYSPATIALREGULARIZATION(bathy)
%

x = bathy.xm; y = bathy.ym;
nx = length(x); ny = length(y);
dx = bathy.params.dxm; dy = bathy.params.dym;
nF = size(bathy.fDependent.fB, 3);
[A,b] = build_Ab(bathy);
switch bathy.params.regx
    case 2
        Dx = build_D2x(nx,ny,dx);
    case 3
        Dx = build_D3x(nx,ny,dx);
    case 4
        Dx = build_D4x(nx,ny,dx);
end
switch bathy.params.regy
    case 2
        Dy = build_D2y(nx,ny,dy);
    case 3
        Dy = build_D3y(nx,ny,dy);
    case 4
        Dy = build_D4y(nx,ny,dy);
end
Gamma = vertcat(A,bathy.params.taux*Dx,bathy.params.tauy*Dy);
phi = vertcat(b,zeros(size(Dx,1)+size(Dy,1),1));
h = Gamma\phi;
r = (A*h-b).^2;
hErr = zeros(nx*ny,1);
[~,J] = find(A);
for i=1:nx*ny
    n = nansum(A(:,i))/max(A(:,i));
    variance = nansum(r(J==i).^2)/...
        nansum(A(:,i).^2);
    hErr(i) = sqrt(variance)*tinv(1-.05/2,n-1);
end
bathy.fCombined.h = reshape(h,ny,nx);
bathy.fCombined.hErr = reshape(hErr,ny,nx);
end


function D2x = build_D2x(nx,ny,dx)
stencil = diag([1 -2 1])*ones(3,nx*ny);
diagonals = [0,ny,2*ny];
D2x = (1/dx^2)*spdiags(stencil',diagonals,(nx-2)*ny,nx*ny);
end


function D2y = build_D2y(nx,ny,dy)
stencil = diag([1 -2 1])*ones(3,ny);
diagonals = [0,1,2];
d2y = (1/dy^2)*spdiags(stencil',diagonals,ny-2,ny);
I = speye(nx);
D2y = kron(I,d2y);
end


function D3x = build_D3x(nx,ny,dx)
stencil = diag([-1 3 -3 1])*ones(4,nx*ny);
diagonals = [0,ny,2*ny,3*ny];
D3x = (1/dx^3)*spdiags(stencil',diagonals,(nx-3)*ny,nx*ny);
end


function D3y = build_D3y(nx,ny,dy)
stencil = diag([-1 3 -3 1])*ones(4,ny);
diagonals = [0,1,2,3];
d3y = (1/dy^3)*spdiags(stencil',diagonals,ny-3,ny);
I = speye(nx);
D3y = kron(I,d3y);
end


function D4x = build_D4x(nx,ny,dx)
stencil = diag([1 -4 6 -4 1])*ones(5,nx*ny);
diagonals = [0,ny,2*ny,3*ny,4*ny];
D4x = (1/dx^4)*spdiags(stencil',diagonals,(nx-4)*ny,nx*ny);
end


function D4y = build_D4y(nx,ny,dy)
stencil = diag([1 -4 6 -4 1])*ones(5,ny);
diagonals = [0,1,2,3,4];
d4y = (1/dy^4)*spdiags(stencil',diagonals,ny-4,ny);
I = speye(nx);
D4y = kron(I,d4y);
end
