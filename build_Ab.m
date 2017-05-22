function [A,b] = build_Ab(bathy)
skill = bathy.fDependent.skill; lam = bathy.fDependent.lam1;
k = bathy.fDependent.k; fB = bathy.fDependent.fB;
g = 9.81;
params = bathy.params;
w = skill.*lam;
psi = (1./k).*atanh(4*pi*pi*fB.*fB./(g*k));
indices = (skill>params.QTOL & lam>params.minLam & ...
    psi>params.MINDEPTH & psi<params.MAXDEPTH);
sw = size(w);
w = reshape(w,sw(1)*sw(2),sw(3))';
psi = reshape(psi,sw(1)*sw(2),sw(3))';
indices = reshape(indices,sw(1)*sw(2),sw(3))';
w = w(indices); psi = psi(indices);
nx = length(bathy.xm); ny = length(bathy.ym); nw = length(w);
[~,J] = find(indices); ind = nw*(J-1)+(1:nw)';
A = spalloc(nw,nx*ny,nw); b = w.*psi;
A(ind)=w;
end