% program to find displacements and stresses in a truss us FEM
clear;
Data;
global j;
j = 3;
global nfun;
nfun = 0;

x = Elem(:,3);

[weight,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
[~,grad_w] = fungrad(x)
[~, ~,grad_s,~] = confungrad(x)
%[grad_w, grad_s] = complex_step(@Truss,delx,ndof,nbc,nelem,E,dens,Node,force,bc,Elem)


