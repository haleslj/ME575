function [ grad_weight, grad_stress] = forward_difference(obj, delx, ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
n=length(Elem);
[weight0, stress0] = obj(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

for ii = 1:n
    Elem_new = Elem;
    Elem_new(ii,3) = Elem_new(ii,3) + delx;
    [weight_new, stress_new] = obj(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_new);
    grad_weight(ii,:) = (weight_new - weight0)/delx;
    grad_stress(ii,:) = (stress_new - stress0)/delx;

end
end

