function [ grad_weight, grad_stress ] = complex_sum(obj, delx, ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

n=length(Elem);

for ii = 1:n
    Elem_new = Elem;
    Elem_new(ii,3) = Elem_new(ii,3) + 1i*delx;
    [weight_new, stress_new] = obj(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_new);
    grad_weight(ii,:) = imag(weight_new)/delx;
    grad_stress(ii,:) = imag(stress_new)/delx;
end

end


