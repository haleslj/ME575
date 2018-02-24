function [ grad_weight, grad_stress ] = central_difference( obj, delx, ndof, nbc, nelem, E, dens, Node, force, bc, Elem)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% [stress0, weight0] =  obj(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
n = length(Elem);
    for ii = 1:n
        Elem_left = Elem;
        Elem_left(ii,3) = Elem_left(ii,3) - delx;
        [weight_left, stress_left] = obj(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_left);
       
        Elem_right = Elem;
        Elem_right(ii,3) = Elem_right(ii,3) + delx;
        [weight_right, stress_right] = obj(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_right); 
        
        grad_weight(ii,:) = (weight_right - weight_left)/(2*delx);
        grad_stress(ii,:) = (stress_right - stress_left)/(2*delx);
    end
end

