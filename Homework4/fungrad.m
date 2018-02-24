function [ f, grad_weight ] = fungrad( x )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
global j;
global nfun;

delx = 1e-4;
n = length(x);
Data;
Elem(:,3) = x;
[f,~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
nfun = nfun + 1;

if j == 1
    [weight0, ~] = f;   
    for ii = 1:n
        Elem_new = Elem;
        Elem_new(ii,3) = Elem_new(ii,3) + delx;
        [weight_new, ~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_new);
        nfun = nfun + 1;
        grad_weight(ii,:) = (weight_new - weight0)/delx;
        %grad_stress(ii,:) = (stress_new - stress0)/delx;
    end
    
   
elseif j == 2
    for ii = 1:n
        Elem_left = Elem;
        Elem_left(ii,3) = Elem_left(ii,3) - delx;
        [weight_left, ~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_left);
       
        Elem_right = Elem;
        Elem_right(ii,3) = Elem_right(ii,3) + delx;
        [weight_right, ~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_right); 
        nfun = nfun + 1;
        grad_weight(ii,:) = (weight_right - weight_left)/(2*delx);
%         grad_stress(ii,:) = (stress_right - stress_left)/(2*delx);
    end

else
    for ii = 1:n
        Elem_new = Elem;
        Elem_new(ii,3) = Elem_new(ii,3) + 1i*delx;
        [weight_new, ~] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_new);
        nfun = nfun + 1;
        grad_weight(ii,:) = imag(weight_new)/delx;
%         grad_stress(ii,:) = imag(stress_new)/delx;
    end
end
end




