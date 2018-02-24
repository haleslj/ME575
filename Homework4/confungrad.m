function [ c, ceq, DC, DCeq ] = confungrad (x)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
global j;
global nfun;

delx = 1e-4;
n = length(x);

Data;
Elem(:,3) = x;
[~,c] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);
nfun = nfun + 1; 
ceq = [];
DCeq = [];

if j == 1 % forward difference
    stress0 = c;
    for ii = 1:n
        Elem_new = Elem;
        Elem_new(ii,3) = Elem_new(ii,3) + delx;
        [~, stress_new] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_new);
        nfun = nfun +1; 
%         grad_weight(ii,:) = (weight_new - weight0)/delx;
        DC(ii,:) = (stress_new - stress0)/delx;

    end
    
elseif j == 2 % central difference
        for ii = 1:n
        Elem_left = Elem;
        Elem_left(ii,3) = Elem_left(ii,3) - delx;
        [~, stress_left] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_left);
       
        Elem_right = Elem;
        Elem_right(ii,3) = Elem_right(ii,3) + delx;
        [~, stress_right] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_right); 
        nfun = nfun + 1;
%         grad_weight(ii,:) = (weight_right - weight_left)/(2*delx);
        DC(ii,:) = (stress_right - stress_left)/(2*delx);
    end
        
else % complex step
            
    for ii = 1:n
        Elem_new = Elem;
        Elem_new(ii,3) = Elem_new(ii,3) + 1i*delx;
        [~, stress_new] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem_new);
        %grad_weight(ii,:) = imag(weight_new)/delx;
        DC(ii,:) = imag(stress_new)/delx;
    end
end
DC = cat(1,DC, -1*DC);

end

