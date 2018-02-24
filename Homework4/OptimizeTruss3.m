
% ------------Starting point and bounds------------
%design variables
x0 = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]; %starting point (all areas = 5 in^2)
lb = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; %lower bound
ub = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20]; %upper bound
global nfun;
nfun = 0;
global j;
j = 2;
global delx;
delx = 1e-6;

% ------------Linear constraints------------
A = [];
b = [];
Aeq = [];
beq = [];

% ------------Call fmincon------------
tic;
options = optimoptions(@fmincon,'display','iter-detailed','Diagnostics','on');
options = optimoptions(options, 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',...
    true, 'FiniteDifferenceType', 'central', 'CheckGradients', true);
[xopt, fopt, exitflag, output] = fmincon(@fungrad, x0, A, b, Aeq, beq, lb, ub, @confungrad, options);
eltime = toc;
eltime
xopt    %design variables at the minimum
fopt    %objective function value at the minumum
[f, c, ceq] = objcon(xopt);
c
nfun

% ------------Objective and Non-linear Constraints------------
function [f, c, ceq] = objcon(x)
global nfun;

%get data for truss from Data.m file
Data

% insert areas (design variables) into correct matrix
for i=1:nelem
    Elem(i,3) = x(i);
end

% call Truss to get weight and stresses
[weight,stress] = Truss(ndof, nbc, nelem, E, dens, Node, force, bc, Elem);

%objective function
f = weight; %minimize weight

%inequality constraints (c<=0)
c = zeros(20,1);         % create column vector
for i=1:10
    c(i) = stress(i)/1000 - 25000/1000; % check stress both pos and neg
end
for i=11:20
    c(i) = -25000/1000 - stress(i-10)/1000;
end

%equality constraints (ceq=0)
ceq = [];
nfun = nfun + 1;


end

% ------------Separate obj/con (do not change)------------
function [f] = obj(x)
[f, ~, ~] = objcon(x);
end
function [c, ceq] = con(x)
[~, c, ceq] = objcon(x);
end


function [ f, grad_weight ] = fungrad( x )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
global j;
global delx;

n = length(x);
[f,~,~] = objcon(x);
grad_weight = zeros(n,1);
if j == 1
    weight0 = f;
    for ii = 1:n
        x_new = x;
        x_new(ii) = x(ii) + delx;
        [weight_new,~,~] = objcon(x_new);
        grad_weight(ii) = (weight_new - weight0)/delx;
    end
    
    
elseif j == 2
    for ii = 1:n
        x_left = x;
        x_left(ii) = x(ii) - delx;
        weight_left = obj(x_left);
        
        x_right = x;
        x_right(ii) = x_right(ii) + delx;
        weight_right = obj(x_right);
        grad_weight(ii) = (weight_right - weight_left)/(2*delx);
    end
    
else
    for ii = 1:n
        x_new = x;
        x_new(ii) = x_new(ii) + 1i*delx;
        weight_new = obj(x_new);
        grad_weight(ii) = imag(weight_new)/delx;
    end
end
end


function [ c, ceq, DC, DCeq ] = confungrad (x)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
global j;
global delx;

n = length(x);
[~,c,~] = objcon(x);

ceq = [];
DCeq = [];
DC=zeros(2*n,n);
if j == 1 % forward difference
    stress0 = c;
    for ii = 1:n
        x_new = x;
        x_new(ii) = x_new(ii) + delx;
        [~,stress_new,~] = objcon(x_new);       
        DC(:,ii) = (stress_new - stress0)/delx;   
    end
    DC = DC.';
    
elseif j == 2 % central difference
    for ii = 1:n
        x_left = x;
        x_left(ii) = x_left(ii) - delx;
        [~, stress_left, ~] = objcon(x_left);
        
        x_right = x;
        x_right(ii) = x_right(ii) + delx;
        [~, stress_right,~] = objcon(x_right);
        DC(:,ii) = (stress_right - stress_left)/(2*delx);
    end
    DC = DC.';
    
else % complex step 
    for ii = 1:n
        x_new = x;
        x_new(ii) = x_new(ii) + 1i*delx;
        [~, stress_new, ~] = objcon(x_new);
        DC(:,ii) = imag(stress_new)/delx;
    end
    DC = DC.';
end
end

