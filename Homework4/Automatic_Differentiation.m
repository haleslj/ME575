function [Func,Jacobian]=Automatic_Differentiation(x)
% design variables
        d = valder(x(1), [1,0,0,0]); %wire diameter (in)
        D = valder(x(2), [0,1,0,0]); %coil diameter (in)
        n = valder(x(3), [0,0,1,0]); %number of coils
        h_f = valder(x(4), [0,0,0,1]); %free height(in)
        
        %Other analysis variables
        h_0 = 1.0; %pre-load height (in)
        delta_0 = 0.4; %deflection from pre-loaded height (in)
        S_f = 1.5; %factor of safety (psi?)
        S_e = 45000; %endurance limit (psi)
        G = 12e6; %Sheer modulus of the material (psi)
        w = 0.18; %material property
        Q = 150000; %material property (psi)
        
        %Analysis Functions
        k = G*d^4/(8*D^3*n); %Spring stiffness (lb/in)
        K = (4*D - d)/(4*(D - d)) + 0.62*d/D; %Wahl factor 
        F = k*(h_f - h_0); %Force at preload height (lb)
        F_max = F + delta_0*k; %Force at deflected height (lb)
        tau_min = 8*F*D*K/(pi*d^3); %stress in a spring with axial load of F
        tau_max = 8*F_max*D*K/(pi*d^3); %stress in spring with max deflection
        tau_m = (tau_max + tau_min)/2; % mean shear stress
        tau_a = (tau_max - tau_min)/2; %alternating shear stress
        h_s = n*d; %solid height (in)
        S_y = 0.44*Q/d^w; %yield strength
        c_a = (h_0 - delta_0) - h_s; %clash allowance (in)
        F_solid = k*(h_f - h_s);
        tau_solid = 8*F_solid*D*K/(pi*d^3); %stress in a spring at solid height
        
        Func = (F.val);
        Jacobian = (F.der);
end