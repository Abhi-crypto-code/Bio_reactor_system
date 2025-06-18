% <---------------------------------------------------------------------------->
% Function Num_Jacobian.m
% This function computes gradient /Jacobian matrix for any function vector
% using central difference method
% Input arguments to the function are
% fun_name  : String containing name of the MATLAB function,
%                     which returns function vector F(z) given z 
% z0             : vector at which gradient should be computed 
% Output arguments of the function are
% gradF        : (n x n) gradient matrix  
% ------------------------------------------------------------------------------

function gradF = Num_Jacobian( fun_name, Z0 ) 

nZ0 = length( Z0 ) ;   % Find dimension of vector Z0

for i = 1 : nZ0 
    
    % Generate a small perturbation in i'th variable     
    eps = Z0(i) / 10000 ;   
   
    if ( eps == 0 )
        eps = 1 / 10000 ;
    end     
    
    Zp = Z0 ;
    
    Zp(i) = Zp(i) + eps ;                % Introduce +ve perturbation in i'th varoable 
    
    fp = feval( fun_name, Zp ) ;   % Evaluate Function vector for perturbed vector  
    


    Zn = Z0 ;
    
    Zn(i) = Zn(i) - eps ;                 % Introduce +ve perturbation in i'th varoable
    
    fn = feval( fun_name, Zn ) ;   % Evaluate Function vector for perturbed vector  
    
    
    % Evaluate i'th column of the Jacobian matrix using central difference method
    
    gradF(:,i) = (fp - fn) / (2 * eps) ;
    
end     
    