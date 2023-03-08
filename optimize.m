%the optimizer
clear all;


x.tSkin = [0   0.01     
           1   0.005]; 

x.tWeb  = [0   0.01    
           1   0.005]; 

x.Stringer          = 4;         
x.StringerHeight    = [0   0.05    
                       1   0.02]; 

x.StringerThickness = [0   0.003 ;
                       1   0.002]; 
       
                
                   
                   
x0 = [x.tSkin(1,2),x.tSkin(2,2),x.tWeb(1,2),x.tWeb(2,2),x.Stringer,x.StringerHeight(1,2),x.StringerHeight(2,2),x.StringerThickness(1,2),x.StringerThickness(2,2)];
A = [];
B = [];
lb = [0 0 0 0 0 0 0 0 0];
ub = [0.1 0.1 0.1 0.1 16 0.2 0.2 0.1 0.1];



options = optimoptions(@fmincon,'Display','iter');
[xopt, fval] = fmincon(@optimiseThis,x0,A,B,[],[],lb,ub,@nonlinear,options);



% Define the output function for fmincon
% function stop = outfun(x,optimValues,state)
%     stop = false;
%     if isequal(state,'iter')
%         fprintf('Iteration: %d, Function value: %f, Norm of gradient: %f\n', ...
%             optimValues.iteration, optimValues.fval, norm(optimValues.gradient));
%     end
% end
