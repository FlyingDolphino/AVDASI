




function[xopt] = optimize()


    %the optimizer

    x.tSkin = [0   0.03     
               1   0.005]; 

    x.tWeb  = [0   0.03    
               1   0.005]; 

    x.Stringer          = 8;         
    x.StringerHeight    = [0   0.05    
                           1   0.02]; 

    x.StringerThickness = [0   0.003 ;
                           1   0.002]; 




    x0 = [x.tSkin(1,2),x.tSkin(2,2),x.tWeb(1,2),x.tWeb(2,2),x.Stringer,x.StringerHeight(1,2),x.StringerHeight(2,2),x.StringerThickness(1,2),x.StringerThickness(2,2)];
    lb = [0 0 0 0 2 0 0 0 0];
    ub = [0.1 0.1 0.1 0.1 16 0.2 0.2 0.1 0.1];

    LoadData = xlsread('Load.xlsx'); %#ok<XLSRD> 
    nonlcon = @(x) nonlinear(x, LoadData);
    options = optimoptions(@fmincon,'Display','iter','MaxFunctionEvaluations',4e10,'MaxIterations',5e10);
    [xopt, fval] = fmincon(@optimiseThis,x0,[],[],[],[],lb,ub,nonlcon,options);
    
end

