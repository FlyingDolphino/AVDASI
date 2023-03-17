




function[xopt] = optimize(shear)
    LoadData = xlsread('Load.xlsx'); %#ok<XLSRD> 

    %the optimizer

    x.tSkin = [0   0.03     
               1   0.005]; 

    x.tWeb  = [0   0.03    
               1   0.005]; 

    Stringer          = 4;         
    x.StringerHeight    = [0   0.05    
                           1   0.02]; 

    x.StringerThickness = [0   0.003 ;
                           1   0.002]; 



    results = zeros(16,8); 
    fvalCheck= 1e50;
    smallest = 0;
    for i = 4:20
        Stringer = i;
        

        x0 = [x.tSkin(1,2),x.tSkin(2,2),x.tWeb(1,2),x.tWeb(2,2),x.StringerHeight(1,2),x.StringerHeight(2,2),x.StringerThickness(1,2),x.StringerThickness(2,2)];
        lb = [1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 0.5e-3 0.5e-3];
        ub = [0.1 0.1 0.1 0.1 0.2 0.2 0.1 0.1];

        if shear == false
            nonlcon = @(x) nonlinear(x, LoadData,Stringer);
            options = optimoptions(@fmincon,'Display','iter','MaxFunctionEvaluations',4e10,'MaxIterations',5e10);
            [results(i,:), fval] = fmincon(@(x)optimiseThis(x,Stringer),x0,[],[],[],[],lb,ub,nonlcon,options);
            
            if fvalCheck > fval
                fvalCheck = fval;
                smallest = i;
            end
            
            
        end

        if shear == true
            nonlcon = @(x) fullopt(x, LoadData,Stringer);
            options = optimoptions(@fmincon,'Display','iter','MaxFunctionEvaluations',4e10,'MaxIterations',5e10);
            [results(i,:), fval] = fmincon(@(x)optimiseThis(x,Stringer),x0,[],[],[],[],lb,ub,nonlcon,options);
            
            if fvalCheck > fval
                fvalCheck = fval;
                smallest = i;
            end
        end
        
        
        
    end
    
    xopt = [results(smallest,:) (smallest) fvalCheck];
        
    
end

