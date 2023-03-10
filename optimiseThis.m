function [wingmass] = optimiseThis(x0)

    
    %% Setup 1: Material definitions for Aluminium
    E     = 70e9;              % Young's modulus [Pa] 
    nu    = 0.34;              % Poisson's ratio
    G     = E/(2*(1+nu));       % shear modulus [Pa]
    rho   = 2700;              % density [kg/m^3]
    nSpan = 10;

    %% Setup 2. Design variables and geometry (to be given by the students)
    % BoxGeo: wing box geometry               [SpanLocation [normalised], StartChord [normalised], EndChord [normalised]]
    % tSkin : wingbox skin/spar thickness     [SpanLocation [normalised], sparcap skin thickness [m]]
    % tWeb  : wingbox web thickness           [SpanLocation [normalised], web thickness[m]]
    % Stringer : integer number of flat stringer, uniformly spaced on top and bottom
    % StringerHeight: height of stringers     [SpanLocation [normalised], height [m]]  
    % StringerThickness: stringers thickness  [SpanLocation [normalised], thickness [m]] 

    x.BoxGeo = [0  0.2   0.5;
                0.5 0.2  0.5;
                1  0.2   0.35]; 

    x.tSkin = [0   x0(1)     
               1   x0(2)]; 

    x.tWeb  = [0   x0(3)    
               1   x0(4)]; 
    x.Stringer          = round(x0(5));  
    if mod(x.Stringer,2) ~=0
        x.Stringer = x.Stringer +1;
    end         
    x.StringerHeight    = [0   x0(6)    
                           1   x0(7)]; 

    x.StringerThickness = [0   x0(8) ;
                           1   x0(9)]; 

    % CS is an array of 10 structures containing the geometric information
    % relevant to each of the 10 cross-sections along the span. To view, say
    % the geometry of the 2nd cross-section you can type CS(2). Then to
    % access, say the web thickness of the 2nd cross-section, you can type 
    % CS(2).tWeb
    CS = WingParameterisation(x,nSpan,false); 

    %% Step 1. Discretised cross-sections and approximate cross-section properties 

    EA = zeros(nSpan,1);
    EIz = zeros(nSpan,1);
    EIy = zeros(nSpan,1);
    ye = zeros(nSpan,1);
    ze = zeros(nSpan,1);
    GJ = zeros(nSpan,1);
    LD = zeros(nSpan,1);
    wingboxCenters = zeros(nSpan,4);
    wingboxZs = zeros(nSpan,4);
    stringerCenters = zeros(nSpan,8);
    stringerZs = zeros(nSpan,8);
    ymax = zeros(nSpan,1);
    ymin = zeros(nSpan,1);


    for n = 1:nSpan

        nodes = CS(n).WingBoxCornerXYZ;
        nodes(:,1) = [];        %only interested in the current box

        ymax(n) = max(nodes(:,1));          %max and min y values of each rib. Needed for max and min stresses later
        ymin(n) = min(nodes(:,1));
        t = CS(n).tSkin;
        ne = length(nodes);
        L = zeros(1,ne);
        y = zeros(1,ne);
        z = zeros(1,ne);
        tots = length(CS(n).TopStringerXYZ)+length(CS(n).BotStringerXYZ);
        ys = zeros(1,tots);
        zs = zeros(1,tots);

        %length of each element in box
        for i = 1:ne
            if i ==ne
                L(i) = norm(nodes(1,:)-nodes(i,:));
            else
                L(i) = norm(nodes(i+1,:)-nodes(i,:));
            end
        end

        %midpoints of each box element
        for i = 1:ne
            if i ==ne
                y(i) = 1/2*(nodes(1,1)+nodes(i,1));
                z(i) = 1/2*(nodes(1,2)+nodes(i,2));
                wingboxCenters(n,i) = y(i);
                wingboxZs(n,i) = z(i);
            else
                y(i) = 1/2*(nodes(i+1,1)+nodes(i,1));
                z(i) = 1/2*(nodes(i+1,2)+nodes(i,2));
                wingboxCenters(n,i) = y(i);
                wingboxZs(n,i) = z(i);
            end
        end

        area = 0;
        Izz = 0;
        Iyy = 0;

        %summing the area and Izz, Iyy contributions of the box
        for i =1:ne                             
            area = area + t*L(i);
            Izz = Izz + y(i)^2.*t*L(i);
            Iyy = Iyy + z(i)^2.*t*L(i);
        end

        %top stringer handling
        stringerNodes = CS(n).TopStringerXYZ;
        ts = CS(n).StringerThickness;
        ns = length(stringerNodes);
        Ls = zeros(1,ns);



        for i= 1:ns
            a = stringerNodes{i}(1,:);
            b = stringerNodes{i}(2,:);
            Ls(i) = norm(a-b); %length of each stringer
            mid = (a+b)/2;      %midpoint of each stringer
            ys(i) = mid(2);     %y coordinate of midpoint
            zs(i)= mid(3);      %z coordinate of midpoint
            stringerCenters(n,i) = ys(i);   %saves the y midpoints of each stringer, needed for later
            stringerZs(n,i) = zs(i);
        end

        %bottom stringer handling
        stringerNodes = CS(n).BotStringerXYZ;
        for i= (ns)+1:2*ns
            a = stringerNodes{i-ns}(1,:);
            b = stringerNodes{i-ns}(2,:);
            Ls(i) = norm(a-b);
            mid = (a+b)/2;
            ys(i) = mid(2);
            zs(i)= mid(3);
            stringerCenters(n,i) = ys(i);
            stringerZs(n,i) = zs(i);
        end

        Sarea = 0;
        SIzz = 0;
        SIyy = 0;

        for i =1:tots
            Sarea = Sarea + ts*Ls(i);
            SIzz = SIzz + ys(i)^2.*ts*Ls(i);
            SIyy = SIyy + zs(i)^2.*ts*Ls(i);
        end


        area = area + Sarea;
        Izz = Izz + SIzz;
        Iyy = Iyy + SIyy;
        J = Izz + Iyy;


         EIy(n) = E*Iyy;
         EIz(n) = E*Izz;
         EA(n) = E*area;
         GJ(n) = G*J;
         LD(n) = area*rho;
         

         %centroid and elastic center
         yCS = 0;
         zCS = 0;
         for i = 1:ne
             yCS = yCS+ E*y(i)*L(i)*t;
             zCS = zCS+ E*z(i)*L(i)*t;
         end
         for i = 1:tots
             yCS = yCS + E*ys(i)*Ls(i)*ts;
             zCS = zCS + E*zs(i)*Ls(i)*ts;
         end



         %saving the elastic center (translated back to original coordinates)
         ye(n) = (1/EA(n))*yCS;
         ze(n) = (1/EA(n))*zCS;



    end
    
    wingmass = trapz(1:nSpan,LD);
    
