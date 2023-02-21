% Main starting code
clear; close all; clc; format short g; format compact;
set(0,'DefaultFigureWindowStyle','docked')

nSpan = 10;     % number of span cross-section to model = number of ribs (constant)
PLOT  = false;   % Boolean to turn plots on/off


%% Setup 1: Material definitions for Aluminium
E     = 70e9;              % Young's modulus [Pa] 
nu    = 0.34;              % Poisson's ratio
G     = E/(2*(1+nu));       % shear modulus [Pa]
rho   = 2700;              % density [kg/m^3]
YieldStrength = 276*1e6;   % [Pa]
ShearStrength = 207*1e6;   % [Pa]

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

x.tSkin = [0   0.01     
           1   0.005]; 

x.tWeb  = [0   0.01    
           1   0.005]; 

x.Stringer          = 4;         
x.StringerHeight    = [0   0.05    
                       1   0.02]; 

x.StringerThickness = [0   0.003 ;
                       1   0.002]; 

% CS is an array of 10 structures containing the geometric information
% relevant to each of the 10 cross-sections along the span. To view, say
% the geometry of the 2nd cross-section you can type CS(2). Then to
% access, say the web thickness of the 2nd cross-section, you can type 
% CS(2).tWeb
CS = WingParameterisation(x,nSpan,PLOT); 




%% Step 1. Discretised cross-sections and approximate cross-section properties 

EA = zeros(nSpan,1);
EIz = zeros(nSpan,1);
EIy = zeros(nSpan,1);
ye = zeros(nSpan,1);
ze = zeros(nSpan,1);
GJ = zeros(nSpan,1);
LD = zeros(nSpan,1);
wingboxCenters = zeros(nSpan,4);
stringerCenters = zeros(nSpan,8);
ymax = zeros(nSpan,1);
ymin = zeros(nSpan,1);


%all inertia is calculated about 0,0 of the parameterizaition axes
for n = 1:nSpan

    nodes = CS(n).WingBoxCornerXYZ;
    nodes(:,1) = [];        %only interested in the current box
    
    
    ymax(n) = max(nodes(:,1));
    ymin(n) = min(nodes(:,1));
    t = CS(n).tSkin;
    ne = length(nodes);
    L = zeros(1,ne);
    y = zeros(1,ne);
    z = zeros(1,ne);
    tots = length(CS(n).TopStringerXYZ)+length(CS(n).BotStringerXYZ);
    ys = zeros(1,tots);
    zs = zeros(1,tots);
    %length of each element
    
    for i = 1:ne
        if i ==ne
            L(i) = norm(nodes(1,:)-nodes(i,:));
        else
            L(i) = norm(nodes(i+1,:)-nodes(i,:));
        end
    end
    
    %midpoints
    for i = 1:ne
        if i ==ne
            y(i) = 1/2*(nodes(1,1)+nodes(i,1));
            z(i) = 1/2*(nodes(1,2)+nodes(i,2));
            wingboxCenters(n,i) = y(i);
        else
            y(i) = 1/2*(nodes(i+1,1)+nodes(i,1));
            z(i) = 1/2*(nodes(i+1,2)+nodes(i,2));
            wingboxCenters(n,i) = y(i);
        end
    end

    area = 0;
    Izz = 0;
    Iyy = 0;

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
        Ls(i) = norm(a-b);
        mid = (a+b)/2;
        ys(i) = mid(2);
        zs(i)=mid(3); 
        stringerCenters(n,i) = ys(i);
    end
    
    %bottom stringer handling
    stringerNodes = CS(n).BotStringerXYZ;
    nsB = length(stringerNodes);
    for i= 5:ns+nsB
        a = stringerNodes{i-4}(1,:);
        b = stringerNodes{i-4}(2,:);
        Ls(i) = norm(a-b);
        mid = (a+b)/2;
        ys(i) = mid(2);
        zs(i)=mid(3);
        stringerCenters(n,i) = ys(i);
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
    
    
%inertias were taken about P4 for each wing box
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

%% Step 2. Wing loads, internal shear and bending moments
LoadData = xlsread('Load.xlsx'); %#ok<XLSRD> 


g = 9.81;
loadFactor = 2.5; %load factor
x = LoadData(:,1);
lenx = length(x);
aeroperL = LoadData(:,2)*loadFactor;
massperL = LoadData(:,3);
Q = zeros(lenx,1);
BM = zeros(lenx,1);

%Q and BM distro
%wingbox is the contribution of the wingbox
wingbox= zeros(nSpan,2);

for i = 1:nSpan
    xCoord = CS(i).WingBoxCornerXYZ;
    wingbox(i,:) = [xCoord(1) LD(i)];
    
end

wingboxL = spline(1:10,wingbox(:,2),x);

for i = 1:lenx-1
    Q(i) = trapz(x(i:lenx),aeroperL(i:lenx))-trapz(x(i:lenx),(massperL(i:lenx)+wingboxL(i:lenx))*g);
end


for i = 1:lenx-1
    BM(i) = trapz(x(i:lenx),Q(i:lenx));
end


    


%% Step 3. Compute axial stresses caused by bending 
sigma = zeros(nSpan,1);
EIz = spline(1:nSpan,EIz,x);
ymax = spline(1:nSpan,ymax,x);
ymin = spline(1:nSpan,ymin,x);



sigma = E*(BM./EIz);
sigmaMax = sigma.*ymax;
sigmaMin = sigma.*ymin;

yplot = linspace(ymin(1),ymax(1));
sigmaVariation = yplot.*sigma(1);






%% Step 4. Compute failure caused by axial stresses
% 4.1  Yield (element per element based on max sigma_xx)

gyield = sigmaMax-276e6;

% 4.2  Top and Bottom Skin Buckling
L = 1;
spCrit = zeros(nSpan,1);
for n = 1:nSpan
    nodes = CS(n).WingBoxCornerXYZ;
    b = nodes(2,3)-nodes(1,3);
    h = CS(n).tSkin;
    abox  = b*h;
    Ibox = b*(h^3)/12;
    boxCenter = nodes(1,2);
    
    %stringer contribution
    
    nodes = CS(n).TopStringerXYZ{1};
    b = CS(n).StringerThickness;
    h = nodes(1,2)-nodes(2,2);
    astring = b*h;
    Istring = b*(h^3)/12;
    stringerCenter = stringerCenters(n,1);
    
    %finding centroid of panel
    centroid = ((boxCenter*abox)+4*(stringerCenter*astring))/(abox+astring*4);
    
    I = (Ibox + abox*(boxCenter-centroid)^2)+4*(Istring+astring*(centroid-stringerCenter));
    r = (I/(abox+astring))^0.5;
    spCrit(n) = pi^2*E/(L/r)^2;
end


% 4.3  Plate Buckling in between Stiffeners

plateCrit = zeros(nSpan,1);
for n = 1:nSpan
    tsk = CS(n).tSkin;
    bsk = CS(n).TopStringerXYZ{2}(1,3)-CS(n).TopStringerXYZ{1}(1,3);
    plateCrit(n) = (4*pi^2*E)/(12*(1-nu^2)) * (tsk/bsk)^2;

end

% 4.4  Stiffeners buckling as plates

stringerCrit = zeros(nSpan,1);
for n = 1:nSpan
    tst = CS(n).StringerThickness;
    bst = CS(n).TopStringerXYZ{1}(1,2)-CS(n).TopStringerXYZ{1}(2,2);
    stringerCrit(n) = (0.43*pi^2*E)/(12*(1-nu^2)) * (tst/bst)^2;
end

%check step
spCrit = spline(1:nSpan,spCrit,x);
plateCrit = spline(1:nSpan,plateCrit,x);
stringerCrit = spline(1:nSpan,stringerCrit,x);

gbuckling = [sigmaMax-spCrit, sigmaMax-plateCrit, sigmaMax-stringerCrit];



%% Step 5. Finite element model and wing deflections

%Will use the x coords given by the loadData as our nodes
ne = lenx-1;    %number of elements
nodes = x;
Le = nodes(1)-nodes(2);
eConn = zeros(ne,2);

Ke = zeros(4,4,ne);

for i = 1:ne
    eConn(i,:) = [i i+1];
end


for i=1:ne
    Ke(:,:,i) = EIz(i)/(Le^3)*[...
       12   6*Le    -12     6*Le 
       6*Le     4*Le^2  -6*Le   2*Le^2
       -12  6*Le    12  -6*Le
       6^Le     2*Le^2  -6*Le   4*Le^2]  ;
end

K = zeros(100,100); %50 nodes * 2DoFs

for i=1:ne
    nodeID     = eConn(i,:);                  % index of nodes linked to element i
    id         = sort([nodeID*2-1 nodeID*2]); % index of dofs  linked to element i
    K(id,id)   = K(id,id) + Ke(:,:,i);        % add stiffnes of element i to global stiffness matrix
end

%Forces and boundary conditions
F = zeros(lenx*2,1);
u = zeros(lenx*2,1);

forcePerL = aeroperL - massperL;
xInd = 1:2:100;
   
for i = 1:2:ne
    
    elementForce = forcePerL(i) + forcePerL(i+1);
    forcePerNode = elementForce/2;
    F(xInd(i))=forcePerNode;
    F(xInd(i+1))=forcePerNode;
end


idBC = [1,2];
idFree = [3:100];

Kf          = K(idFree,idFree);
Ff          = F(idFree);
u(idFree)   = Kf^-1*Ff;          % displacement solution caused by F




%% Step 6. Design the lightest wing that does not fail due to axial stresses, 
% and remains below a given maximum deflection constraint





%% Step 7. Compute shear stresses caused by bending



%% Step 8. Compute failure caused by combined shear and axial stresses


%% Step 9. Design the lightest wing able to resist failure due to axial and shear stresses,  
% and remains below a given maximum deflection constraints





