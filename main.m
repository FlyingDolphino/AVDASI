% Main starting code
clear; close all; clc; format short g; format compact;
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultAxesFontSize',20)
nSpan = 10;     % number of span cross-section to model = number of ribs (constant)
PLOT  = true;   % Boolean to turn plots on/off
optimiser = false; %turns the optimiser on and off


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
wingboxZs = zeros(nSpan,4);
stringerCenters = zeros(nSpan,8);
stringerZs = zeros(nSpan,8);
ymax = zeros(nSpan,1);
ymin = zeros(nSpan,1);
zmax = zeros(nSpan,1);
zmin = zeros(nSpan,1);
secondMoment = zeros(nSpan,1);



for n = 1:nSpan

    nodes = CS(n).WingBoxCornerXYZ;
    nodes(:,1) = [];        %only interested in the current box

    ymax(n) = max(nodes(:,1));          %max and min y values of each rib. Needed for max and min stresses later
    ymin(n) = min(nodes(:,1));
    zmax(n) = max(nodes(:,2));
    zmin(n) = min(nodes(:,2));
    tSkin = CS(n).tSkin;
    tWeb = CS(n).tWeb;
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

    Izz = 0;
    Iyy = 0;


    %summing the area and Izz, Iyy contributions of the box
    
    area = 2*max(L)*tSkin + 2*min(L)*tWeb;
    for i =1:ne
        if mod(i,2) == 0 
            Izz = Izz + y(i)^2.*tWeb*L(i);
            Iyy = Iyy + z(i)^2.*tWeb*L(i);
        else
            Iyy = Iyy + z(i)^2.*tSkin*L(i);
            Izz = Izz + y(i)^2.*tSkin*L(i);
        end
        
    end
    
    secondMoment(n)=Izz;
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
         yCS = yCS+ E*y(i)*L(i)*tSkin;
         zCS = zCS+ E*z(i)*L(i)*tWeb;
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
loadFactor = [-1,2.5]; %load factor
x = LoadData(:,1);
lenx = length(x);
aeroperL = LoadData(:,2);
massperL = LoadData(:,3);
Q = zeros(lenx,length(loadFactor));
BM = zeros(lenx,length(loadFactor));

%Q and BM distro
%wingbox is the contribution of the wingbox
wingbox= zeros(nSpan,2);

for i = 1:nSpan
    xCoord = CS(i).WingBoxCornerXYZ;
    wingbox(i,:) = [xCoord(1) LD(i)];
    
end

wingboxL = spline(1:10,wingbox(:,2),x);



for i = 1:lenx-1
    for j = 1:length(loadFactor)
        Q(i,j) = trapz(x(i:lenx),aeroperL(i:lenx)*loadFactor(j))-trapz(x(i:lenx),(massperL(i:lenx)+wingboxL(i:lenx))*g);
    end
end


for i = 1:lenx-1
    for j = 1:length(loadFactor)
        BM(i,j) = -trapz(x(i:lenx),Q(i:lenx,j));
    end
end


    


%% Step 3. Compute axial stresses caused by bending 
EIz = spline(1:nSpan,EIz,x);
ymax = spline(1:nSpan,ymax,x);
ymin = spline(1:nSpan,ymin,x);



sigma = E*(BM./EIz); %normalized to y.
sigmaMax = sigma.*ymax;
sigmaMin = sigma.*ymin;

yplot = linspace(ymin(1),ymax(1));
sigmaVariation = yplot.*sigma(1);


%%transverse stresses caused by bending
EIy = spline(1:nSpan,EIy,x);
zmax = spline(1:nSpan,zmax,x);
zmin = spline(1:nSpan,zmin,x);


sigmayy = E*(BM./EIy);
sigmayyMax = sigmayy.*zmax;
sigmayyMin = sigmayy.*zmin;



%% Step 4. Compute failure caused by axial stresses
% 4.1  Yield (element per element based on max sigma_xx)

gyield = abs(sigmaMax(:,2))-276e6;


% 4.2  Top and Bottom Skin Buckling
L = 1;
spCrit = zeros(nSpan,1);
Iplates = zeros(nSpan,1);
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
    I = (Ibox + abox*(boxCenter-centroid)^2)+4*(Istring+astring*(centroid-stringerCenter)^2);
    Iplates(n) = I;
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

maxsigma = spline(x,sigmaMax(:,2),(1:nSpan).');
gbuckling = [abs(maxsigma)-spCrit, abs(maxsigma)-plateCrit, abs(maxsigma)-stringerCrit];



%% Step 5. Finite element model and wing deflections

%Will use the x coords given by the loadData as our nodes
ne = lenx-1;    %number of elements
nodes = x;
Le = nodes(2) - nodes(1);
eConn = zeros(ne,2);

Ke = zeros(4,4,ne);

for i = 1:ne
    eConn(i,:) = [i i+1];
end


for i=1:ne
    Ke(:,:,i) = EIz(i)/(Le^3)*[...
       12   6*Le    -12     6*Le 
       6*Le     4*Le^2  -6*Le   2*Le^2
       12  6*Le    -12  6*Le
       -6^Le     -2*Le^2  6*Le   4*Le^2]  ;
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


%distro load needs to change

distroLoad = (aeroperL*loadFactor)-(massperL*g)-(wingboxL*g);


pointLoads = zeros(lenx,2);
for i = 1:2:lenx
    for j = 1:2
        pointLoads(i,j) = 0.5*(distroLoad(i,j)+distroLoad(i+1,j));
        pointLoads(i+1,j) = pointLoads(i,j);
    end
end

uY = zeros(50,2);
ind = 1;
for j = 1:2
    F(2:2:100) = pointLoads(:,j);
    idBC = [1,2];
    idFree = (3:100);

    Kf          = K(idFree,idFree);
    Ff          = F(idFree);
    u(idFree)   = Kf^-1*Ff;          % displacement solution caused by F
    displ = reshape(u,2,[])';
    uY(:,j) = displ(:,2);
end


%% Step 6. Design the lightest wing that does not fail due to axial stresses, 
% and remains below a given maximum deflection constraint

%see after plotting


    


%% Step 7. Compute shear stresses caused by bending

Q1 = spline(x,Q(:,1),1:nSpan);
Q2 = spline(x,Q(:,2),1:nSpan);
V = [Q1;Q2]';
qstress = zeros(nSpan,298);
arclength = zeros(nSpan,298);
maxShear = zeros(nSpan,2);
minShear = zeros(nSpan,2);

for j = 1:2
    for n = 1:nSpan
        
        tSkin = CS(n).tSkin;
        tWeb = CS(n).tWeb;
        
        z_length = linspace(0,zmax(n));
        y_length = linspace(0,ymax(n));
        

        Qtop = tSkin*ymax(n)*z_length;
        qtop = V(n,j).*Qtop./(secondMoment(n)*tSkin);
        
        Qside = tWeb*y_length.*(ymax(n)-y_length);
        qside = qtop(end) + V(n,j).*Qside./(secondMoment(n)*tWeb);
        
        %join together
        qstress(n,:) = [qtop(1:end-1) qside(1:end-1) fliplr(qtop)];

        arc = [z_length(1:end-1) y_length+z_length(end)];
        arclength(n,:) = [arc(1:end-1) z_length+arc(end)];
        
        
        maxShear(n,j) = max(abs(qstress(n,:)));
        minShear(n,j) = -maxShear(n,j);
   
    end
end

%% Step 8. Compute failure caused by combined shear and axial stresses

%principal stresses
tau_max = zeros(length(x),2);
tau_max(:,1) = spline(1:nSpan,maxShear(:,1),x);
tau_max(:,2)  = spline(1:nSpan,maxShear(:,2),x);

sigma1 = (sigmaMax + sigmaMin) / 2 + sqrt(((sigmaMax - sigmaMin) / 2).^2 + tau_max.^2);
sigma2 = (sigmaMax + sigmaMin) / 2 - sqrt(((sigmaMax - sigmaMin) / 2).^2 + tau_max.^2);

vonMises = sqrt(0.5*((sigmaMax - sigmaMin).^2 + 3*tau_max.^2));
vonMises2 = sqrt(sigma1.^2 - sigma1.*sigma2 + sigma2.^2 + 3*tau_max.^2);


vonFail = vonMises2 - YieldStrength;


%% Step 9. Design the lightest wing able to resist failure due to axial and shear stresses,  
% and remains below a given maximum deflection constraints
%see after plotting


%%plotting

if PLOT==true
    %Example of discretised section
    figure(2);
    subplot(2,3,6);
    hold on;
    scatter(CS(1).WingBoxCornerXYZ(:,3),CS(1).WingBoxCornerXYZ(:,2),'filled')   %corners
    line(CS(1).WingBoxCornerXYZ(:,3),CS(1).WingBoxCornerXYZ(:,2));
    scatter(wingboxZs(1,:),wingboxCenters(1,:),'filled');                       %midpoints of elements
    scatter(stringerZs(1,:),stringerCenters(1,:),'filled');
    line([CS(1).WingBoxCornerXYZ(4,3),CS(1).WingBoxCornerXYZ(1,3)],[CS(1).WingBoxCornerXYZ(4,2),CS(1).WingBoxCornerXYZ(1,2)]);
    for i = 1:length(CS(1).TopStringerXYZ)
        line(CS(1).TopStringerXYZ{i}(:,3),CS(1).TopStringerXYZ{i}(:,2))
        line(CS(1).BotStringerXYZ{i}(:,3),CS(1).BotStringerXYZ{i}(:,2))
    end
    legend('Wingbox Element Ends','Elements','Wingbox Element Midpoints','Stringer Element Midpoints','Location','best');
    ylabel('Y position');
    xlabel('X position');
    grid on;
    title('Example of descretized section');
    hold off;
    
    %variation of properties
    
    EA = spline(1:nSpan,EA,x);
    GJ = spline(1:nSpan,GJ,x);
    LD = spline(1:nSpan,LD,x);
    

    hold on;
    %EIzz plot
    subplot(2,3,1);
    plot(x,EIz);
    grid on;
    xlabel('Span location (m)');
    ylabel('EIzz (Nm^4)');
    title('Bending Stiffness about z-z');
    
    
    %EIyy plot
    subplot(2,3,2);
    plot(x,EIy);
    grid on;
    xlabel('Span location (m)');
    ylabel('EIyy (Nm^4)');
    title('Bending Stiffness about y-y')
    
    %EA plot
    subplot(2,3,3);
    plot(x,EA);
    grid on;
    xlabel('Span Location (m)');
    ylabel('EA (N)');
    title('Axial Stiffness');
    
    %GJ plot
    subplot(2,3,4);
    plot(x,GJ);
    grid on;
    xlabel('Span Location');
    ylabel('GJ (Nm^2)')
    title('Torsional Stiffness');
    
    %LD plot
    subplot(2,3,5);
    plot(x,LD);
    grid on;
    xlabel('Span Location (m)');
    ylabel('Linear Density (kg/m)');
    title('Linear Density');
    f = figure;

    %Question 2 plot
    figure(3);
    %shear force plot
    hold on;
    plot(x,Q(:,1));
    plot(x,Q(:,2));
    grid on;
    xlabel('Span Location (m)');
    ylabel('Internal Shear Force (N)');
    legend('-1g','2.5g');


    figure(4)
    hold on;
    plot(x,BM(:,1));
    plot(x,BM(:,2));
    grid on;
    xlabel('Span Location (m)');
    ylabel('Internal Bending Moment (Nm)');
    legend('-1g','2.5g');

    
    %question 3 plot
    figure(5)
    plot(sigmaVariation,yplot);
    grid on;
    xline(0)
    xlabel('Stress (Pa)');
    ylabel('Y coordinate in cross section');
    title('Stress Variation with Y, 1st Rib -1g Loading');
    
    figure(6)
    hold on;
    grid on;
    plot(x,sigmaMax(:,1),'k');
    plot(x,sigmaMax(:,2),'r');
    plot(x,sigmaMin(:,1),'k');
    plot(x,sigmaMin(:,2),'r');
    legend('-1g load','2.5g load');
    title('Maximum and Minimum Axial Stress Along The Span')
    xlabel('Span (m)')
    ylabel('Stress (Pa)');
    
    
    
    %question 4 plots
    spCrit = gbuckling(:,1);
    plateCrit = gbuckling(:,2);
    stringCrit = gbuckling(:,3);
    gyield = spline(x,gyield,1:nSpan);
    
    figure(7)
    hold on
    grid on
    plot(1:nSpan,spCrit,'k','LineStyle','--','LineWidth',1.5)
    plot(1:nSpan,plateCrit,'k','LineStyle',':','LineWidth',1.5)
    plot(1:nSpan,stringCrit,'k','LineStyle','-.','LineWidth',1.5)
    plot(1:nSpan,gyield,'k','LineWidth',1.5)
    ylabel('Failiure Indices (Pa)')
    xlabel('Span location (m)')
    ylim([-0.5e9 0.25e9])
    yline(0,'r')
    legend('Top and Bottom Column Buckling','Plate Buckling Between Stringers','Stringers Buckling','Yield');
    title('Failiure Indices, Positive Means Failiure')
    
    
   
    
    %Question 5 plotting
    %plot of the forces and the point load equivalents
    figure(8)
    hold on;
    plot(x,distroLoad);
    stem(x,pointLoads);
    yline(0);
    grid on;
    xlabel('Spanwise Location (m)');
    ylabel('Resultant Force (N)');
    title('Plot of the distributed load and the point load equivalent')
    legend('-1g','2.5g');
    
    figure(9);
    hold on;
    yyaxis left
    plot(x,uY(:,1),'r');
    plot(x,uY(:,2),'b','LineStyle','-');
    ylabel('Deflection(m)');
    ylim([-12e-4,12e-4]);
    
    
    
    yyaxis right
    plot(x,gradient(uY(:,1)),'r','LineStyle','--');
    plot(x,gradient(uY(:,2)),'b','LineStyle','--');
    
    xlabel('Spanwise Location(m)')
    ylabel('Deflection rate');
    ylim([-1.5e-4,1.5e-4]);
    yline(0);
    legend('-1g Load Factor', '2.5g Load factor','Location','NorthWest')
    title('Wing deflection (Solid lines) and Gradient (Dashed lines)')
    
    
    figure(10)
    plot(x,vonFail);
    yline(0);
    title('Von Mises Failure Indices for Both Load Cases')
    ylabel('Von Mises Failure Indices (Pa)')
    xlabel('Span location (m)')
    legend('-1g Load', '2.5g Load')
    
    
    
    figure(11)
    plot(arclength(1,:),qstress(1,:))
    title('Shear Stress Variation of top and bottom surface of 1st rib');
    xlabel('Arc length, increasing clockwise (m)')
    ylabel('Magnitude of Shear Stress (Pa)')
    xline(zmax(1));
    xline(zmax(1)+ymax(1));
    ylim([0 18e6]);
    xlim([0 2.06]);
    txt = 'Half of Top surface';
    text(0.6,17e6,txt,'HorizontalAlignment','right','FontSize',20)
    txt = 'Right Surface';
    text(0.95,17e6,txt,'FontSize',20);
    txt = 'Half of Bottom Surface';
    text(1.8,17e6,txt,'HorizontalAlignment','right','FontSize',20)
    
    
    figure(12)
    
    hold on;
    grid on;
    plot(1:nSpan,maxShear(:,1),'k');
    plot(1:nSpan,maxShear(:,2),'r');
    plot(1:nSpan,minShear(:,1),'k');
    plot(1:nSpan,minShear(:,2),'r');
    legend('-1g load','2.5g load');
    title('Maximum and Miminum Shear Stress Along The Span')
    xlabel('Span (m)')
    ylabel('Shear Stress (Pa)');
    
    figure(13)
    vonPlot = spline(x,vonFail(:,2),1:nSpan);
    
    
    hold on
    grid on
    plot(1:nSpan,spCrit,'k','LineStyle','--','LineWidth',1.5)
    plot(1:nSpan,plateCrit,'k','LineStyle',':','LineWidth',1.5)
    plot(1:nSpan,stringCrit,'k','LineStyle','-.','LineWidth',1.5)
    plot(1:nSpan,gyield,'k','LineWidth',1.5)
    plot(1:nSpan,vonPlot,'b','LineWidth',1.5);
    ylabel('Failiure Indices (Pa)')
    xlabel('Span location (m)')
    ylim([-0.5e9 0.25e9])
    yline(0,'r')
    legend('Top and Bottom Column Buckling','Plate Buckling Between Stringers','Stringers Buckling','Yield','Von Mises');
    title('Failiure Indices, Positive Means Failiure')

    

end

%optimising
if optimiser == true
    xopt = optimize(false);
    xoptWithShear = optimize(true);
end


