% Created by João Meneses for ICNAAM Conference
% CDRSP - IPL (2019)

%--------------------------------------------------------------------------
% Simulation Geometry Constants & Coordinates (Arround Scaffold ROI)
%--------------------------------------------------------------------------
wl = 0.010; %Wire Length - Meters
wr = 0.0002; %Wire Radious - Meters
roi = 0.002; %ROI Square Size - Meters
z = (-roi/2:0.0001:roi/2); %Cilyndrical Coordinate Z - Meters 
r = (wr:0.0001:roi); %Cilyndrical Coordinate R - Meters
theta = (-1*pi:0.1:pi); %Cilyndrical Coordinate Theta - Radians

%--------------------------------------------------------------------------
% SINGLE WIRE  
%--------------------------------------------------------------------------
% >> Electrically insulated - (ex:PCL polymer)
% >> Finite Uniformly Charged Wire
% >> Constant Static Potencial Applyed to Wire (ex: 1 Coulomb)
% >> E Field with radial(r) and length(z) components
%--------------------------------------------------------------------------
lambda = 1; %Coulomb per Meter (charge linear density)
e0 = 8.8541878128e-12; %Permittivity of free space
eR = 2.20; %PCL Permittivity

rElectricFieldWirePotencial = zeros(size(z,2),size(r,2));
zElectricFieldWirePotencial = zeros(size(z,2),size(r,2));
totalVectorElements = size(rElectricFieldWirePotencial,1)*size(rElectricFieldWirePotencial,2);
rVector = zeros(1,totalVectorElements);
zVector = zeros(1,totalVectorElements);
rElectricVector = zeros(1,totalVectorElements);
zElectricVector = zeros(1,totalVectorElements);
i = 1;
for ri = 1:1:size(r,2)
    for zi = 1:1:size(z,2)
        % ROI Elements
        rVector(1,i) = r(1,ri);
        zVector(1,i) = z(1,zi);
        % Electric Field Calculation
        thetaL = atan(abs(-1*(wl/2)-z(1,zi))/r(1,ri));
        thetaR = atan(abs((wl/2)-z(1,zi))/r(1,ri));
        rElectricFieldWirePotencial(zi,ri) = ((eR*lambda)/r(1,ri))*(sin(thetaR)+sin(thetaL));
        zElectricFieldWirePotencial(zi,ri) = ((eR*lambda)/r(1,ri))*(cos(thetaR)+cos(thetaL));
        rElectricVector(1,i) = rElectricFieldWirePotencial(zi,ri);
        zElectricVector(1,i) = zElectricFieldWirePotencial(zi,ri);
        i = i + 1;
    end
end

figure
subplot(1,2,1);
scatter3(rVector,zVector,rElectricVector);
xlabel('(r) radial position - meter') 
ylabel('(z) axial position - meter') 
zlabel('E(r) Electric Field Intensity - volt/meter')

subplot(1,2,2);
scatter3(rVector,zVector,zElectricVector);
xlabel('(r) radial position - meter') 
ylabel('(z) axial position - meter') 
zlabel('E(z) Electric Field Intensity - volt/meter')