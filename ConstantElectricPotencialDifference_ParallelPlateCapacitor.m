% Created by João Meneses for ICNAAM Conference
% CDRSP - IPL (2019)

%--------------------------------------------------------------------------
% Simulation Geometry Constants & Coordinates (Arround Scaffold ROI)
%--------------------------------------------------------------------------
wl = 0.010; %Plate Length - Meters
wr = 0.0002; %Plate Radious - Meters
roi = 0.002; %ROI Square Size - Meters
z = (-roi/2:0.0001:roi/2); %Cilyndrical Coordinate Z - Meters 
r = (wr:0.0001:roi); %Cilyndrical Coordinate R - Meters
theta = (-1*pi:0.1:pi); %Cilyndrical Coordinate Theta - Radians

%--------------------------------------------------------------------------
% PARALLEL PLATE CAPACITOR  
%--------------------------------------------------------------------------
% >> High dielectric medium - (ex:PCL polymer)
% >> Constant Static Potencial Difference Applyed to terminals (ex:100mV)
% >> E(r) Field with radial(r) component
% >> Plate geometry strong dependence
%--------------------------------------------------------------------------
potencialDifferenceCap = 0.1; %Volt 
relativePermitivity = 2.20; %Dielectric Constant PCL Polymer (article based doi:10.1109/TBME.2011.2157918)
parallelPlateCapDistance = 0.002; %Meter

electricFieldCap = zeros(size(z,2),size(r,2));
totalVectorElements = size(electricFieldCap,1)*size(electricFieldCap,2);
rVector = zeros(1,totalVectorElements);
zVector = zeros(1,totalVectorElements);
rElectricVector = zeros(1,totalVectorElements);
u = zeros(1,totalVectorElements);
v = zeros(1,totalVectorElements);
w = zeros(1,totalVectorElements); 
i = 1;

for ri = 1:1:size(r,2)
    for zi = 1:1:size(z,2)
        % ROI Elements
        rVector(1,i) = r(1,ri);
        zVector(1,i) = z(1,zi);
        % Electric Field Calculation
        electricFieldCap(zi,ri) = (potencialDifferenceCap*relativePermitivity)/parallelPlateCapDistance; %Volt/Meter
        rElectricVector(1,i) = electricFieldCap(zi,ri); 
        i = i + 1;
    end
end

figure
scatter3(rVector,zVector,rElectricVector);
xlabel('(r) radial position - meter') 
ylabel('(z) axial position - meter') 
zlabel('E(r) Electric Field Intensity - volt/meter')