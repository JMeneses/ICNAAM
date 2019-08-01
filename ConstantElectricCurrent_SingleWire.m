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
% >> Conductor - (ex:PPY polymer)
% >> Constant Electric Current in Finite Wire (ex: 1mA)
% >> Constant Static Potencial Difference Applyed to Finite Wire
% >> E(r,z) Field with radial(r) and length(z) components
% >> B(theta) Field with theta(t) component
%--------------------------------------------------------------------------
R = ((1/110)*wl)/(pi*wr); %Wire Resistance - Ohms [PPY 1.1 S/cm 110 S/m]
I = 0.001; %Current Intensity - Ampere 
pL = 0; %Wire Potencial Left - Volt
pR = R*I; %Wire Potencial Rigth - Volt
u0 = 4*pi*1.00000000082e-7; %Free Space Permeabilitty

rElectricFieldWireCurrent = zeros(size(z,2),size(r,2));
zElectricFieldWireCurrent = zeros(size(z,2),size(r,2));
tMagneticFieldWireCurrent = zeros(size(z,2),size(r,2));
rVector = zeros(1,totalVectorElements);
zVector = zeros(1,totalVectorElements);
rElectricVector = zeros(1,totalVectorElements);
zElectricVector = zeros(1,totalVectorElements);
tMagneticVector = zeros(1,totalVectorElements);
i = 1;
for ri = 1:1:size(r,2)
    for zi = 1:1:size(z,2)
        % ROI Elements
        rVector(1,i) = r(1,ri);
        zVector(1,i) = z(1,zi);
        % Electric and Magnetic Fields Calculation
        thetaL = atan(r(1,ri)/abs(-1*(wl/2)-z(1,zi)));
        thetaR = atan(r(1,ri)/abs((wl/2)-z(1,zi)));
        rElectricFieldWireCurrent(zi,ri) = ((-1/(log(wl/wr)))*(((R*I)/wl)*z(1,zi)-(R*I+2*pR)/2)/r(1,ri)); %Volt/Meter
        zElectricFieldWireCurrent(zi,ri) = ((R*I)/wl)*(log(wl/r(1,ri))/log(wl/wr)); %Volt/Meter
        tMagneticFieldWireCurrent(zi,ri) = ((u0*I)/(4*pi*r(1,ri)))*(cos(thetaR)+cos(thetaL));
        rElectricVector(1,i) = rElectricFieldWireCurrent(zi,ri);
        zElectricVector(1,i) = zElectricFieldWireCurrent(zi,ri);
        tMagneticVector(1,i) = tMagneticFieldWireCurrent(zi,ri);
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

figure
scatter3(rVector,zVector,tMagneticVector);
xlabel('(r) radial position - meter') 
ylabel('(z) axial position - meter') 
zlabel('H(theta) Magnetic Field Intensity - ampere/meter')
