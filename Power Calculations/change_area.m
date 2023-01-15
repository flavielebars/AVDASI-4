clear all;
close all;

%Author: Flavie Le Bars

%files outputed by ANSYS STK after uploading a CAD model of the satellite
%with the correct sizing of the solar panels. Very long process to upload a
%new CAD model with different panels sizing, hence this code recalculates
%the power produced by the solar panels without having to upload the
%updated CAD model.

%Power formula used in STK:
    %POWER PRODUCED = efficiency solar cells x solar intensity x effective area x solar irradiance

    % where effective area = area of the solar panels x cos(incidence angle(ie:theta))
    % solar intensity can be directly outputed by STK
    % solar irradiance is not directly ouputed by STK, hence it is calculated below
    % remain cst even when area changes are: solar irradiance, solar intensity, efficiency, incidence angle
    % change with the area: effective area & power produced
load("solar_intensity.mat");
solar_intensity = solar_intensity(1:length(solar_intensity)/2); %get the variable to be the correct length

  
%outputed by STK
load("AREA_effective.mat");
Ae = AREA(1:length(solar_intensity));
load("POWER.mat");
efficiency = 0.295; %efficieny of the Triple GAas solar cells
A = 17.89; %initial area calculated, the one wih which the STK simulation ran
theta = acosd(Ae/A);

denom = efficiency.* solar_intensity.* Ae;

%calculate solar irradiance
for m = 1:length(denom)
    solar_irradiance(m) = POWER(m)/denom(m);
    if solar_intensity(m)==0
        solar_irradiance(m) = 0;
    end
end

b=0;
k=0;
for i=1:length(theta)
    b = b+1;
    
    %stores theta when he incidence angle is below 30°; 0° means perfect
    %incidence with the sun, max power produced but also worst thermal
    %conditions so the solar panels are never orientated below 30°
    %from the sun
    if theta(b)<30
        theta_30(b)=30;
    end
end


%% CALCULATION POWER GENERATED WITH NEW AREA
area_new = 10.81782; %given in m^2


%calculate the new effective area
Area_effective= cosd(theta).*area_new;
T = load('timev4.mat'); %data file with the power modes for each instrument of the payload
load("total_req.mat"); % total power requirement calculated in "power_calculations.mat"

ECLIPSE = T.timelinev2.Eclipse;
d=1;
Ld = 0.4017;%lifetime degradation of the solar cells; calculated in "power_calculations.mat"
for k= 1:length(ECLIPSE)
        if (ECLIPSE(k)==0 || ECLIPSE(k)==1) %full or partial eclipse
            power_neww(k) = 0;
        else %full sun
            power_neww(k) = solar_irradiance(k).*efficiency.*Area_effective(k).*solar_intensity(k)*Ld;
            power_day(d) = power_neww(k);
            req_day(d)= payload_total_req(k);
            d = d+1;
        end
end