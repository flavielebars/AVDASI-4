close all;
clear all;

% Author: Flavie Le Bars
% Main outputs of this code:
    %1 - power requirements of the satellite (payload and subsystems) in
    %eclipse and daylight
    %2 - Calculation of the cell end of life degradation (Ld)
    %3 - Cell area required on the solar panels to meet the power
    %requirement
    %4 - Sizing of the battery packs


%Instruments of the satellite (ie: Payload): HGA, RSE, SRS, VENSAR, VENSPEC_H, VENSPEC_M
%& VENSPEC_U. The mode of each of these instruments is associated with a
%power mode (every 5 minute for 1 year) - this code calculates the power
%consumption of the payload and then all the subsystems every 5 minute

% Subsystems consuming power = AOCS, Thermal, Comms, OBDH, Power,
% Mechanisms, Propulsion

% all power given in Watts

% Read the power modes for each instrument of the satellite
T = load('timev4.mat');
HGA=T.timelinev2.HGA;
RSE=T.timelinev2.RSE;
SRS=T.timelinev2.SRS;
VENSAR=T.timelinev2.VENSAR;
VENSPEC_H=T.timelinev2.VENSPEC_H;
VENSPEC_M=T.timelinev2.VENSPEC_M;
VENSPEC_U=T.timelinev2.VENSPEC_U;

ECLIPSE = T.timelinev2.Eclipse;
time=T.timelinev2.Earth_Days;
date=T.timelinev2.Date;
count=length(HGA);



%% SRS
for (i=1:count)
    if    ismember(SRS(i), 'STANDBY')
        pow_SRS(i,:) = 19;
        
    elseif ismember(SRS(i), 'HIGH_DENSITY')
        pow_SRS(i,:) = 115; 

    elseif ismember(SRS(i), 'LOW_DENSITY')
        pow_SRS(i,:) = 115; 

    elseif ismember(SRS(i), 'WARMUP')
        pow_SRS(i,:) = 19; 
    end
end

%% RSE
for (i=1:count)
    if    ismember(RSE(i), 'GRAVITY')
        pow_RSE(i,:) = 11.7;
        
    elseif ismember(RSE(i), 'OCCULTATION')
        pow_RSE(i,:) = 5.15;

    elseif ismember(RSE(i), 'OFF')
        pow_RSE(i,:) = 0;
    end
end


%% VENSAR
for (i=1:count)
    if    ismember(VENSAR(i), 'ALTRAD')
        pow_VENSAR(i,:) = 250.5;
        
    elseif ismember(VENSAR(i), 'OFF')
        pow_VENSAR(i,:) = 0;

    elseif ismember(VENSAR(i), 'RAD')
        pow_VENSAR(i,:) = 147.2;

    elseif ismember(VENSAR(i), 'RAD_OFFPOINT')
        pow_VENSAR(i,:) = 147.2;

     elseif ismember(VENSAR(i), 'STANDBY')
        pow_VENSAR(i,:) = 127.7;

     elseif ismember(VENSAR(i), 'VN1')
        pow_VENSAR(i,:) = 1364.2;
    end
end


%% VENSPEC H
for (i=1:count)
    if    ismember(VENSPEC_H(i), 'OFF')
        pow_VENSPEC_H(i,:) = 0;
        
    elseif ismember(VENSPEC_H(i), 'PRECOOL')
        pow_VENSPEC_H(i,:) = 28.8;

    elseif ismember(VENSPEC_H(i), 'SCIENCE_DAY')
        pow_VENSPEC_H(i,:) = 27.5;
    
    elseif ismember(VENSPEC_H(i), 'SCIENCE_NIGHT')
        pow_VENSPEC_H(i,:) = 27.5;
    
    elseif ismember(VENSPEC_H(i), 'SCIENCE_OFFPOINT')
        pow_VENSPEC_H(i,:) = 27.5;
    
    elseif ismember(VENSPEC_H(i), 'STANDBY')
        pow_VENSPEC_H(i,:) = 15.4;
    end
end


%% VENSPECM
for (i=1:count)
    if    ismember(VENSPEC_M(i), 'OFF')
        pow_VENSPEC_M(i,:) = 0;
        
    elseif ismember(VENSPEC_M(i), 'SCIENCE')
        pow_VENSPEC_M(i,:) = 16.6; %TBF

    elseif ismember(VENSPEC_M(i), 'STANDBY')
        pow_VENSPEC_M(i,:) = 8.4;
    end
end


%% VENSPECU
for (i=1:count)
    if    ismember(VENSPEC_U(i), 'OFF')
        pow_VENSPEC_U(i,:) = 0;
        
    elseif ismember(VENSPEC_U(i), 'PRECOOL')
        pow_VENSPEC_U(i,:) = 13.8;

    elseif ismember(VENSPEC_U(i), 'STANDBY')
        pow_VENSPEC_U(i,:) = 7.4;
    
    elseif ismember(VENSPEC_U(i), 'VENUS_OBSERVATION')
        pow_VENSPEC_U(i,:) = 13.8;
    
    elseif ismember(VENSPEC_U(i), 'VENUS_OBSERVATION_OFFPOINT')
        pow_VENSPEC_U(i,:) = 13.8; 
    end
end

%% ALL PAYLOAD

%sum payload requirement every 5 minute
pow_tot_5min = (pow_SRS + pow_VENSAR + pow_VENSPEC_H + pow_VENSPEC_M + pow_VENSPEC_U + pow_RSE); 
pow_tot_5min_margin = pow_tot_5min*1.2; %Add a 20% margin following the requirements guidelines

peak_payload = max(pow_tot_5min_margin);
mean_payload = mean(pow_tot_5min_margin);


%% HGA (communication subsystem)
for (i=1:count)
    if  ismember(HGA(i), 'OFF')
        pow_HGA(i,:) = 0;
        
    elseif ismember(HGA(i), 'HGA_DUMP')
        pow_HGA(i,:) = 410;

    elseif ismember(HGA(i), 'HGA_RSE_OCCULT')
        pow_HGA(i,:) = 410; 
    end
end

%% SOLAR ARRAY AREA REQ
% Eclipse: 0 full eclipse, 1 partial eclipse, 2  sun (as shown on the
% provided timeline)
i = 1;
e= 1;

for (k=1:count)
    if  ECLIPSE(k)==2 %no eclipse
        Pd(i) = pow_tot_5min_margin(k) + pow_HGA(k)*1.05; %power required by the payload and the communication system when no eclipse
        i = i+1;
        
    else %full or partial eclipse
        Pe(e) = pow_tot_5min_margin(k) + pow_HGA(k)*1.05; %power required by the payload and the communication system when full or partial eclipse
        e = e+1;
    end
end

%% Eclipse & daylight length calculation
j = 0;
m = 0;
k=1;
length_data = 64316; %number of rows in the data file; ie: length given is 64316*5 in minutes

while k<length_data
    w=0;
    b=0;
        while (ECLIPSE(k)==0 || ECLIPSE(k)==1) && k<length_data  %full or partial eclipse
            w = w+1;
            k=k+1;
        end
       j=j+1;
       eclipse(j,:)=w;
       eclipse_index(j,:)=k;
       k=k+1;  
       while (ECLIPSE(k)==2) && k<length_data  %no eclipse
            b=b+1;
            k=k+1;
       end
    m=m+1;   
    no_eclipsee(m,:)=b;
    no_eclipse_index(m,:)=k;

    %The two longest periods of daylight last 6573 min and 6971 min; the
    %code below outputs the payload power mode for the length of these full
    %sun periods
        if(no_eclipsee(m,:) == 6573)  
             ind = k-no_eclipsee(m,:);
             long_eclipse1 = T.timelinev2(ind:k,25:31);
        end
        if(no_eclipsee(m,:) == 6971)  
             ind = k-no_eclipsee(m,:);
             long_eclipse2 = T.timelinev2(ind:k,25:31);
        end
end


%% POWER REQUIREMENTS INCLUDING PAYLOAD AND SUBSYSTEMS

%5% magin added to subsystems (all subsystems except thermal and comms)
subsytem = 175; 

% thermal power requirements (5% margin)
thermal_sun = 222.78*1.05;
thermal_eclipse_peak = 911.47*1.05;
thermal_eclipse = 545*1.05;

e=1;
d=1;

%theta is the incidence angle between the solar panels and the sun (angle
%greatly affecting the power produced by the solar panels). It is
%calculated in the "change_area" matlab code
load("theta.mat");


for k= 1:64316
    if (ECLIPSE(k)==0 ||ECLIPSE(k)==1) %full or partial eclipse
        total_req(k) = subsytem*1.1 + pow_tot_5min_margin(k) + pow_HGA(k)*1.05 + thermal_eclipse_peak; %thermal peak to accoun for worst case scenario
        Pe_all(e)=total_req(k); %total power required
        e=e+1;
        
   
    elseif (ECLIPSE(k)==2) %no eclipse - full sun
        total_req(k) = subsytem*1.1 + pow_tot_5min_margin(k) + pow_HGA(k)*1.05 + thermal_sun;
        Pd_all(d)=total_req(k);%total power required
        d=d+1;
        theta_day(d) = theta(k);%incidence angle full sun

    end
end


% timea = 1:length(theta);
% timed = 1:length(theta_day);
% plot(timed,theta_day)
% figure
% plot(time,total_req)
% xlabel("time")
% ylabel("Power generated (W)")
% legend ("Power requirement of the payload and subsytem for the timeline")


%% Array sizing

%Pe power req eclipse
%Pd power req daylight
%Xe efficiency of the paths from solar arrays through the batteries to
%individual load paths 
%Xd same as Xe without the batteries
Xe = 0.65;
Xd = 0.85;

%Average time per orbit in minute
To = 94.8877;

%calculate average time in daylight and eclipse
Td = sum(no_eclipsee,'all')*5;
Te = sum(eclipse,'all')*5;
time_orbit = count*5/To;
Td_mean = Td/time_orbit;
Te_mean = Te/time_orbit;




Pe_peak = (max(Pe) + thermal_eclipse + subsytem*1.1)*1.2;
Pd_peak = (max(Pd) + thermal_sun + subsytem*1.1)*1.2;

%Psa is the power that the solar array must provide during daylight to
%power the spacecraft for the entire orbit
Psa = ((mean(Pe_all)*Te/Xe)+(mean(Pd_all)*Td/Xd))/Td;


%Po cell output performance per unit area W/m^2
irrad = 2620; %solar irradiation near Venus in W/m^2
efficiency = 0.295; %Mulijunction GaAs
Po = irrad*efficiency;


Id = 0.77;%inherent degradation
degradation = 0.32375*0.5; %Based on Magellan values
year = 4+7/12; %science mission and 7 months around venus
year2 = 2; %2 years before arriving near venus
degradation2 = 0.05;
Pbol = Po*Id*cosd(31.36); %arrays power per unit area  --- beginning of life
Ld = (1-degradation2)^year2 * (1-degradation)^year;
Peol = Pbol*Ld; %Power produced at the end of life; taking the cells degradation into account

Asa = Psa/Peol; % array sizing in m^2


%% Batterie sizing
DOD = 0.2; %determined on the number of cycles throughout the mission
N=2; %nb of batteries
n=0.75; %batterie to load efficiency


Vd = 28;
C = (max(Pe_all)*(max(eclipse)*5)/60)/(n*Vd*DOD);

load("Billy_Eclipse_Summary.mat");%file with data on eclipse obtained from STK
max_eclipse = max(BillyEclipseSummary)/60; % calculate max eclipse time bc the batteries are sized for the worst case scenario

Pe_peak = 2453.3;
Cr = ((Pe_peak)*((max_eclipse)/60))/(DOD*N*n); %batt capacity Whr

vbus = 28;
Cr_amp = Cr/vbus;
vcell = 3.65; %datasheet
Ns = vbus/vcell;%nb of batteries in series

ccell = 7.5;
Nb = Cr_amp/ccell;
one_batt = 0.38; %mass one battery in kg
Batt_mass = 2*32*one_batt*1.2; %add a 20% margin to account for the box to store the batteries
