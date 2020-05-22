%% Clear ALL
clear all
close all
clc

%% Import Data
% NEDC
nedc = xlsread("NEDC_WLTP_data.xlsx","NEDC");
time_NEDC = nedc(:, 1); % [s]
speed_NEDC = nedc(:, 2); % [km/h]
gear_NEDC = nedc(:, 3); % [-]
figure(1)
yyaxis left
plot(time_NEDC, speed_NEDC, 'Linewidth', 2)
ylim([0 130])
yyaxis right
plot(time_NEDC, gear_NEDC)
ylim([0 6])
title("NEDC Speed/Gear Profile")
legend("Speed", "Gear")
grid on

% Vehicle 1 WLTP
V1_wltp = xlsread("NEDC_WLTP_data.xlsx","WLTP_Vehicle_1");
time_V1_wltp = V1_wltp(:, 1); % [s]
speed_V1_wltp = V1_wltp(:, 2); % [km/h]
gear_V1_wltp = V1_wltp(:, 3); % [-]

% Vehicle 2 WLTP
V2_wltp = xlsread("NEDC_WLTP_data.xlsx","WLTP_Vehicle_2");
time_V2_wltp = V2_wltp(:, 1); % [s]
speed_V2_wltp = V2_wltp(:, 2); % [km/h]
gear_V2_wltp = V2_wltp(:, 3); % [-]

figure(2)
yyaxis left
plot(time_V1_wltp, speed_V1_wltp, 'Linewidth', 2)
ylim([0 150])
yyaxis right
plot(time_V1_wltp, gear_V1_wltp, 'r')
hold on
plot(time_V1_wltp, gear_V2_wltp, 'g')
hold off
ylim([0 6])
title("WLTP Speed/Gear Profile")
legend("Speed", "Idea Gear", "Punto Gear")
grid on

V1_characteristics
V2_characteristics

%% rpm
% Vehicle 1 NEDC
gear_ratio = zeros(length(gear_NEDC), 1);
gear_ratio(gear_NEDC ~= 0) = V1_tau_g(gear_NEDC(gear_NEDC ~= 0)); % [-]
V1_rpm = V1_omega_min*ones(length(nedc), 1);
V1_rpm(gear_NEDC ~= 0) = speed_NEDC(gear_NEDC ~= 0)*60/2/pi/V1_wheel_radius/3.6*V1_finalDrive_ratio.*gear_ratio(gear_NEDC ~= 0);
V1_rpm(V1_rpm < V1_omega_min) = V1_omega_min; % Engine Speed [rpm]

% Vehicle 2 NEDC
V2_rpm = V2_omega_min*ones(length(nedc), 1);
V2_rpm(gear_NEDC ~= 0) = speed_NEDC(gear_NEDC ~= 0)*60/2/pi/V2_wheel_radius/3.6*V2_finalDrive_ratio.*gear_ratio(gear_NEDC ~= 0);
V2_rpm(V2_rpm < V2_omega_min) = V2_omega_min; % Engine Speed [rpm]

% Vehicle 1 WLTP
V1_wltp_gear_ratio = zeros(length(gear_V1_wltp), 1);
V1_wltp_gear_ratio(gear_V1_wltp ~= 0) = V1_tau_g(gear_V1_wltp(gear_V1_wltp ~= 0)); % [-]
V1_rpm_wltp = V1_omega_min*ones(length(V1_wltp), 1);
V1_rpm_wltp(gear_V1_wltp ~= 0) = speed_V1_wltp(gear_V1_wltp ~= 0)*60/2/pi/...
    V1_wheel_radius/3.6*V1_finalDrive_ratio.*V1_wltp_gear_ratio(gear_V1_wltp ~= 0);
V1_rpm_wltp(V1_rpm_wltp < V1_omega_min) = V1_omega_min; % Engine Speed [rpm]

% Vehicle 2 WLTP
V2_wltp_gear_ratio = zeros(length(gear_V2_wltp), 1);
V2_wltp_gear_ratio(gear_V2_wltp ~= 0) = V1_tau_g(gear_V2_wltp(gear_V2_wltp ~= 0)); % [-]
V2_rpm_wltp = V2_omega_min*ones(length(V2_wltp), 1);
V2_rpm_wltp(gear_V2_wltp ~= 0) = speed_V2_wltp(gear_V2_wltp ~= 0)*60/2/pi/...
    V2_wheel_radius/3.6*V2_finalDrive_ratio.*V2_wltp_gear_ratio(gear_V2_wltp ~= 0);
V2_rpm_wltp(V2_rpm_wltp < V2_omega_min) = V2_omega_min; % Engine Speed [rpm]
%% BMEP
% Vehicle 1 NEDC
V1_Fres_nedc = V1_NEDC_F0 + V1_NEDC_F1*speed_NEDC + V1_NEDC_F2*(speed_NEDC).^2; % Resistive Force [N]
V1_m_trasl_nedc = V1_NEDC_mass + V1_wheels_inertia/V1_wheel_radius^2 + V1_engine_inertia/...
    V1_wheel_radius^2*V1_finalDrive_ratio^2.*gear_ratio.^2; % Apparent Mass [kg]
V1_acc_nedc = [0; diff(speed_NEDC/3.6)]; % Acceleration [m/s^2]
V1_power_nedc = (V1_Fres_nedc + V1_m_trasl_nedc.*V1_acc_nedc).*speed_NEDC/3.6; % Vehicle Motive Power [W]
V1_bmep_nedc = 1200*V1_power_nedc*1e-3./(V1_rpm*V1_engine_displacement*V1_gear_eff); % Engine BMEP [bar]
V1_bmep_nedc(V1_bmep_nedc < 0) = 0;

% Vehicle 2 NEDC
V2_Fres_nedc = V2_NEDC_F0 + V2_NEDC_F1*speed_NEDC + V2_NEDC_F2*(speed_NEDC).^2; % Resistive Force [N]
V2_m_trasl_nedc = V2_NEDC_mass + V2_wheels_inertia/V2_wheel_radius^2 + V2_engine_inertia/...
    V2_wheel_radius^2*V2_finalDrive_ratio^2.*gear_ratio.^2; % Apparent Mass [kg]
V2_acc_nedc = [0; diff(speed_NEDC/3.6)]; % Acceleration [m/s^2]
V2_power_nedc = (V2_Fres_nedc + V2_m_trasl_nedc.*V2_acc_nedc).*speed_NEDC/3.6; % Vehicle Motive Power [W]
V2_bmep_nedc = 1200*V2_power_nedc*1e-3./(V2_rpm*V1_engine_displacement*V2_gear_eff); % Engine BMEP [bar]
V2_bmep_nedc(V2_bmep_nedc < 0) = 0;

% Vehicle 1 WLTP
V1_Fres_wltp = V1_WLTP_F0 + V1_WLTP_F1*speed_V1_wltp + V1_WLTP_F2*(speed_V1_wltp).^2; % Resistive Force [N]
V1_m_trasl_wltp = V1_WLTP_mass + V1_wheels_inertia/V1_wheel_radius^2 + V1_engine_inertia/...
    V1_wheel_radius^2*V1_finalDrive_ratio^2.*V1_wltp_gear_ratio.^2; % Apparent Mass [kg]
V1_acc_wltp = [0; diff(speed_V1_wltp/3.6)]; % Acceleration [m/s^2]
V1_power_wltp = (V1_Fres_wltp + V1_m_trasl_wltp.*V1_acc_wltp).*speed_V1_wltp/3.6; % Vehicle Motive Power [W]
V1_bmep_wltp = 1200*V1_power_wltp*1e-3./(V1_rpm_wltp*V1_engine_displacement*V1_gear_eff); % Engine BMEP [bar]
V1_bmep_wltp(V1_bmep_wltp < 0) = 0;

% Vehicle 2
V2_Fres_wltp = V2_WLTP_F0 + V2_WLTP_F1*speed_V2_wltp + V2_WLTP_F2*(speed_V2_wltp).^2; % Resistive Force [N]
V2_m_trasl_wltp = V2_WLTP_mass + V2_wheels_inertia/V2_wheel_radius^2 + V2_engine_inertia/...
    V2_wheel_radius^2*V2_finalDrive_ratio^2.*V2_wltp_gear_ratio.^2; % Apparent Mass [kg]
V2_acc_wltp = [0; diff(speed_V2_wltp/3.6)]; % Acceleration [m/s^2]
V2_power_wltp = (V2_Fres_wltp + V2_m_trasl_wltp.*V2_acc_wltp).*speed_V2_wltp/3.6; % Vehicle Motive Power [W]
V2_bmep_wltp = 1200*V2_power_wltp*1e-3./(V2_rpm_wltp*V1_engine_displacement*V2_gear_eff); % Engine BMEP [bar]
V2_bmep_wltp(V2_bmep_wltp < 0) = 0;

%% Definition of WOT performance
rpm_max = [0.8; 1; 1.25; 1.5; 1.75; 2; 2.5; 3; 3.5; 4; 4.5; 5]*1e3; % [rpm]
bmep_max = [6.05; 8.73; 12.08; 16.2; 18.12; 18.12; 17.12; 15.61; 13.97; 12.38; 10.52; 8.57]; % [bar]

% Calibration to go below limit
for i = 1:length(V1_bmep_wltp)
    for j = 1:length(bmep_max)
       if V1_rpm_wltp(i) > rpm_max(j) && V1_rpm_wltp(i) < rpm_max(j+1) &&...
               V1_bmep_wltp(i) > interp1(rpm_max, bmep_max, V1_rpm_wltp(i))
          V1_bmep_wltp(i) = interp1(rpm_max, bmep_max, V1_rpm_wltp(i));
          break;
       end
    end
end

for i = 1:length(V2_bmep_wltp)
    for j = 1:length(bmep_max)
       if V2_rpm_wltp(i) > rpm_max(j) && V2_rpm_wltp(i) < rpm_max(j+1) &&...
               V2_bmep_wltp(i) > interp1(rpm_max, bmep_max, V2_rpm_wltp(i))
          V2_bmep_wltp(i) = interp1(rpm_max, bmep_max, V2_rpm_wltp(i)); 
          break;
       end
    end
end
%% Calculation of fuel consumption
% Get fuel consumption map
pq = dlmread("fuel_cons_kg_h_vers_matlab.pqm", "\t");
rpm_pq = pq(1, 2:end); % [rpm]
bmep_pq = pq(2:end, 1); % [bar]
fc_pq = pq(2:end, 2:end); % [kg/h]

% Interpolation to get the engine fuel rate
% Vehicle 1 NEDC
V1_fc_nedc = interp2(rpm_pq, bmep_pq,fc_pq, V1_rpm, V1_bmep_nedc); % [kg/h]
% Map Limitations
V1_fc_nedc(V1_rpm < 850 | V1_bmep_nedc < 0.5) = V1_fuel_const_idle*...
    V1_rpm(V1_rpm < 850 | V1_bmep_nedc < 0.5)/V1_omega_min*1e-3; % [kg/h]
% Fuel Cut-Off
V1_fc_nedc(V1_power_nedc < 0 & V1_rpm > 1000) = 0; % [kg/h]
% Stop-Start
V1_fc_nedc(speed_NEDC <= 0 & gear_ratio <= 0) = 0; % [kg/h]

% Vehicle 2 NEDC
V2_fc_nedc = interp2(rpm_pq, bmep_pq,fc_pq, V2_rpm, V2_bmep_nedc); % [kg/h]
% Map Limitations
V2_fc_nedc(V2_rpm < 850 | V2_bmep_nedc < 0.5) = V2_fuel_const_idle*...
    V2_rpm(V2_rpm < 850 | V2_bmep_nedc < 0.5)/V2_omega_min*1e-3; % [kg/h]
% Fuel Cut-Off
V2_fc_nedc(V2_power_nedc < 0 & V2_rpm > 1000) = 0; % [kg/h]
% Stop-Start
V2_fc_nedc(speed_NEDC == 0 & gear_ratio == 0) = 0; % [kg/h]

% Vehicle 1 WLTP
V1_fc_wltp = interp2(rpm_pq, bmep_pq,fc_pq, V1_rpm_wltp, V1_bmep_wltp); % [kg/h]
% Map Limitations
V1_fc_wltp(V1_rpm_wltp < 850 | V1_bmep_wltp < 0.5) = V1_fuel_const_idle*...
    V1_rpm_wltp(V1_rpm_wltp < 850 | V1_bmep_wltp < 0.5)/V1_omega_min*1e-3; % [kg/h]
V1_fc_wltp(V1_bmep_wltp > 18.12) = 4.7074148; % [kg/h]
% Fuel Cut-Off
V1_fc_wltp(V1_power_wltp < 0 & V1_rpm_wltp > 1000) = 0; % [kg/h]
% Stop-Start
V1_fc_wltp(speed_V1_wltp == 0 & V1_wltp_gear_ratio == 0) = 0; % [kg/h]

% Vehicle 2 WLTP
V2_fc_wltp = interp2(rpm_pq, bmep_pq,fc_pq, V2_rpm_wltp, V2_bmep_wltp); % [kg/h]
% Map Limitations
V2_fc_wltp(V2_rpm_wltp < 850 | V2_bmep_wltp < 0.5) = V2_fuel_const_idle*...
    V2_rpm_wltp(V2_rpm_wltp < 850 | V2_bmep_wltp < 0.5)/V2_omega_min*1e-3; % [kg/h]
% Fuel Cut-Off
V2_fc_wltp(V2_power_wltp < 0 & V2_rpm_wltp > 1000) = 0; % [kg/h]
% Stop-Start
V2_fc_wltp(speed_V2_wltp == 0 & V2_wltp_gear_ratio == 0) = 0; % [kg/h]

% NOx rate
nox_map = dlmread("nox_g_h.txt", "\t");
rpm_nox = nox_map(1, 2:end); % [rpm]
bmep_nox = nox_map(2:end, 1); % [bar]
nox_pq = nox_map(2:end, 2:end); % [g/h]

% Vehicle 1 NEDC
V1_nox_nedc = interp2(rpm_nox, bmep_nox, nox_pq, V1_rpm, V1_bmep_nedc); % [g/h]
% Map Limitations
V1_nox_nedc(V1_rpm < 1000 | V1_bmep_nedc < 1) = 1.2; % [g/h]
% Fuel Cut-Off
V1_nox_nedc(V1_power_nedc < 0 & V1_rpm > 1000) = 0; % [g/h]
% Stop-Start
V1_nox_nedc(speed_NEDC == 0 & gear_ratio == 0) = 0; % [g/h]

% Vehicle 2 NEDC
V2_nox_nedc = interp2(rpm_nox, bmep_nox, nox_pq, V2_rpm, V2_bmep_nedc); % [g/h]
% Map Limitations
V2_nox_nedc(V2_rpm < 1000 | V2_bmep_nedc < 1) = 1.2; % [g/h]
% Fuel Cut-Off
V2_nox_nedc(V1_power_nedc < 0 & V2_rpm > 1000) = 0; % [g/h]
% Stop-Start
V2_nox_nedc(speed_NEDC == 0 & gear_ratio == 0) = 0; % [g/h]

% Vehicle 1 WLTP
V1_nox_wltp = interp2(rpm_nox, bmep_nox, nox_pq, V1_rpm_wltp, V1_bmep_wltp); % [g/h]
% Map Limitations
V1_nox_wltp(V1_rpm_wltp < 1000 | V1_bmep_wltp < 1) = 1.2; % [g/h]
% Fuel Cut-Off
V1_nox_wltp(V1_power_wltp < 0 & V1_rpm_wltp > 1000) = 0; % [g/h]
% Stop-Start
V1_nox_wltp(speed_V1_wltp == 0 & V1_wltp_gear_ratio == 0) = 0; % [g/h]

% Vehicle 2 WLTP
V2_nox_wltp = interp2(rpm_nox, bmep_nox, nox_pq, V2_rpm_wltp, V2_bmep_wltp); % [g/h]
% Map Limitations
V2_nox_wltp(V2_rpm_wltp < 1000 | V2_bmep_wltp < 1) = 1.2; % [g/h]
% Fuel Cut-Off
V2_nox_wltp(V1_power_wltp < 0 & V2_rpm_wltp > 1000) = 0; % [g/h]
% Stop-Start
V2_nox_wltp(speed_V2_wltp == 0 & V2_wltp_gear_ratio == 0) = 0; % [g/h]

% Cumulative fuel consumption
% Vehicle 1
V1_fc_tot_nedc = cumtrapz(time_NEDC, V1_fc_nedc)/3.6; % [g]
V1_fc_tot_wltp = cumtrapz(time_V1_wltp, V1_fc_wltp)/3.6; % [g]

% V1_fc_tot_wltp = cumtrapz(time_V1_wltp, V1_fc_wltp)/3.6; % [g]

% Vehicle 2
V2_fc_tot_nedc = cumtrapz(time_NEDC, V2_fc_nedc)/3.6; % [g]
V2_fc_tot_wltp = cumtrapz(time_V2_wltp, V2_fc_wltp)/3.6; % [g]

% Cumulative NOx emissions
% Vehicle 1
V1_nox_tot_nedc = cumtrapz(time_NEDC, V1_nox_nedc)/3600; % [g]
V1_nox_tot_wltp = cumtrapz(time_V1_wltp, V1_nox_wltp)/3600; % [g]

% Vehicle 2
V2_nox_tot_nedc = cumtrapz(time_NEDC, V2_nox_nedc)/3600; % [g]
V2_nox_tot_wltp = cumtrapz(time_V2_wltp, V2_nox_wltp)/3600; % [g]

% Distance travelled
km_tot_nedc = cumtrapz(time_NEDC, speed_NEDC)/3600; % [km]
V1_km_tot_wltp = cumtrapz(time_V1_wltp, speed_V1_wltp)/3600; % [km]
V2_km_tot_wltp = cumtrapz(time_V2_wltp, speed_V2_wltp)/3600; % [km]

% Fuel economy
% Vehicle 1
V1_l_100km_nedc = V1_fc_tot_nedc(end)/V1_fuel_density/km_tot_nedc(end)*100; % [l/100km]
V1_l_100km_wltp = V1_fc_tot_wltp(end)/V1_fuel_density/V1_km_tot_wltp(end)*100; % [l/100km]

% Vehicle 2
V2_l_100km_nedc = V2_fc_tot_nedc(end)/V2_fuel_density/km_tot_nedc(end)*100; % [l/100km]
V2_l_100km_wltp = V2_fc_tot_wltp(end)/V2_fuel_density/V2_km_tot_wltp(end)*100; % [l/100km]

% NOx emisisons
% Vehicle 1
V1_nox_spec_nedc = V1_nox_tot_nedc(end)/km_tot_nedc(end); % [g/km]
V1_nox_spec_wltp = V1_nox_tot_wltp(end)/V1_km_tot_wltp(end); % [g/km]

% Vehicle 2
V2_nox_spec_nedc = V2_nox_tot_nedc(end)/km_tot_nedc(end); % [g/km]
V2_nox_spec_wltp = V2_nox_tot_wltp(end)/V2_km_tot_wltp(end); % [g/km]

% CO2 emissions
% Vehicle 1
V1_mCO2_nedc = V1_fuel_density*V1_l_100km_nedc/0.0315*1e-3; % [g/km]
V1_mCO2_wltp = V1_fuel_density*V1_l_100km_wltp/0.0315*1e-3; % [g/km]

% Vehicle 2
V2_mCO2_nedc = V2_fuel_density*V2_l_100km_nedc/0.0315*1e-3; % [g/km]
V2_mCO2_wltp = V2_fuel_density*V2_l_100km_wltp/0.0315*1e-3; % [g/km]

% Sumarizing Data in Table
FC_NEDC = [V1_fc_tot_nedc(end); V2_fc_tot_nedc(end)];
FC_WLTP = [V1_fc_tot_wltp(end); V2_fc_tot_wltp(end)];
NOX_NEDC = [V1_nox_tot_nedc(end); V2_nox_tot_nedc(end)];
NOX_WLTP = [V1_nox_tot_wltp(end); V2_nox_tot_wltp(end)];
FE_NEDC = [V1_l_100km_nedc; V2_l_100km_nedc];
FE_WLTP = [V1_l_100km_wltp; V2_l_100km_wltp];
NOX_SPEC_NEDC = [V1_nox_spec_nedc; V2_nox_spec_nedc];
NOX_SPEC_WLTP = [V1_nox_spec_wltp; V2_nox_spec_wltp];
CO2_NEDC = [V1_mCO2_nedc; V2_mCO2_nedc];
CO2_WLTP = [V1_mCO2_wltp; V2_mCO2_wltp];

Data = table(FC_NEDC, FC_WLTP, NOX_NEDC, NOX_WLTP, FE_NEDC, FE_WLTP,...
    NOX_SPEC_NEDC, NOX_SPEC_WLTP, CO2_NEDC, CO2_WLTP)

%% Mechanical Energy and Maximum Energy Recovered Through Regenerative Braking
% NEDC
for i = 1:length(V1_power_nedc)
   if V1_power_nedc(i) > 0
       V1_pwr_pos_nedc(i) = V1_power_nedc(i);
   elseif V1_power_nedc(i) < 0
       V1_pwr_neg_nedc(i) = V1_power_nedc(i);
   end
   
   if V2_power_nedc(i) > 0
       V2_pwr_pos_nedc(i) = V2_power_nedc(i);
   elseif V2_power_nedc(i) < 0
       V2_pwr_neg_nedc(i) = V2_power_nedc(i);
   end
end

% WLTP
for i = 1:length(V1_power_wltp)
   if V1_power_wltp(i) > 0
       V1_pwr_pos_wltp(i) = V1_power_wltp(i);
   elseif V1_power_wltp(i) < 0
       V1_pwr_neg_wltp(i) = V1_power_wltp(i);
   end
   
   if V2_power_wltp(i) > 0
       V2_pwr_pos_wltp(i) = V2_power_wltp(i);
   elseif V2_power_wltp(i) < 0
       V2_pwr_neg_wltp(i) = V2_power_wltp(i);
   end
end


% Vehicle 1
V1_mechanical_energy_nedc = zeros(length(V1_pwr_pos_nedc), 1);
for i = 2:length(V1_pwr_pos_nedc)
   V1_mechanical_energy_nedc(i) = V1_mechanical_energy_nedc(i-1) + V1_pwr_pos_nedc(i); % [J]
end
V1_specific_mechanical_energy_nedc = V1_mechanical_energy_nedc(end)/(1000*km_tot_nedc(end)); % [kJ/km]

V1_brake_energy_nedc = zeros(length(V1_pwr_neg_nedc), 1);
for i = 2:length(V1_pwr_neg_nedc)
   V1_brake_energy_nedc(i) = V1_brake_energy_nedc(i-1) + V1_pwr_neg_nedc(i);
end
V1_specific_brake_energy_nedc = V1_brake_energy_nedc(end)/(1000*km_tot_nedc(end)); % [kJ/km]

V1_mechanical_energy_wltp = zeros(length(V1_pwr_pos_wltp), 1);
for i = 2:length(V1_pwr_pos_wltp)
   V1_mechanical_energy_wltp(i) = V1_mechanical_energy_wltp(i-1) + V1_pwr_pos_wltp(i); % [J]
end
V1_specific_mechanical_energy_wltp = V1_mechanical_energy_wltp(end)/(1000*V1_km_tot_wltp(end)); % [kJ/km]

V1_brake_energy_wltp = zeros(length(V1_pwr_neg_wltp), 1);
for i = 2:length(V1_pwr_neg_wltp)
   V1_brake_energy_wltp(i) = V1_brake_energy_wltp(i-1) + V1_pwr_neg_wltp(i);
end
V1_specific_brake_energy_wltp = V1_brake_energy_wltp(end)/(1000*V1_km_tot_wltp(end)); % [kJ/km]

% Vehicle 2
V2_mechanical_energy_nedc = zeros(length(V2_pwr_pos_nedc), 1);
for i = 2:length(V2_pwr_pos_nedc)
   V2_mechanical_energy_nedc(i) = V2_mechanical_energy_nedc(i-1) + V2_pwr_pos_nedc(i); % [J]
end
V2_specific_mechanical_energy_nedc = V2_mechanical_energy_nedc(end)/(1000*km_tot_nedc(end)); % [kJ/km]

V2_brake_energy_nedc = zeros(length(V2_pwr_neg_nedc), 1);
for i = 2:length(V2_pwr_neg_nedc)
   V2_brake_energy_nedc(i) = V2_brake_energy_nedc(i-1) + V2_pwr_neg_nedc(i);
end
V2_specific_brake_energy_nedc = V2_brake_energy_nedc(end)/(1000*km_tot_nedc(end)); % [kJ/km]

V2_mechanical_energy_wltp = zeros(length(V2_pwr_pos_wltp), 1);
for i = 2:length(V2_pwr_pos_wltp)
   V2_mechanical_energy_wltp(i) = V2_mechanical_energy_wltp(i-1) + V2_pwr_pos_wltp(i); % [J]
end
V2_specific_mechanical_energy_wltp = V2_mechanical_energy_wltp(end)/(1000*V2_km_tot_wltp(end)); % [kJ/km]

V2_brake_energy_wltp = zeros(length(V2_pwr_neg_wltp), 1);
for i = 2:length(V2_pwr_neg_wltp)
   V2_brake_energy_wltp(i) = V2_brake_energy_wltp(i-1) + V2_pwr_neg_wltp(i);
end
V2_specific_brake_energy_wltp = V2_brake_energy_wltp(end)/(1000*V2_km_tot_wltp(end)); % [kJ/km]

%% Post-Processing
% Engine Speed vs Time
figure(3)
subplot(2, 1, 1)
plot(time_NEDC, V1_rpm, time_NEDC, V2_rpm, 'LineWidth', 1.5)
grid on
xlabel("Time [s]")
ylabel("Engine Speed [rpm]")
legend("Idea", "Punto")
title("Engine RPM vs NEDC Time")

subplot(2, 1, 2)
plot(time_V1_wltp, V1_rpm_wltp, time_V2_wltp, V2_rpm_wltp, 'LineWidth', 1.5)
grid on
xlabel("Time [s]")
ylabel("Engine Speed [rpm]")
legend("Idea", "Punto")
title("Engine RPM vs WLTP Time")

% BMEP vs Time
figure(4)
subplot(2, 1, 1)
plot(time_NEDC, V1_bmep_nedc, time_NEDC, V2_bmep_nedc, 'LineWidth', 1.5)
grid on
xlabel("Time [s]")
ylabel("PME [bar]")
legend("Idea", "Punto")
title("BMEP vs NEDC Time")

subplot(2, 1, 2)
plot(time_V1_wltp, V1_bmep_wltp, time_V2_wltp, V2_bmep_wltp, 'LineWidth', 1.5)
grid on
xlabel("Time [s]")
ylabel("PME [bar]")
legend("Idea", "Punto")
title("BMEP vs WLTP Time")

% Instantaneous Fuel Consumption vs Time
figure(5)
subplot(2, 1, 1)
yyaxis left
plot(time_NEDC, V1_fc_nedc, '-b', time_NEDC, V2_fc_nedc, '-r', 'LineWidth', 1.5)
yyaxis right
plot( time_NEDC, V1_fc_tot_nedc, '-y', time_NEDC, V2_fc_tot_nedc, '-g', 'LineWidth', 1.5)
grid on
xlabel("Time [s]")
ylabel("Fuel Consumption [kg/h]")
legend("Idea Instantaneous", "Punto Instantaneous", "Idea Cumulative", ...
    "Punto Cumulative")
title("Fuel Consumption vs NEDC Time")

subplot(2, 1, 2)
yyaxis left
plot(time_V1_wltp, V1_fc_wltp, '-b', time_V2_wltp, V2_fc_wltp, '-r', 'LineWidth', 1.5)
yyaxis right
plot(time_V1_wltp, V1_fc_tot_wltp, '-y', time_V2_wltp, V2_fc_tot_wltp, '-g', 'LineWidth', 1.5)
grid on
xlabel("Time [s]")
ylabel("Fuel Consumption [kg/h]")
legend("Idea Instantaneous", "Punto Instantaneous", "Idea Cumulative", ...
    "Punto Cumulative")
title("Fuel Consumption vs WLTP Time")

% NOx Instantaneous Emissions vs Time
figure(6)
subplot(2, 1 ,1)
yyaxis left
plot(time_NEDC, V1_nox_nedc, '-b', time_NEDC, V2_nox_nedc, '-r', 'LineWidth', 1.5)
yyaxis right
plot(time_NEDC, V1_nox_tot_nedc, '-y', time_NEDC, V2_nox_tot_nedc, '-g', 'LineWidth', 1.5')
grid on
xlabel("Time [s]")
ylabel("NOx [g/h]")
legend("Idea Instantaneous", "Punto Instantaneous", "Idea Cumulative", ...
    "Punto Cumulative")
title("NOx vs NEDC Time")

subplot(2, 1, 2)
yyaxis left
plot(time_V1_wltp, V1_nox_wltp, '-b', time_V2_wltp, V2_nox_wltp, '-r', 'LineWidth', 1.5)
yyaxis right
plot(time_V1_wltp, V1_nox_tot_wltp, '-y', time_V2_wltp, V2_nox_tot_wltp, '-g', 'LineWidth', 1.5)
grid on
xlabel("Time [s]")
ylabel("NOx [g/h]")
legend("Idea Instantaneous", "Punto Instantaneous", "Idea Cumulative", ...
    "Punto Cumulative")
title("NOx vs WLTP Time")

%% Operating Map Creation
% Power corresponsing to fuel consumption data
P = bmep_pq*rpm_pq*V1_engine_displacement/1200; % [kW]
bsfc = fc_pq./P*1000; % [g/kWh]

%% Limiting BSFC values to the WOT curve
bmep_wot = interp1(rpm_max, bmep_max, rpm_pq, 'linear', 'extrap');
bmep_mtx = repmat(bmep_wot, length(bmep_pq),1);
bmep_pq_mtx = repmat(bmep_pq, 1, length(rpm_pq));
bsfc(bmep_pq_mtx > bmep_mtx) = NaN;

%% Creation of BSFC Contour
V = [200 210 220 230 240 250 275 300 325 350 375 400];
figure(7), [c,h] = contour(rpm_pq, bmep_pq, bsfc, V, 'LineWidth', 3);
clabel(c, h)
xlabel("Engine Speed [rpm]")
ylabel("bmep [bar]")
hold all
plot(rpm_max, bmep_max, '-r', 'LineWidth', 3)
plot(V1_rpm, V1_bmep_nedc, 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 1)
plot(V2_rpm, V2_bmep_nedc, 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1)
hold off
title("Engine Map NEDC");
legend("BSFC", "WOT", "Idea OP", "Punto OP")

figure(8), [c,h] = contour(rpm_pq, bmep_pq, bsfc, V, "LineWidth", 3);
clabel(c, h)
xlabel("Engine Speed [rpm]")
ylabel("bmep [bar]")
hold all
plot(rpm_max, bmep_max, "--r", "LineWidth", 3)
plot(V1_rpm_wltp, V1_bmep_wltp, 'o', 'MarkerEdgeColor', 'r', "LineWidth", 1)
plot(V2_rpm_wltp, V2_bmep_wltp, 'o', 'MarkerEdgeColor', 'b', "LineWidth", 1)
title("Engine Map WLTP")
legend("BSFC", "WOT", "Idea OP", "Punto OP")
hold off