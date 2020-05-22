%% Punto
%% Gears
V2_gears = 5; % [-]
V2_finalDrive_ratio = 3.563; % [-]
V2_tau_g = [3.909; 2.238; 1.444; 1.029; 0.767]; % [-]
V2_finalDrive_eff = 1; % [-]
V2_gear_eff = 0.94; % [-]

%% Inertia
V2_engine_inertia = 0.183; % [kg.m^2]
V2_wheels_inertia = 2.7794; % [kg.m^2]

%% NEDC
V2_NEDC_mass = 1063; % [kg]
V2_NEDC_F0 = 114.2; % [N]
V2_NEDC_F1 = 0; % [N/(km/h)]
V2_NEDC_F2 = 0.0344; % [N/(km/h)^2]

%% WLTP TMH
V2_WLTP_mass = 1210; % [kg]
V2_WLTP_F0 = 166; % [N]
V2_WLTP_F1 = 0; % [N/(km/h)]
V2_WLTP_F2 = 0.039; % [N/(km/h)^2]

%% Other
V2_fuel_const_idle = 315; % [g/h]
V2_omega_min = 800; % [rpm]
V2_fuel_density = 835; % [g/l]
V2_wheel_radius = 285.3648*1e-3; % [m]
V2_engine_displacement = 1.248; % [dm^3]