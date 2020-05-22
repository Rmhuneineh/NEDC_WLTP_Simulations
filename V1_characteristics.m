%% Idea
%% Gears
V1_gears = 5; % [-]
V1_finalDrive_ratio = 3.563; % [-]
V1_tau_g = [3.909; 2.238; 1.444; 1.029; 0.767]; % [-]
V1_finalDrive_eff = 1; % [-]
V1_gear_eff = 0.94; % [-]

%% Inertia
V1_engine_inertia = 0.183; % [kg.m^2]
V1_wheels_inertia = 2.7794; % [kg.m^2]

%% NEDC
V1_NEDC_mass = 1168; % [kg]
V1_NEDC_F0 = 124.7; % [N]
V1_NEDC_F1 = 0; % [N/(km/h)]
V1_NEDC_F2 = 0.0364; % [N/(km/h)^2]

%% WLTP TMH
V1_WLTP_mass = 1360; % [kg]
V1_WLTP_F0 = 186; % [N]
V1_WLTP_F1 = 0; % [N/(km/h)]
V1_WLTP_F2 = 0.0419; % [N/(km/h)^2]

%% Other
V1_fuel_const_idle = 315; % [g/h]
V1_omega_min = 800; % [rpm]
V1_fuel_density = 835; % [g/l]
V1_wheel_radius = 289.3436865*1e-3; % [m]
V1_engine_displacement = 1.248; % [dm^3]