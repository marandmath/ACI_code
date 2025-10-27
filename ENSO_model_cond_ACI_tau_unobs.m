% Code corresponding to the ENSO model in Sections 3.2 and SI.4.3 - "A 
% Stochastic Model Capturing the El Niño-Southern Oscillation (ENSO) Diversity" 
% of the paper "Assimilative Causal Inference".
%
% Authors: Marios Andreou, Nan Chen, Erik Bollt.
%
% Code Information: Application of conditional Assimilative Causal Inference 
% (ACI) to a six-dimensional stochastic dynamical system that captures the 
% diverse behaviors and non-Gaussian statistics of ENSO, known as ENSO diversity 
% or complexity:
%
%       du/dt   = - ru - δ_u(T_C + T_E)/2 + β_u(I)τ + σ_u·dot(W)_u,
%       dh_W/dt = - rh_W - δ_h(T_C + T_E)/2 + β_h(I)τ + σ_h·dot(W)_h,
%       dT_C/dt = (r_C - c_1(t,T_C))T_C + ζ_CT_E + γ_Ch_W + σ(I)u + C_u + β_C(I)τ + σ_C·dot(W)_C,
%       dT_E/dt = (r_E - c_2(t))T_E - ζ_ET_C + γ_Eh_W + β_E(I)τ + σ_E·dot(W)_E,
%       dτ/dt   = -d_ττ + σ_τ(t,T_C)·dot(W)_τ,
%       dI/dt   = -λ(I-m) + σ_I(I)·dot(W)_I.
%
% This is a reduced-order conceptual model for capturing the ENSO diversity with
% temporal multiscale features and a three-region resolution of the equatorial 
% Pacific ocean (averaged over 5°S–5°N): Western Pacific (WP), Central Pacific 
% (CP), and Eastern Pacific (EP). It was first developed by Chen, Feng and Yu: 
%   10.1038/s41612-022-00241-x
% It is non-dimensionalized. The state variables are the following:
%
%   • u:   Ocean zonal current in the CP (interannual climatological anomaly)
%   • h_W: Thermocline depth in the WP (interannual climatological anomaly)
%   • T_C: Sea surface temperature (SST) in the CP (interannual climatological 
%          anomaly)
%   • T_E: SST in the EP (interannual climatological anomaly)
%   • τ:   Atmospheric westerly wind burst amplitude, including the MJO
%          (intraseasonal)
%   • I:   Background Walker circulation and the zonal SST difference between 
%          the WP and CP regions that directly determines the strength of the 
%          zonal advective feedback, σ(I) (decadal)
%
% All the prognostic state variables represent the deviations from their 
% corresponding monthly climatology during the analysis period. The two SST 
% variables, T_C (Niño4; 5°S–5°N, 160°E–150°W) and T_E (Niño3; 5°S–5°N, 
% 150°W–90°W), allow reconstruction of spatiotemporal SST patterns across the 
% equatorial Pacific, specifically the formation of an SST anomaly index for the
% Niño3.4 region (T_C+T_E; 5°S–5°N, 170°W–120°W) through a bivariate linear 
% regression model (10.1038/s41612-022-00241-x; lon_sst.mat & 
% nino3.4ssta_regression_constant_T_E_T_C_h_W_coeffs.mat). Niño3.4 represents 
% the average of the equatorial SSTs across the Pacific from about the dateline 
% to the South American coast. ENSO exhibits remarkable diversity in its spatial 
% patterns, temporal evolution, and impacts, which can be broadly 
% categorized into two main events: EP and CP El Niños, where an anomalous 
% warming center occurs in the eastern and central Pacific, respectively. The 
% opposite phases with anomalous cooling SSTs are called La Niña. The way these 
% different ENSO events are rigorously defined in this study/script are as 
% follows:
%
%   The definitions of the different ENSO events are based on the average SST 
%   anomalies (SSTas) during boreal winter (December-January-February; DJF). 
%   Using the definitions of Kug et al. (2009; 10.1175/2008JCLI2624.1), when the 
%   EP is warmer than the CP and the EP SSTa T_E is greater than 0.5°C, it is 
%   classified as an EP El Niño event. Based on the classification in Wang et 
%   al. (2019; 10.1073/pnas.1911130116), an extreme EP El Niño event corresponds 
%   to when the maximum of the EP SSTas from April to the following March is 
%   larger than 2.5°C. Accordingly, when the CP is warmer than the EP and the CP 
%   SSTa T_C is larger than 0.5°C, it is defined as a CP El Niño event. Finally, 
%   when either the T_C or T_E SSTas are cooler than −0.5°C, it is defined as a 
%   La Niña event.
%
% This model is a conditional Gaussian nonlinear system (CGNS) in u, h_W, T_E
% and τ, i.e., when these variables are considered to be part of the unobserved
% process; See Section 2.1 of the Supplementary Information. It can also be 
% considered as a CGNS in T_C, by approximating the nonlinear damping factor 
% c_1(t,T_C) around the climatology T_C = 0 through a zeroth-order accurate 
% Taylor approximation, thus transforming the feedback of T_C in T_C to be state 
% independent, r_C-c_1(t,0), which achieves conditional linearity in T_C.
% 
% In this script, the conditional ACI framework is employed to study the
% conditional causal relationship τ(t) → (T_C,T_E,I) | (u,h_W) over time 
% t∈[0,T], as well as any of its other conditional variants, i.e., where we
% additionally resolve the observational contributions of any subset of one or 
% two of the observed variables (T_C,T_E,I) from this conditional causal 
% relationship through the conditional ACI framework; See Section 1.4 of the 
% Supplementary Information. Therefore, in this script we consider 
% (u,h_W,T_C,T_E,I) to be the observables, while τ is the unobserved variable. 
% This code uses the same parameter values as those cited in the following 
% paper: 
%   10.1175/JCLI-D-24-0017.1
%
% Written and tested in MATLAB R2024b.
%
% MATLAB Toolbox and M-file Requirements:
% 
% Code used to obtain the required m-file scripts and MATLAB toolboxes:
% [fList, pList] = matlab.codetools.requiredFilesAndProducts('ENSO_model_cond_ACI_u_h_W_tau_unobs.m');
%
% M-file Scripts:
%   ➤ ENSO_model_cond_ACI_u_h_W_tau_unobs.m
%   ➤ progress_bar.m
%   ➤ legendUnq.m (https://www.mathworks.com/matlabcentral/fileexchange/67646-legendunq)
%
% Data: The user has two options for which observations to use in this script or
% ACI analysis: Either use synthetic data generated from the ENSO model 
% simulation, or instead assimilate real data. The user can simply load the 
% observational data as needed from the datasets provided in the "ENSO_DATA" 
% subdirectory and instead use these observations in the filter, smoother, and 
% online smoother algorithms required for the ACI analysis. This script takes a 
% more direct approach and for its ACI analysis it instead assimilates the 
% values generated for the state variables u, h_W, T_C, T_E, τ, I from running 
% the six-dimensional model forward via the Milstein numerical integration
% scheme implemented in this script. A detailed description of how the real data
% are obtained and curated is provided below:
% 
%   The monthly ocean temperature and flow data are from the NCEP Global Ocean 
%   Data Assimilation System reanalysis dataset (GODAS; Behringer and Xue 2004):
%       https://www.esrl.noaa.gov/psd/data/gridded/data.godas.html
%   The thermocline depth along the equatorial Pacific is approximated from the 
%   potential temperature as the depth of the 20°C isotherm contour. The 
%   analysis period of the GODAS dataset that is included in this repository is 
%   40 years: 01/1980–12/2019. The anomalies are calculated by removing the 
%   monthly mean climatology over the whole period. The T_C (Niño4) and T_E 
%   (Niño3) are the averages of the SSTas over the CP and EP regions, 
%   respectively. The h_W index is the mean thermocline depth anomaly over the 
%   WP region (5°S–5°N, 120°E–180°), while the u index is the mean mixed-layer 
%   zonal current in the CP region. If needed by the user, the NOAA Extended 
%   Reconstructed Sea Surface Temperature Version 5 dataset (ERSST.v5; Huang 
%   et al. 2017) can be used instead:
%       https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
%   The analysis period of the ERSST.v5 dataset that is included in this 
%   repository is almost 151 years: 01/1870–10/2020. Finally, the daily zonal 
%   wind data are measured at 850 hPa and are taken from the NCEP–NCAR 
%   reanalysis 1 project (Kalnay et al. 1996):
%       https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html
%   It is used to describe the wind bursts in the intraseasonal scale by 
%   removing the daily mean climatology and averaging the anomalies over the WP
%   to create the wind burst index. The wind process evolves on a faster time 
%   scale than all other variables (daily than monthly), and although a single 
%   daily value of τ describing the wind anomalies has a minor effect on the 
%   SSTas, the accumulated wind effect over time will modulate the SST
%   variations. Finally, the Walker circulation strength index, I, data are 
%   included to illustrate the modulation of the decadal variation on the 
%   interannual ENSO characters. It is defined as the sea level pressure 
%   difference over the CP/EP (5°S–5°N, 160°W–80°W) and the Indian Ocean/WP 
%   (5°S–5°N, 80°E–160°E) (Kang et al. 2020; 10.1126/sciadv.abd3021). The 
%   monthly zonal SST gradient between the WP and CP regions highly correlates 
%   with this Walker circulation strength index (with a simultaneous correlation 
%   coefficient of ≈0.85), suggesting the significance of the air–sea 
%   interactions over the equatorial Pacific. Since the latter is more directly 
%   related to the zonal advective feedback strength σ(I) over the CP region, 
%   the decadal state variable I mainly illustrates this quantity. 
% 
%   The dataset files included in this repository, in either .mat or .txt form, 
%   are found in the "ENSO_DATA" subdirectory:
%
%   ➤ u (MONTHLY DATA; m/s): 
%       • GODAS_u_h_W_T_C_T_E.mat\u_obs (01/1980–12/2018)
%   ➤ h_W (MONTHLY DATA; m): 
%       • GODAS_u_h_W_T_C_T_E.mat\h_W_obs (01/1980–12/2019)
%   ➤ T_C (Niño4; MONTHLY DATA; °C): 
%       • ERSST.V5_nino4.txt (01/1870–10/2020)
%       • GODAS_u_h_W_T_C_T_E.mat\T_C_obs (01/1950–12/2019)
%   ➤ T_E (Niño3; MONTHLY DATA; °C): 
%       • ERSST.V5_nino3.txt (01/1870–10/2020)
%       • GODAS_u_h_W_T_C_T_E.mat\T_E_obs (01/1950–12/2019)
%   ➤ τ (DAILY DATA; m/s): 
%       • NCEP-NCAR_tau.mat\tau_obs (01/01/1982–29/02/2020)   
%   ➤ I (MONTHLY DATA):
%       • GODAS_I.mat\I_obs (01/1980–12/2017)
%
%   Note: The most recent versions of these datasets can be downloaded from the 
%   links provided above.
%
% Toolboxes:
%   ➤ N/A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  LOADING THE OBSERVATIONAL DATA  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The user should uncomment and run the code in this section if they would like 
% to use real data in their ACI analysis:
%   ➤ u (MONTHLY DATA; m/s): 
%       • GODAS_u_h_W_T_C_T_E.mat\u_obs (01/1980–12/2018)
%   ➤ h_W (MONTHLY DATA; m): 
%       • GODAS_u_h_W_T_C_T_E.mat\h_W_obs (01/1980–12/2019)
%   ➤ T_C (Niño4; MONTHLY DATA; °C): 
%       • ERSST.V5_nino4.txt (01/1870–10/2020)
%       • GODAS_u_h_W_T_C_T_E.mat\T_C_obs (01/1950–12/2019)
%   ➤ T_E (Niño3; MONTHLY DATA; °C): 
%       • ERSST.V5_nino3.txt (01/1870–10/2020)
%       • GODAS_u_h_W_T_C_T_E.mat\T_E_obs (01/1950–12/2019)
%   ➤ τ (DAILY DATA; m/s): 
%       • NCEP-NCAR_tau.mat\tau_obs (01/01/1982–29/02/2020)   
%   ➤ I (MONTHLY DATA):
%       • GODAS_I.mat\I_obs (01/1980–12/2017)

% % Overlapping observation years (among all 6 of the state variables).
% start_year = 1982;
% end_year = 2017;
% 
% % Loading the u, h_W, T_C, T_E, and I GODAS observations into the workspace.
% load('ENSO_DATA\\GODAS_u_h_W_T_C_T_E.mat')
% load('ENSO_DATA\\GODAS_I.mat', 'I_obs')
% % Loading the smoothed version of the westerly wind burst τ NCEP-NCAR data, 
% % obtained from a moving average of the true data using a window width of one 
% % month instead of the true observations tau_obs.
% load('ENSO_DATA\\NCEP-NCAR_tau.mat', 'tau_smoothed', 'tau_day')
% 
% % Loading the ERSST.V5 data into the workspace.
% ERSST_V5_data_year_start = 1950;
% data3 = importdata('ENSO_DATA\\ERSST.V5_nino3.txt');
% nino3 = reshape(data3(:, 2:end)', [], 1);
% nino3 = nino3((ERSST_V5_data_year_start-1870)*12+1:end-12).';
% data4 = importdata('ENSO_DATA\\ERSST.V5_nino4.txt');
% nino4 = reshape(data4(:, 2:end)', [], 1);
% nino4 = nino4((ERSST_V5_data_year_start-1870)*12+1:end-12).';
% data34 = importdata('ENSO_DATA\\ERSST.V5_nino3.4.txt');
% nino34 = reshape(data34(:, 2:end)', [], 1);
% nino34 = nino34((ERSST_V5_data_year_start-1870)*12+1:end-12).';
% 
% % Defining the temporal monthly grid over the intersecting observational period.
% t_obs = start_year:1/12:end_year+11/12;
% obs_period_years = end_year-start_year+1;
% obs_period_months = obs_period_years*12;
% 
% u_obs_start_year = 1980;
% u_obs_end_year = 2018;
% h_W_obs_start_year = 1980;
% h_W_obs_end_year = 2019;
% T_C_obs_start_year = 1950;
% T_C_obs_end_year = 2019;
% T_E_obs_start_year = 1950;
% T_E_obs_end_year = 2019;
% tau_obs_start_year = 1982;
% tau_obs_end_year = 2020;
% I_obs_start_year = 1980;
% I_obs_end_year = 2017;
% 
% % Extracting the mean monthly data (climatology) from the daily westerly wind 
% % burst smoothed observation data.
% tau_obs = tau_smoothed; 
% tau_t = tau_day;
% tau_obs_daily = tau_obs(tau_t >= tau_obs_start_year & tau_t < end_year+1);
% tau_obs_daily = reshape(tau_obs_daily, 365, []);
% daysInMonths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
% monthEnds = [0, cumsum(daysInMonths)];
% tau_obs = zeros(obs_period_years, 12);
% for j = 1:obs_period_years
%     for m = 1 : length(monthEnds)-1
%       firstDay = monthEnds(m)+1;
%       lastDay = monthEnds(m+1);
%       tau_obs(j, m) = mean(tau_obs_daily(firstDay:lastDay, j));
%     end
% end
% tau_obs = reshape(tau_obs.', 1, []);
% 
% % Scaling the observational data to make it dimensionless, since the model 
% % equations and state variables are non-dimensional. Characteristic scales are
% % in the original paper:
% %   10.1038/s41612-022-00241-x
% u_obs_ndim = u_obs(1+(start_year-u_obs_start_year)*12:end-(u_obs_end_year-end_year)*12) / 1.5;
% h_W_obs_ndim = h_W_obs(1+(start_year-h_W_obs_start_year)*12:end-(h_W_obs_end_year-end_year)*12) / 150;
% T_C_obs_ndim = T_C_obs(1+(start_year-T_C_obs_start_year)*12:end-(T_C_obs_end_year-end_year)*12) / 7.5;
% T_E_obs_ndim = T_E_obs(1+(start_year-T_E_obs_start_year)*12:end-(T_E_obs_end_year-end_year)*12) / 7.5;
% tau_obs_ndim = tau_obs / 5;
% I_obs_ndim = I_obs(1+(start_year-I_obs_start_year)*12:end-(I_obs_end_year-end_year)*12);
% 
% % Creating the matrices needed to identify the specific historical ENSO events 
% % which can then be used to draw these events in the associated plots.
% % Extreme EP El Niño events.
% eep_events = zeros(2, obs_period_years);
% % (Moderate) EP El Niño events.
% ep_events = zeros(2, obs_period_years);
% % CP El Niño events.
% cp_events = zeros(2, obs_period_years);
% % La Niña events.
% ln_events = zeros(2, obs_period_years);
% 
% % Checking the first year for extreme EP El Niño manually.
% if max(T_E_obs_ndim(4:12+3) * 7.5) >= 2.5
%    eep_events(:, 1) = [t_obs(1); t_obs(12)]; 
% end
% 
% for year_obs = 1:obs_period_years-1
% 
%     if year_obs <= obs_period_years-2
% 
%         if max(T_E_obs_ndim(year_obs*12+4:(year_obs+1)*12+3) * 7.5) >= 2.5
%             eep_events(:, year_obs+1) = [t_obs(1+year_obs*12); t_obs((year_obs+1)*12)];
%         elseif mean(T_E_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) > 0.5 || mean(T_C_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) > 0.5
%             if mean(T_E_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) > mean(T_C_obs_ndim(year_obs*12:year_obs*12+2) * 7.5)
%                ep_events(:, year_obs+1) = [t_obs(1+year_obs*12); t_obs((year_obs+1)*12)];
%             else    
%                cp_events(:, year_obs+1) = [t_obs(1+year_obs*12); t_obs((year_obs+1)*12)];
%             end
%         elseif mean(T_E_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) < -0.5 || mean(T_C_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) < -0.5
%              ln_events(:, year_obs+1) = [t_obs(1+year_obs*12); t_obs((year_obs+1)*12)];
%         end
% 
%     else
% 
%         if mean(T_E_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) > 0.5 || mean(T_C_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) > 0.5
%             if mean(T_E_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) > mean(T_C_obs_ndim(year_obs*12:year_obs*12+2) * 7.5)
%                ep_events(:, year_obs+1) = [t_obs(1+year_obs*12); t_obs((year_obs+1)*12)];
%             else    
%                cp_events(:, year_obs+1) = [t_obs(1+year_obs*12); t_obs((year_obs+1)*12)];
%             end
%         elseif mean(T_E_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) < -0.5 || mean(T_C_obs_ndim(year_obs*12:year_obs*12+2) * 7.5) < -0.5
%              ln_events(:, year_obs+1) = [t_obs(1+year_obs*12); t_obs((year_obs+1)*12)];
% 
%         end
% 
%     end
% 
% end
% 
% % 2017 was a La Niña year, and as such we denote 2017 in the associated plots 
% % with a blue indicator  manually. This is hard-coded in this script for 
% % convenience, since we do not import data for 2018 to be able to do this
% % automatically in the for-loop above (like all other historical ENSO events).
% ln_events(:, obs_period_years) = [t_obs(end-11); t_obs(end)];
% 
% % These flags are used to appropriately color the ENSO events in the associated 
% % plots:
% %   • Magenta = (Moderate) EP El Niño event
% %   • Orange = CP El Niño event
% %   • Red = Extreme EP El Niño event
% %   • Blue = La Niña event
% events_timeline = zeros(1,obs_period_months);
% eep_flags = any(eep_events, 1);
% ep_flags = any(ep_events, 1);
% cp_flags = any(cp_events, 1);
% ln_flags = any(ln_events, 1);
% for year_obs = 1:obs_period_years
% 
%     if eep_flags(year_obs)
%         events_timeline(1+(year_obs-1)*12:year_obs*12) = 1;
%     end
%     if ep_flags(year_obs)
%         events_timeline(1+(year_obs-1)*12:year_obs*12) = 2;
%     end
%     if cp_flags(year_obs)
%         events_timeline(1+(year_obs-1)*12:year_obs*12) = 3;
%     end
%     if ln_flags(year_obs)
%         events_timeline(1+(year_obs-1)*12:year_obs*12) = 4;
%     end
% 
% end
% 
% % Removing unnecessary zero columns from the arrays used to create the
% % indicators for the ENSO events.
% eep_events(:, ~any(eep_events, 1)) = [];
% ep_events(:, ~any(ep_events, 1)) = [];
% cp_events(:, ~any(cp_events, 1)) = [];
% ln_events(:, ~any(ln_events, 1)) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODEL SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fixing the random number seed for reproducibility across each simulation.
rng(12)

% Numerical integration time step ≈ 7.3hrs = 438 mins, since the characteristic
% time scale of the model is 2 months.
dt = 0.005;
% Total simulation time = 2*T months = 15000*2 = 30000 months = 2500 years.
T = 15000;
% Total number of time steps/observations within the given time interval.
N = round(T/dt);

% Window length needed for calculating the monthly mean data. Characteristic 
% time scale is [t] = 2 months (10.1038/s41612-022-00241-x), so k_dt=1/([t]dt).
% This is essentially the number of array indices spanning a single month.
k_dt = 0.5/dt;

% Model parameters used in the deterministic three-region skeleton system; See
% the following for details:
%   • 10.1038/s41612-022-00241-x
%   • 10.1175/JCLI-D-24-0017.1

% ENSO diversity/complexity modulating parameter, used to modulate the # of 
% extreme EP El Niños, # multi-year La Niñas > # multi-year El Niños, # EP El 
% Niños > # CP El Niños, among other factors.
factor = 0.65;

% Classical recharge-oscillator model's non-dimensionalized deterministic
% parameters.

% High-end estimation of the thermocline tilt, which is in balance with the
% zonal wind stress produced by the SSTas.
b_0 = 2.5;
%  Relative coupling coefficient.
mu = 0.5;  
% These parameters ensure that for a given steady zonal wind stress forcing, the 
% zonal mean thermocline depth anomaly of the recharge-oscillator model is about 
% zero at the equilibrium state, i.e., h_W+h_E=0 (this shows why α_2 is r/2 and 
% α_1 is α_2/2 = r/4).
alpha_2 = 0.125 * factor; 
alpha_1 = alpha_2/2 * factor;

% Average SSTa feedback in u.
delta_u = alpha_1 * b_0 * mu; 
% Average SSTa feedback in h_W.
delta_h = alpha_2 * b_0 * mu; 
% Ocean zonal advection's and thermocline's damping coefficient (collective
% damping rate in the ocean adjustment).
r = 0.25 * factor; 
% Thermocline feedback in T_C.
gamma_C = 0.75 * factor;
% Thermocline feedback in T_E.
gamma_E = 0.75 * factor; 
% T_C's mean damping coefficient.
r_C = gamma_C*b_0*mu/2; 
% T_E's mean damping coefficient.
r_E = 3*gamma_E*b_0*mu/2;
% T_E feedback in T_C.
zeta_C = gamma_C*b_0*mu/2;
% T_C feedback in T_E.
zeta_E = gamma_E*b_0*mu/2;
% A-posteriori calculated term C_u in T_C's stochastic differential equation 
% used to enforce zero climatology in the anomaly model.
C_u = 0.03  * factor; 
% Wind burst's damping coefficient.
d_tau = 2;
% Walker circulation's damping coefficient.
lambda = 2/60; 
% Walker circulation's target equilibrium mean (forcing coefficient).
m = 2; 

% Additive noise feedbacks.
sigma_u = 0.04 * sqrt(factor);
sigma_h = 0.02 * sqrt(factor);
sigma_C = 0.04 * sqrt(factor);
% In the original model, T_E is actually a random ordinary differential 
% equation, i.e., σ_E is set equal to zero. Per the requirements of conditional 
% ACI, when considering T_E to be part of the observables it is necessary to 
% include some additive noise in its evolution equation, akin to noise 
% inflation, as to have a well-defined inverse for the observational noise 
% feedback matrix; See Sections 1.2.1 and 1.4 of the Supplementary Information. 
% As to not impair the baseline model's skill in capturing the ENSO diversity, 
% we empirically set the smallest possible additive noise feedback σ_Ε while 
% still maintaining the stability of the filter and (online) smoother algorithms 
% for the underlying CGNS.
sigma_E = sqrt(5) * 1e-2 * sqrt(factor);
% sigma_E = 0;

% Ocean zonal current in the CP.
u = zeros(1, N+1);
% Thermocline depth in the WP.
h_W = zeros(1, N+1);
% SST in the CP.
T_C = zeros(1, N+1);
% SST in the EP.
T_E = zeros(1, N+1);
% Westerly wind burst amplitude.
tau = zeros(1, N+1);
% Background Walker circulation strength.
I = zeros(1, N+1);

% These ICs are the averages from the GODAS data between 01/1982–12/2017. Real 
% observational data from the datasets in the "ENSO_DATA" subdirectory can be 
% used if preferred by the user.
u(1) = 6.9136e-04;
h_W(1) = -0.0028;
T_C(1) = 0.0039;
T_E(1) = 0.0051;
tau(1) = -0.0256;
I(1) = 1.5841;

% IN THIS SCRIPT:
%   • Observables: x = (u,h_W,T_C,T_E,I)
%   • Unobserved variables: y = τ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NOTE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN THE FOLLOWING, WE CAN EITHER CONSIDER u AND h_W AS PART OF THE OBSERVABLES, 
% AS ALREADY OUTLINED, OR WE CAN INSTEAD APPLY CONDITIONAL ACI BY PROXY VIA AN 
% OPERATIONAL SHORTCUT, AS DESCRIBED IN SECTION 2.3 OF THE MAIN TEXT:
%   "For dynamical systems with explicit governing equations, a shortcut to 
%    exclude the influence in the uncertainty of x_B is to treat x_B as a 
%    prescribed [deterministic] forcing term in the reduced system defined by 
%    (x_A, y) during Bayesian inference, with values defined by its observed 
%    time series[.]"
% THEREFORE, IN THE FOLLOWING IT IS EQUIVALENT TO SIMPLY ASSUME x = (T_C,T_E,I)
% AND HAVE u AND h_W TO BE DETERMINISTIC FORCINGS DEFINED BY THEIR OBSERVED 
% VALUES (IN THIS CASE THOSE GENERATED BY THE MODEL SIMULATION) DURING THE 
% BAYESIAN UPDATE WHEN DERIVING THE POSTERIOR FILTER AND SMOOTHER DISTRIBUTIONS; 
% THIS STRATEGY IS MADE SIMPLE BY THE FACT THERE IS ALSO NO NOISE 
% CROSS-INTERACTING TERMS TO CONSIDER FOR THIS MODEL. NOTE THAT BOTH APPROACHES, 
% EITHER THIS ONE OR THE ONE WITH x = (u,T_C,T_E,τ,I), WILL YIELD THE SAME 
% RESULTS DURING CONDITIONAL ACI ANALYSIS, BUT THE FORMER IS MUCH MORE 
% COMPUTATIONALLY EFFICIENT.

% Observable coefficient matrix: Feedback of y in x.
L_x = zeros(3, N+1);
% Forcing in the observable process.
f_x = zeros(3, N+1);
% Noise feedback matrix in the observable process.
S_x = zeros(3, 3, N+1);
% NOISE CROSS-INTERACTION TERMS ARE ABSENT FROM THIS MODEL SO WE DO NOT DEFINE 
% THE ASSOCIATED MATRIX-VALUED FUNCTIONALS IN THIS SCRIPT FOR SIMPLICITY.

% Unobservable coefficient matrix: Feedback of y in y.
L_y = -d_tau;
% Forcing in the unobservable process.
f_y = 0;
% Noise feedback matrix in the unobservable process.
S_y = zeros(1, N+1);
% NOISE CROSS-INTERACTION TERMS ARE ABSENT FROM THIS MODEL SO WE DO NOT DEFINE 
% THE ASSOCIATED MATRIX-VALUED FUNCTIONALS IN THIS SCRIPT FOR SIMPLICITY.

L_x(:, 1) = [0.8 * (1 + (1 - I(1)/5))*0.15 * sqrt(factor); 1 * (1 + (1 - I(1)/5))*0.15 * sqrt(factor); 0];
f_x(:, 1) = [
    (r_C - (25 * (T_C(1) + 0.75/7.5).^2 + 0.9) * (1 + 0.3*sin(0*dt*2*pi/6 - pi/6)) * factor) * T_C(1) + zeta_C * T_E(1) + gamma_C * h_W(1) + I(1)/5 * u(1) * factor + C_u;
    (r_E - 1.4 * factor * (1 + 0.3*sin(0*dt*2*pi/6 + 2*pi/6) + 0.25*sin(2*0*dt*2*pi/6 + 2*pi/6))) * T_E(1) - zeta_E * T_C(1) + gamma_E * h_W(1);
    - lambda * (I(1) - m)
];
S_x(:, :, 1) = [
    sigma_C, 0, 0;
    0, sigma_E, 0;
    0, 0, sqrt(lambda * (4-I(1)) * I(1))
];

S_y(1) = 0.9 * (tanh(7.5*T_C(1)) + 1) * (1 + 0.3*cos(0*dt*2*pi/6 + 2*pi/6));

% Used for the text-based progress bar.
start_time = tic;

for j = 2:N+1

    progress_bar('Simulating the ENSO Model Using a Milstein Scheme', j-1, N, start_time);

    % Milstein scheme for the decadal Walker circulation (I). Since I is an 
    % Ornstein-Uhlenbeck process, with this multiplicative or state-dependent 
    % noise feedback σ_I(I) we have that the target equilibrium distribution of 
    % I is Uniform([0, 4]).
    if I(j-1) > 4 || I(j-1) < 0
        sigma_I = 0;
        der_sigma_I = 0;
        if I(j-1) > 4
            I(j-1) = 4;
        end
        if I(j-1) < 0
            I(j-1) = 0;
        end
    else
        % I's multiplicative noise coefficient.
        sigma_I = sqrt(lambda * (4-I(j-1)) * I(j-1));
        % Derivative of I's multiplicative noise coefficient.
        der_sigma_I = lambda * (2-I(j-1)) / sigma_I; 
    end

    dW_I = randn;
    % This numerical scheme contains the Milstein correction term which is
    % nonzero since the noise coefficient σ_I(I) is state-dependent.
    I(j) = I(j-1) - lambda * (I(j-1) - m) * dt + sigma_I * sqrt(dt) * dW_I + 0.5 * sigma_I * der_sigma_I * (dt * dW_I^2 - dt);

    % Model parameters appearing in the equations of the interannual state 
    % variables.

    % Zonal advection feedback in T_C.
    sigma = I(j-1)/5 * factor; 
    % Quadratic damping which induces a cubic non-linearity in the collective 
    % mean term of T_C's evolution equation through its mean damping 
    % deviation/discharge coefficient.
    c_1 = 25 * (T_C(j-1) + 0.75/7.5).^2 + 0.9;
    c_1 = c_1 * (1 + 0.3*sin((j-1)*dt*2*pi/6 - pi/6)) * factor;
    % T_E's mean deviation/discharge coefficient.
    c_2 = 1.4 * factor * (1 + 0.3*sin((j-1)*dt*2*pi/6 + 2*pi/6) + 0.25*sin(2*(j-1)*dt*2*pi/6 + 2*pi/6));
    % Wind burst feedback in u.
    beta_u = -0.2 * (1 + (1 - I(j-1)/5))*0.15 * sqrt(factor);
    % Wind burst feedback in h_W.
    beta_h = -0.4 * (1 + (1 - I(j-1)/5))*0.15 * sqrt(factor);
    % Wind burst feedback in T_C.
    beta_C = 0.8 * (1 + (1 - I(j-1)/5))*0.15 * sqrt(factor);
    % Wind burst feedback in T_E.
    beta_E = 1 * (1 + (1 - I(j-1)/5))*0.15 * sqrt(factor);
    
    dW_u = randn;
    dW_h = randn;
    dW_C = randn;
    dW_E = randn;

    % Euler-Maruyama numerical integration scheme for the interannual variables 
    % since their noise feedbacks are constant or state-independent.
    u(j) = u(j-1) + (-r * u(j-1) - delta_u * (T_C(j-1) + T_E(j-1)) / 2) * dt + beta_u * tau(j-1) * dt + sigma_u * sqrt(dt) * dW_u;
    h_W(j) = h_W(j-1) + (-r * h_W(j-1) - delta_h * (T_C(j-1) + T_E(j-1)) / 2) * dt + beta_h .* tau(j-1) * dt + sigma_h * sqrt(dt) * dW_h;
    T_C(j) = T_C(j-1) + ((r_C - c_1) * T_C(j-1) + zeta_C * T_E(j-1) + gamma_C * h_W(j-1) + sigma * u(j-1) + C_u) * dt + beta_C * tau(j-1) * dt + sigma_C * sqrt(dt) * dW_C;
    T_E(j) = T_E(j-1) + ((r_E - c_2) * T_E(j-1) - zeta_E * T_C(j-1) + gamma_E * h_W(j-1)) * dt + beta_E * tau(j-1) * dt + sigma_E * sqrt(dt) * dW_E; 
    
    % Milstein scheme for the intraseasonal westerly wind bursts (τ). This 
    % numerical scheme contains the Milstein correction term which is nonzero 
    % since the noise coefficient σ_τ(t,T_C) is state-dependent.
    % τ's multiplicative noise coefficient.
    sigma_tau = 0.9 * (tanh(7.5*T_C(j-1)) + 1) * (1 + 0.3*cos((j-1)*dt*2*pi/6 + 2*pi/6));
    %  Derivative of τ's multiplicative noise coefficient.
    der_sigma_tau = 0.9 * 7.5 * sech(7.5*T_C(j-1))^2 * (1 + 0.3*cos((j-1)*dt*2*pi/6 + 2*pi/6)); %
    dW_tau = randn;
    tau(j) = tau(j-1) - d_tau * tau(j-1) * dt + sigma_tau * sqrt(dt) * dW_tau + 0.5 * sigma_tau * der_sigma_tau * (dt * dW_tau^2 - dt);

    L_x(:, j) = [0.8 * (1 + (1 - I(j)/5))*0.15 * sqrt(factor); 1 * (1 + (1 - I(j)/5))*0.15 * sqrt(factor); 0];   
    f_x(:, j) = [
        (r_C - (25 * (T_C(j) + 0.75/7.5).^2 + 0.9) * (1 + 0.3*sin(j*dt*2*pi/6 - pi/6)) * factor) * T_C(j) + zeta_C * T_E(j) + gamma_C * h_W(j) + I(j)/5 * u(j) * factor + C_u;
        (r_E - 1.4 * factor * (1 + 0.3*sin(j*dt*2*pi/6 + 2*pi/6) + 0.25*sin(2*j*dt*2*pi/6 + 2*pi/6))) * T_E(j) - zeta_E * T_C(j) + gamma_E * h_W(j);
        - lambda * (I(j) - m)
    ];
    S_x(:, :, j) = [
        sigma_C, 0, 0;
        0, sigma_E, 0;
        0, 0, sqrt(lambda * (4-I(j)) * I(j))
    ];
    
    S_y(j) = 0.9 * (tanh(7.5*T_C(j)) + 1) * (1 + 0.3*cos(j*dt*2*pi/6 + 2*pi/6));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  OBTAINING MONTHLY CLIMATOLOGY FROM SIMULATION DATA  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtaining mean monthly (climatological) data from the values of the state 
% variables generated by the model simulation. These are used for plotting 
% purposes but they can also be used for monthly (conditional) ACI analysis.

% Smoothing the response data using movmean() to obtain a moving mean time 
% series which smoothes-out short time-scale variations in the interannual and 
% decadal state variables. This is different from the smooth() and filter() 
% functions (see MATLAB docs); movmean() is preferable for our purposes.
u_monthly = movmean(u, k_dt);
u_monthly = u_monthly(1:k_dt:end);
h_W_monthly = movmean(h_W, k_dt);
h_W_monthly = h_W_monthly(1:k_dt:end);
T_C_monthly = movmean(T_C, k_dt);
T_C_monthly = T_C_monthly(1:k_dt:end);
T_E_monthly = movmean(T_E, k_dt);
T_E_monthly = T_E_monthly(1:k_dt:end);
% Not needed for τ since the intraseasonal wind bursts are already in the daily 
% time scale, which is in a much smaller time scale than the other variables.
% tau_monthly = movmean(tau, k_dt); 
tau_monthly = tau(1:k_dt:end);
I_monthly = movmean(I, k_dt);
I_monthly = I_monthly(1:k_dt:end);

% Transforming the interannual SST state variables to mean anomaly variables 
% (i.e., deviations from their monthly climatology).
u_monthly = u_monthly-mean(u_monthly);
h_W_monthly = h_W_monthly-mean(h_W_monthly);
T_C_monthly = T_C_monthly-mean(T_C_monthly);
T_E_monthly = T_E_monthly-mean(T_E_monthly);

% Rescaling the data based on their characteristic scales (since the model is
% non-dimensionalized).
u_monthly = 1.5 * u_monthly;
h_W_monthly = 150 * h_W_monthly;
T_C_monthly = 7.5 * T_C_monthly;
T_E_monthly = 7.5 * T_E_monthly;
tau_monthly = 5 * tau_monthly;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILTERING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inverse of the observational noise coefficient's Grammian.
S_xoS_x_inv = zeros(3, 3, N+1);
% Calculating the pseudoinverse for stability concerns.
S_xoS_x_inv(:, :, 1) = pinv(S_x(:, :, 1) * S_x(:, :, 1)');

% THE FOLLOWING CODE IS FOR IMPLEMENTING THE CONDITIONAL ACI FRAMEWORK TO THIS 
% CASE STUDY: USE THIS DEFINITION OF S_xoS_x_inv INSTEAD FOR CONDITIONAL ACI.
% Use the following to study h_W(t) → T_C | (u,T_E,τ,I) over time t∈[0,T]:
S_xoS_x_inv(1, 1, :) = 1/sigma_C^2;
% Use the following to study h_W(t) → T_E | (u,T_C,τ,I) over time t∈[0,T]:
%   S_xoS_x_inv(2, 2, :) = 1/sigma_E^2;
% Use the following to study h_W(t) → I | (u,T_C,T_E,τ) over time t∈[0,T]:
%   sigma_I = sqrt(lambda .* (4-I) .* I);
%   S_xoS_x_inv(3, 3, :) = 1./sigma_I.^2;
%   S_xoS_x_inv(3, 3, sigma_I == 0) = 0;

% Grammian of the unobservable process' noise feedback.
S_yoS_y = S_y.^2;
% NOISE CROSS-INTERACTION TERMS ARE ABSENT FROM THIS MODEL SO WE DO NOT DEFINE 
% THE ASSOCIATED MATRIX-VALUED FUNCTIONALS IN THIS SCRIPT FOR SIMPLICITY.

% Posterior filter mean of the latent variable y.
filter_mean = zeros(1, N+1);
% Initial value of the posterior filter mean.
filter_mean(1) = tau(1);
% Posterior filter covariance matrix of the latent variable y.
filter_cov = zeros(1, N+1);
% Initial value of the posterior filter covariance. Choosing a positive definite 
% matrix as to preserve the positive-definiteness of the posterior covariance 
% matrices over time.
filter_cov(1) = 0.1;

% Used for the text-based progress bar.
start_time = tic;

for j = 2:N+1

    progress_bar('Filter Algorithm', j-1, N, start_time);

    dx = [T_C(j)-T_C(j-1); T_E(j)-T_E(j-1); I(j)-I(j-1)];

    % Update the posterior filter mean and posterior filter covariance using the
    % optimal nonlinear filter state estimation equations for CGNSs; See Section 
    % 2.1.2 of the Supplementary Information.
    % NOISE CROSS-INTERACTION TERMS ARE ABSENT FROM THIS MODEL SO WE DO NOT 
    % INCLUDE THEM IN THE STATE ESTIMATION EQUATIONS IN THIS SCRIPT FOR 
    % SIMPLICITY.
    filter_mean(j) = filter_mean(j-1) ...
                     + (L_y * filter_mean(j-1) + f_y) * dt ...
                     + filter_cov(j-1) * L_x(:, j-1)' * S_xoS_x_inv(:, :, j-1) * (dx - (L_x(:, j-1) * filter_mean(j-1) + f_x(:, j-1)) * dt);
    filter_cov(j) = filter_cov(j-1) ...
                    + (L_y * filter_cov(j-1) + filter_cov(j-1) * L_y' + S_yoS_y(j-1)) * dt ...
                    - (filter_cov(j-1) * L_x(:, j-1)' * S_xoS_x_inv(:, :, j-1) * L_x(:, j-1) * filter_cov(j-1)) * dt;

    % IF IMPLEMENTING THE CONDITIONAL ACI FRAMEWORK, COMMENT-OUT THE FOLLOWING
    % LINE.
    % S_xoS_x_inv(:, :, j) = pinv(S_x(:, :, j) * S_x(:, :, j)');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SMOOTHING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Posterior smoother mean of the latent variable y.
smoother_mean = zeros(1, N+1); 
% Posterior smoother covariance matrix of the latent variable y.
smoother_cov = zeros(1, N+1);

% Smoother runs backwards: "Initial" values of the smoother statistics (i.e., at
% the last time instant) are the corresponding posterior filter statistics.
smoother_mean(N+1) = filter_mean(N+1);
smoother_cov(N+1) = filter_cov(N+1);

% Auxiliary matrices used for the calculation of the online smoother for this 
% CGNS. The online smoother is required for the calculation of the objective 
% causal influence range (CIR) length for y(t) → x. Notation used is consistent 
% with that of the original CGNS online smoother work and the accompanying 
% martingale-free introduction to CGNSs paper: 
%   10.48550/arXiv.2411.05870 and 10.48550/arXiv.2410.24056
% SINCE THIS MODEL DOES NOT INCLUDE NOISE CROSS-INTERACTION TERMS, THE ONLINE 
% SMOOTHER AUXILIARY MATRICES SIMPLIFY SIGNIFICANTLY.
E_j_matrices = zeros(1, N+1);
F_j_matrices = zeros(1, 3, N+1);
G_y_j = L_y + S_yoS_y(N+1) / filter_cov(N+1);
C_jj = 1 - G_y_j * dt;
E_j_matrices(N+1) = C_jj;
F_j_matrices(:, :, N+1) = G_y_j * filter_cov(N+1) * L_x(:, N+1)' * S_xoS_x_inv(:, :, N+1) * dt;

muT = smoother_mean(N+1);
RT = smoother_cov(N+1);

% Used for the text-based progress bar.
start_time = tic;

for j = N:-1:1

    progress_bar('Smoother Algorithm', N-j+1, N, start_time);

    % Calculation of the online smoother auxiliary matrices.
    G_y_j = L_y + S_yoS_y(j) / filter_cov(j);
    C_jj = 1 - G_y_j * dt;
    E_j_matrices(j) = C_jj;
    F_j_matrices(:, :, j) = G_y_j * filter_cov(j) * L_x(:, j)' * S_xoS_x_inv(:, :, j) * dt;

    A_j = L_y;
    B_j = S_yoS_y(j);

    % Update the posterior smoother mean and posterior smoother covariance using
    % the optimal nonlinear smoother state estimation backward equations for 
    % CGNSs; See Section 2.1.2 of the Supplementary Information.
    % NOISE CROSS-INTERACTION TERMS ARE ABSENT FROM THIS MODEL SO WE DO NOT 
    % INCLUDE THEM IN THE STATE ESTIMATION EQUATIONS IN THIS SCRIPT FOR 
    % SIMPLICITY.
    mu = muT - (L_y * muT + f_y - B_j / filter_cov(j) * (filter_mean(j) - muT)) * dt;
    R = RT - ((A_j + B_j / filter_cov(j)) * RT + RT * (A_j + B_j / filter_cov(j))' - B_j) * dt;

    smoother_mean(j) = mu;
    smoother_cov(j) = R;
    muT = mu;
    RT = R;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ACI ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NOTE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since the total model simulation period is very long, for the ACI analysis we 
% instead focus on a much smaller period of years, controlled by the parameter 
% ACI_period_years. In this period we always include 4 more years than what we 
% actually want as to accommodate for the fact that we will initiate the online
% smoother CGNS algorithm not from the start of the simulation but instead from 
% a time instant usually much later than that. As a result, to account for this
% approach for computational efficiency, we add two years before the period of 
% interest and two years after it. Therefore, if for example we want to study 
% the period of 20 model simulation years between 01/350–12/369, we should set:
%   • sim_year_start = 348 (2 years before)
%   • ACI_period_years = 24 (4 years more than needed)

sim_year_start = 1980;
ACI_period_years = 14;
sim_month_start = sim_year_start*12*k_dt;
sim_month_end = sim_month_start + ACI_period_years*12*k_dt;

% Shortening all variables of interest in the MATLAB workspace to the period of 
% interest (± the 2 years added as mentioned at the start of this section). We
% add '_s' to the names of these variables to denote this fact.
N_s = ACI_period_years*12*k_dt;
u_s = u(sim_month_start:sim_month_end);
h_W_s = h_W(sim_month_start:sim_month_end);
T_C_s = T_C(sim_month_start:sim_month_end);
T_E_s = T_E(sim_month_start:sim_month_end);
tau_s = tau(sim_month_start:sim_month_end);
I_s = I(sim_month_start:sim_month_end);
L_x_s = L_x(:, sim_month_start:sim_month_end);
f_x_s = f_x(:, sim_month_start:sim_month_end);
S_x_s = S_x(:, :, sim_month_start:sim_month_end);
S_yoS_y_s = S_yoS_y(sim_month_start:sim_month_end);
S_xoS_x_inv_s = S_xoS_x_inv(:, :, sim_month_start:sim_month_end);
E_j_matrices_s = E_j_matrices(sim_month_start:sim_month_end);
F_j_matrices_s = F_j_matrices(:, :, sim_month_start:sim_month_end);

filter_mean_s = filter_mean(sim_month_start:sim_month_end);
filter_cov_s = filter_cov(sim_month_start:sim_month_end);
smoother_mean_s = smoother_mean(sim_month_start:sim_month_end);
smoother_cov_s = smoother_cov(sim_month_start:sim_month_end);

% Calculating the ACI metric for y(t) → x at each time t∈[0,T].
signal_smoother_filter = 0.5 * (smoother_mean_s - filter_mean_s).^2 ./ filter_cov_s;
cov_ratio_smoother_filter = smoother_cov_s ./ filter_cov_s;
dispersion_smoother_filter = 0.5 * (-log(cov_ratio_smoother_filter) + cov_ratio_smoother_filter - 1);
ACI_metric = signal_smoother_filter + dispersion_smoother_filter;

% Implementation of the fixed-lag online smoother for CGNSs: 
%   10.48550/arXiv.2411.05870
% This is required for the calculation of the objective CIR length for y(t) → x
% at each time t∈[0,T]. The fixed-lag parameter is set equal to the total number
% of observations, N = ⌈T/Δt⌉, such that at each time instant the full backward 
% algorithm is carried out. This is because each online smoother distribution,
% pₙ(yʲ|x), is needed for the calculation of the subjective and objective CIRs; 
% See Section 2.3 of the Supplementary Information.
fixed_lag = N_s+1;

% Saving the online smoother mean and covariance matrices, and update matrices 
% in a cell array where each row is another cell array with as many columns as 
% the index of the current row. Using such nested cell arrays efficiently 
% simulates staggered arrays in MATLAB. This approach preserves space in memory 
% without defining unnecessarily large high-order tensors to store the online 
% smoother estimations and update matrices. In these nested cell arrays, the 
% first/parent index corresponds to n∈{j,j+1,...,N}, for the current observation 
% xⁿ, while the second/child index corresponds to j∈{0,1,...,N}, the time 
% instant tⱼ at which we carry out the online smoother state estimation for 
% yʲ=y(tⱼ).
online_fixed_mean = cell(N_s, 1);
online_fixed_cov = cell(N_s, 1);
update_matrices_fixed = cell(N_s-2, 1);
for n = 1:(N_s-2)
    update_matrices_fixed{n} = zeros(1, n+1);
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end
for n = N_s-1:N_s
    online_fixed_mean{n} = zeros(1, n);
    online_fixed_cov{n} = zeros(1, n);
end

% Details of the online smoother algorithm for CGNSs are also briefly reviewed 
% in Section 2.2 of the Supplementary Information.

% Need to do the first two observations manually.

% A single observation (n=1).
online_fixed_mean{1}(1) = filter_mean_s(1);
online_fixed_cov{1}(1) = filter_cov_s(1);

% Two observations (n=2).
online_fixed_mean{2}(2) = filter_mean_s(2);
online_fixed_cov{2}(2) = filter_cov_s(2);
if fixed_lag == 0

    online_fixed_mean{2}(1) = online_fixed_mean{1}(1);
    online_fixed_cov{2}(1) = online_fixed_cov{1}(1);

else

    aux_vec = filter_mean_s(1) ...
              - E_j_matrices_s(1) * ((1 + L_y * dt) * filter_mean_s(1) + f_y_s(1) * dt) ...
              + F_j_matrices_s(:, :, 1) * ([T_C_s(2)-T_C_s(1); T_E_s(2)-T_E_s(1); I_s(2)-I_s(1)] - (L_x_s(:, 1) * filter_mean_s(1) + f_x_s(:, 1)) * dt);
    online_fixed_mean{2}(1) = E_j_matrices_s(1) * filter_mean_s(2) + aux_vec;
    aux_mat = filter_cov_s(1) ...
              - E_j_matrices_s(1) * (1 + L_y * dt) * filter_cov_s(1) ... 
              - F_j_matrices_s(:, :, 1) * L_x_s(:, 1) * filter_cov_s(1) * dt;
    online_fixed_cov{2}(1) = E_j_matrices_s(1) * filter_cov_s(2) * E_j_matrices_s(1)' + aux_mat;

end

% Used for the text-based progress bar.
start_time = tic;

for n = 3:N_s+1

    progress_bar('Online Smoother Algorithm', n-2, length(3:N_s+1), start_time);

    online_fixed_mean{n}(n) = filter_mean_s(n);
    online_fixed_cov{n}(n) = filter_cov_s(n);

    if fixed_lag == 0

        online_fixed_mean{n}(n-1) = online_fixed_mean{n-1}(n-1);
        online_fixed_cov{n}(n-1) = online_fixed_cov{n-1}(n-1);

    else

        aux_vec = filter_mean_s(n-1) ...
                  - E_j_matrices_s(n-1) * ((1 + L_y * dt) * filter_mean_s(n-1) + f_y_s(n-1) * dt) ...
                  + F_j_matrices_s(:, :, n-1) * ([T_C_s(n)-T_C_s(n-1); T_E_s(n)-T_E_s(n-1); I_s(n)-I_s(n-1)] - (L_x_s(:, n-1) * filter_mean_s(n-1) + f_x_s(:, n-1)) * dt);
        online_fixed_mean{n}(n-1) = E_j_matrices_s(n-1) * filter_mean_s(n) + aux_vec;
        aux_mat = filter_cov_s(n-1) ...
                  - E_j_matrices_s(n-1) * (1 + L_y * dt) * filter_cov_s(n-1) ... 
                  - F_j_matrices_s(:, :, n-1) * L_x_s(:, n-1) * filter_cov_s(n-1) * dt;
        online_fixed_cov{n}(n-1) = E_j_matrices_s(n-1) * filter_cov_s(n) * E_j_matrices_s(n-1)' + aux_mat;

    end

    for j = (n-1):-1:1

        if  (1 <= j) && (j <= n-1-fixed_lag)

            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j);

        elseif (n-fixed_lag <= j) && (j <= n-1)

            if j == n-1
                update_matrices_fixed{n-2}(n-1) = 1;
            elseif j == n-2
                update_matrices_fixed{n-2}(n-2) = E_j_matrices_s(n-2);
            else
                update_matrices_fixed{n-2}(j) = update_matrices_fixed{n-3}(j) * E_j_matrices_s(n-2);
            end
            online_mean_inov = online_fixed_mean{n}(n-1) - filter_mean_s(n-1);
            online_fixed_mean{n}(j) = online_fixed_mean{n-1}(j) + update_matrices_fixed{n-2}(j) * online_mean_inov;
            online_cov_inov = online_fixed_cov{n}(n-1) - filter_cov_s(n-1);
            online_fixed_cov{n}(j) = online_fixed_cov{n-1}(j) + update_matrices_fixed{n-2}(j) * online_cov_inov * update_matrices_fixed{n-2}(j)';

        end        
    end
end

% Actual plotting interval for the ACI analysis (excluding the 2 years before 
% and 2 years after the period of interest).
time_start_plot = 2;
time_end_plot = ACI_period_years-2;

% Calculating the subjective CIR length for y(t) → x at each time t∈[0,T] and 
% for various orders O(10⁻ᵏ) of ε values. The associated objective CIR length is 
% also calculated using its computationally efficient underestimating 
% approximation. The theory behind the subjective and objective CIR length is 
% given in Section 1.5 of the Supplementary Information, while their 
% computational details for CGNSs are given in Section 2.3. 

% Letting 10⁻¹⁰ ≤ ε ≤ 10⁰ with a resolution of 128 points.
epsilon_resolution = 128;
lowest_order = -10;
highest_order = 0;
eps_ord_values = flip(linspace(lowest_order, highest_order, epsilon_resolution));

% Calculating the subjective and objective CIRs over the plotting time interval 
% of choice. We add a lookahead tolerance for the lagged observational time: 
%   T'∈[t,time_end_plot+lookahead_tolerance],
% to avoid observational saturation as t approaches time_end_plot.
lookahead_tolerance = 2;
first_idx = round(time_start_plot*12*k_dt);
if time_end_plot+lookahead_tolerance < T
    last_idx = round((time_end_plot+lookahead_tolerance)*12*k_dt);
else
    last_idx = round(time_end_plot*12*k_dt);
end
plot_len = length(first_idx:last_idx);

subjective_CIR = zeros(length(eps_ord_values), plot_len);
approx_objective_CIR = zeros(1, plot_len);
% CIR relative entropy metric δ(T';t) (See Section 1.5 of the Supplementary 
% Information) used to calculate the approximate objective CIR via a time 
% integral over the lagged observational time T' instead of integrating 
% the associated subjective CIR over ε as in the definition to get the exact 
% objective CIR length. Using the notation from Section 2.3 of the Supplementary 
% Information, RE_metric is Pₙʲ, with the rows of RE_metric corresponding to the 
% natural time t (j∈{first_idx,first_idx+1,...,last_idx} index) while the 
% columns correspond to the lagged observational time T' (n∈{j,j+1,...,last_idx} 
% index).
RE_metric = zeros(plot_len, plot_len); 
max_RE_metric = zeros(1, plot_len);

% Used for the text-based progress bar.
start_time = tic;

for eps_idx = 1:length(eps_ord_values)

    epsilon = 10^eps_ord_values(eps_idx);
    
    for j = first_idx:last_idx

        progress_bar('Calculation of the CIRs', j-first_idx+1+(eps_idx-1)*plot_len, length(eps_ord_values)*plot_len, start_time);
        
        % Calculation of the objective CIR length approximation.
        if eps_idx == 1

            % Calculating and storing Pₙʲ over n∈{j,j+1,...,last_idx} for a 
            % fixed j∈{first_idx,first_idx+1,...,last_idx}.
            RE_n = zeros(1, length(j:last_idx));
            for obs = j:last_idx
                cov_ratio = online_fixed_cov{end}(j) / online_fixed_cov{obs}(j);
                RE_n(obs-j+1) = 0.5 * (online_fixed_mean{end}(j) - online_fixed_mean{obs}(j))^2 / online_fixed_cov{obs}(j) ...
                                + 0.5 * (-log(abs(cov_ratio)) + cov_ratio - 1);
            end
            max_RE_metric(j-first_idx+1) = max(RE_n);
            RE_metric(j-first_idx+1, 1:length(RE_n)) = RE_n;
            % If it is essentially zero then do not calculate the CIR and set 
            % equal to 0 instead.
            RE_metric_threshold = 1e-8;
            if max(RE_n) > RE_metric_threshold 
                approx_objective_CIR(j-first_idx+1) = trapz(RE_n)*dt/max(RE_n);
            else
                approx_objective_CIR(j-first_idx+1) = 0;
            end

        end
        
        % Calculation of the subjective CIR length for this ε value.
        RE_n = RE_metric(j-first_idx+1, 1:length(j:last_idx));
        subj_CIR_idx = find(RE_n > epsilon, 1, 'last');
        if isempty(subj_CIR_idx)
            subj_CIR_idx = 0;
        end
        subjective_CIR(eps_idx, j-first_idx+1) = subj_CIR_idx*dt;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  NOTE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% When choosing a sufficiently large epsilon_resolution and a sufficiently large 
% k where 10⁻ᵏ ≤ ε, then the objective CIR at each time tⱼ can be calculated 
% through the definition by averaging the corresponding subjective CIR over ε 
% via a numerical quadrature method using the following command (the flip 
% operations are needed to put the ε (see eps_ord_values variable) interval in 
% ascending order): 
defn_objective_CIR = trapz(10.^flip(eps_ord_values), flipud(subjective_CIR(:, 1:end-lookahead_tolerance*12*k_dt)), 1)./max_RE_metric(1:end-lookahead_tolerance*12*k_dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  PLOTTING ACI ANALYSIS RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time series of the state variables over the period of interest.
figure('WindowState', 'maximized');

subplot(1, 7, 1)
% Modeling the Niño3.4 index using a time-dependent bivariate linear model with 
% T_C and T_E as the explanatory variables:
%   SSTa3.4(t,x) = B_0(x) + B_C(x)T_C(t) + B_E(x)T_E(t).
% x (in lon_sst.mat) is the longitude independent variable for 120°–280° (with 
% 180° being roughly the dateline) while t is time.

% nino3.4ssta_regression_constant_T_E_T_C_h_W_coeffs.mat contains the regression
% coefficients of this linear model over x: The first column is B_0(x), the 
% second column is B_C(x), and the third one is B_E(x). These coefficients were 
% determined via bivariate linear regression using observational data of T_C 
% (Niño4) and T_E (Niño3) at each longitude grid point x over multiple years 
% using the ERSST.V5 dataset. The fourth and final column corresponds to the 
% regression coefficient of h_W(t) from a multivariate linear model which
% includes h_W as an additional explanatory variable, i.e., the fourth column 
% corresponds to B'_h  for:
%   SST3.4(t,x) = B'_0(x) + B'_C(x)T_C(t) + B'_E(x)T_E(t) + B'_h(x)h_W(t).
load('ENSO_DATA\\nino3.4ssta_regression_constant_T_E_T_C_h_W_coeffs.mat')
load('ENSO_DATA\\lon_sst.mat')
range_months = 1+12*(sim_year_start+2):12*(sim_year_start+ACI_period_years-2);
range_years = range_months/12;
range_idxs = range_months(1)*k_dt:range_months(end)*k_dt;
range_sim_years = range_idxs/k_dt/12;
defn_window = 12; 
total_loop = (ACI_period_years-4-1)*12/defn_window;
Hov = zeros(length(lon_sst), length(range_months));
xx = [ones(size(T_C_monthly(range_months)')), T_E_monthly(range_months)', T_C_monthly(range_months)'];
for i = 1:length(lon_sst)
    Hov(i, :) = xx * reg_te_tc_hw_ssta(i, 1:3)';
end
[xx, yy] = meshgrid(range_years, lon_sst);
contourf(yy, xx, Hov, 30, Linestyle='none', DisplayName='Niño3.4')
hold on
colormap('jet')
plot([180, 180], [range_years(1), range_years(end)], 'k--', Linewidth=2, DisplayName='Dateline');
temp_tau = range_idxs(1:round(k_dt/3):end);
plot(180 + 10*tau(temp_tau), range_sim_years(1:round(k_dt/3):end), 'color', [0.2, 0.49, 0.2], Linewidth=1.5, DisplayName='WWB (τ)');
plot(180 + 10*I(temp_tau), range_sim_years(1:round(k_dt/3):end), 'color', [0.72, 0.27, 1], Linewidth=1.5, DisplayName='Walker Cell (I)');
for k = 1:total_loop
    % Extreme EP El Niño (EEP EN).
    if mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > 1.0 && mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window))
        plot([120, 120], [range_years(1-4+k*defn_window), range_years(1+1+defn_window*k)], 'r', Linewidth=10, DisplayName='EEP EN') 
    % Moderate EP El Niño (MEP EN).
    elseif mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > 0.5 && mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window))
        plot([120, 120], [range_years(1-4+defn_window*k), range_years(1+1+defn_window*k)], 'm', Linewidth=10, DisplayName='MEP EN') 
    % CP El Niño (CP EN).
    elseif mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > 0.5 && mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window))
        plot([120, 120], [range_years(1-4+defn_window*k), range_years(1+1+defn_window*k)], 'color', [1, 0.38, 0], Linewidth=10, DisplayName='CP EN') 
    % La Niña (LN).
    elseif mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) < -0.5 || mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) < -0.5
        plot([120, 120], [range_years(1-4+defn_window*k), range_years(1+1+defn_window*k)], 'b', Linewidth=10, DisplayName='LN')
    end
end
xlim([120, 280])
clim([-3, 3])
xlabel('Longitude');
if j == 1
    ylabel('Model Simulated Year');
end
xticks(120:60:280)
xticklabels({'120°', '180°', '240°'})
yticks([range_years(1):1:range_years(end), range_years(end)])
yticklabels(string(round([range_years(1):1:range_years(end), range_years(end)])))
unqHandles = legendUnq(gca);
legend(unqHandles);
title('Hovmöller Diagram');
cb = colorbar('westoutside', Ticks=-3:3, TickLabels=arrayfun(@(x) sprintf('%d°C', x), -3:3, 'UniformOutput', false));
cb.Position(1) = cb.Position(1) - 0.05;
fontsize(16, 'points')

subplot(1, 7, 2)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(1.5 * u(time_start_plot:time_end_plot), time_start_plot:time_end_plot, 'b', Linewidth=2)
ylim([time_start_plot, time_end_plot])
title('Time Series of u')
xlabel('u (m/s)')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
grid on
fontsize(16, 'points')

subplot(1, 7, 3)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(150 * h_W(time_start_plot:time_end_plot), time_start_plot:time_end_plot, 'm', Linewidth=2)
ylim([time_start_plot, time_end_plot])
title('Time Series of h_W')
xlabel('h_W (m)')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
grid on
fontsize(16, 'points')

subplot(1, 7, 4)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(7.5 * T_C(time_start_plot:time_end_plot), time_start_plot:time_end_plot, 'color', [1, 0.38, 0], Linewidth=2)
ylim([time_start_plot, time_end_plot])
title('Time Series of T_C')
xlabel('T_C (°C)')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
grid on
fontsize(16, 'points')

subplot(1, 7, 5)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(7.5 * T_E(time_start_plot:time_end_plot), time_start_plot:time_end_plot, 'r', Linewidth=2)
ylim([time_start_plot, time_end_plot])
title('Time Series of T_E')
xlabel('T_E (°C)')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
grid on
fontsize(16, 'points')

subplot(1, 7, 6)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(5 * tau(time_start_plot:time_end_plot), time_start_plot:time_end_plot, 'color', [0.2, 0.49, 0.2], Linewidth=2)
ylim([time_start_plot, time_end_plot])
title('Time Series of τ')
xlabel('τ (m/s)')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
grid on
fontsize(16, 'points')

subplot(1, 7, 7)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(I(time_start_plot:time_end_plot), time_start_plot:time_end_plot, 'color', [0.72, 0.27, 1], Linewidth=2)
ylim([time_start_plot, time_end_plot])
title('Time Series of I')
xlabel('I')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
grid on
fontsize(16, 'points')

% Time series of the ACI metric and objective CIR length for y(t) → x, as well 
% as heatmap of the subjective CIR length over time and ε, all over the period
% of interest.
figure('WindowState', 'maximized');

subplot(1, 7, 1)
range_months = 1+12*(sim_year_start+2):12*(sim_year_start+ACI_period_years-2);
range_years = range_months/12;
range_idxs = range_months(1)*k_dt:range_months(end)*k_dt;
range_sim_years = range_idxs/k_dt/12;
defn_window = 12; 
total_loop = (ACI_period_years-4-1)*12/defn_window;
Hov = zeros(length(lon_sst), length(range_months));
xx = [ones(size(T_C_monthly(range_months)')), T_E_monthly(range_months)', T_C_monthly(range_months)'];
for i = 1:length(lon_sst)
    Hov(i, :) = xx * reg_te_tc_hw_ssta(i, 1:3)';
end
[xx, yy] = meshgrid(range_years, lon_sst);
contourf(yy, xx, Hov, 30, Linestyle='none', DisplayName='Niño3.4')
hold on
colormap('jet')
plot([180, 180], [range_years(1), range_years(end)], 'k--', Linewidth=2, DisplayName='Dateline');
temp_tau = range_idxs(1:round(k_dt/3):end);
plot(180 + 10*tau(temp_tau), range_sim_years(1:round(k_dt/3):end), 'color', [0.2, 0.49, 0.2], Linewidth=1.5, DisplayName='WWB (τ)');
plot(180 + 10*I(temp_tau), range_sim_years(1:round(k_dt/3):end), 'color', [0.72, 0.27, 1], Linewidth=1.5, DisplayName='Walker Cell (I)');
for k = 1:total_loop
    % Extreme EP El Niño (EEP EN).
    if mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > 1.0 && mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window))
        plot([120, 120], [range_years(1-4+k*defn_window), range_years(1+1+defn_window*k)], 'r', Linewidth=10, DisplayName='EEP EN') 
    % Moderate EP El Niño (MEP EN).
    elseif mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > 0.5 && mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window))
        plot([120, 120], [range_years(1-4+defn_window*k), range_years(1+1+defn_window*k)], 'm', Linewidth=10, DisplayName='MEP EN') 
    % CP El Niño (CP EN).
    elseif mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > 0.5 && mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) > mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window))
        plot([120, 120], [range_years(1-4+defn_window*k), range_years(1+1+defn_window*k)], 'color', [1, 0.38, 0], Linewidth=10, DisplayName='CP EN') 
    % La Niña (LN).
    elseif mean(T_E_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) < -0.5 || mean(T_C_monthly(range_months(1)-1+k*defn_window:range_months(1)+1+k*defn_window)) < -0.5
        plot([120, 120], [range_years(1-4+defn_window*k), range_years(1+1+defn_window*k)], 'b', Linewidth=10, DisplayName='LN')
    end
end
xlim([120, 280])
clim([-3, 3])
xlabel('Longitude');
if j == 1
    ylabel('Model Simulated Year');
end
xticks(120:60:280)
xticklabels({'120°', '180°', '240°'})
yticks([range_years(1):1:range_years(end), range_years(end)])
yticklabels(string(round([range_years(1):1:range_years(end), range_years(end)])))
unqHandles = legendUnq(gca);
legend(unqHandles);
title('Hovmöller Diagram');
cb = colorbar('westoutside', Ticks=-3:3, TickLabels=arrayfun(@(x) sprintf('%d°C', x), -3:3, 'UniformOutput', false));
cb.Position(1) = cb.Position(1) - 0.05;
fontsize(16, 'points')

subplot(1, 7, 2)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(ACI_metric(2*12*k_dt:12*12*k_dt), time_start_plot:time_end_plot, 'k', Linewidth=2)
ylim([time_start_plot, time_end_plot])
title('ACI Metric for y(t) \rightarrow x')
xlabel('Relative Entropy')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
grid on
fontsize(16, 'points')

subplot(1, 7, 3)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(5 * tau(time_start_plot:time_end_plot), time_start_plot:time_end_plot, 'color', [0.2, 0.49, 0.2], Linewidth=2)
hold on
temp = approx_objective_CIR(1:end-lookahead_tolerance*12*k_dt);
for j = time_start_plot:k_dt:time_end_plot
    plot([5*tau(j), 5*tau(j)], [j, j+approx_objective_CIR(j-time_start_plot+1)/dt], 'color', [0, 0, 0], Linewidth=2)
    plot(5*tau(j), j+approx_objective_CIR(j-time_start_plot+1)/dt, 'color', 'k', 'marker', '_', Linewidth=1)
end
ylim([time_start_plot, time_end_plot])
title('Time Series and Objective CIR Length for y(t) \rightarrow x')
xlabel('τ (m/s)')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
legend('τ', 'Objective CIR Length')
set(gca,'XGrid', 'off', 'YGrid', 'on')
fontsize(16, 'points')

subplot(1, 7, 4)
time_start_plot = (sim_year_start+2)*12*k_dt;
time_end_plot = (sim_year_start+ACI_period_years-2)*12*k_dt;
plot(approx_objective_CIR(1:end-lookahead_tolerance*12*k_dt), time_start_plot:time_end_plot, 'r', Linewidth=2)
hold on
plot(defn_objective_CIR, time_start_plot:time_end_plot, 'b', Linewidth=2)
ylim([time_start_plot, time_end_plot])
title({'Comparison of Objective', 'CIRs for y(t) \rightarrow x'})
xlabel('Objective CIR Length')
ylabel('Model Simulated Year')
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
legend('Underestimating Approximation', 'Definition')
grid on
fontsize(16, 'points')

subplot(1, 7, 5:7)
xx = eps_ord_values;
yy = time_start_plot:time_end_plot;
[X, Y] = meshgrid(xx, yy);
% pcolor(X, Y, (log10(subjective_CIR(:, 1:end-lookahead_tolerance*12*k_dt)./max_RE_metric(1:end-lookahead_tolerance*12*k_dt))).');
pcolor(X, Y, log10(subjective_CIR(:, 1:end-lookahead_tolerance*12*k_dt)).');
shading interp
colormap("jet")
clim([-2, 0.5])
cb = colorbar('eastoutside', Ticks=-2:0.5:0.5, TickLabels=arrayfun(@(x) sprintf('10^{%.1f}', x), -2:0.5:0.5, 'UniformOutput', false));
cb.Position(1) = cb.Position(1) + 0.1;
cb.Position(3) = 0.015;
xlim([lowest_order/2, highest_order])
set(gca, 'XDir','reverse')
title('Subjective CIR Length for y(t) \rightarrow x (Logarithmic Scale)')
xlabel('ε')
ylabel('Model Simulated Year')
xticks(round(lowest_order/2):highest_order)
xticklabels(arrayfun(@(x) sprintf('10^{%d}', x), round(lowest_order/2):highest_order, 'UniformOutput', false))
yticks(time_start_plot:12*k_dt:time_end_plot)
yticklabels(string((time_start_plot:12*k_dt:time_end_plot)/1200))
fontsize(16, 'points')