
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Vestibular Afferent Simulation Pipeline
%
%  PURPOSE
%  --------
%  Convert raw six-degree-of-freedom head-mounted-IMU data into vestibular
%  afferent firing-rate estimates and save the results as a CSV file.
%
%  WHAT THE SCRIPT DOES
%  1.  Loads a CSV file containing linear acceleration (m/s²) and angular
%      velocity (deg s⁻¹) measured at an off-centre IMU.
%  2.  Transforms the motion signals to the head-centre reference frame:
%          • compensates for the IMU lever arm *r*  
%          • converts right-ear semicircular-canal axes  
%          • high-passes everything above 0 Hz and low-passes at 20 Hz
%  3.  Applies fractional-order transfer functions (from
%      **VestibularTFs**) to obtain:
%          – Otolith Regular and Irregular afferents (3 axes each)  
%          – Canal Regular and Irregular afferents  (3 axes each)
%  4.  Trims the first 50 samples (filter transients) and writes the
%      resulting 12-column matrix to  
%          ‘C:\Users\Data_Science_Team\Vestibular_Afferents.csv’
%
%  INPUTS
%  ------
%      C:\Users\Data_Science_Team\Head_IMU.csv   — [N × 6] numeric CSV  
%          [accel_x  accel_y  accel_z  vel_x  vel_y  vel_z]
%
%  KEY PARAMETERS
%      Fs      … sampling rate (500 Hz)  
%      r       … IMU position relative to head centre [m]  
%      R_matrix… right-ear canal orientation unit-vectors  
%      Butter  … 4-pole, 20 Hz low-pass (applied with *filtfilt*)
%
%  OUTPUT
%  ------
%      Vestibular_Afferents.csv  — [N-50 × 12]  
%          Columns: OtoReg(1:3) | OtoIrreg(1:3) | CanReg(1:3) | CanIrreg(1:3)
%
%  DEPENDENCIES
%      • VestibularTFs_DS.m     (transfer-function container)  
%      • Signal Processing Toolbox (for *butter*, *filtfilt*)  
%      • FOMCON toolbox (fractional-order TF support used by
%        VestibularTFs_DS)
%
%  USAGE
%      Simply run the script — edit the *data* and *output* paths or the
%      sampling rate if your dataset differs.
%
%  AUTHOR / DATE
%      Data-Science Team · 2025-04-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load IMU data (45000 x 6) array
data = readmatrix('C:\Users\A02236463\Box\001 PROJECTS\A022 - Vestibular Engagement\Data_Science_Team\Head_IMU.csv');  

% Sampling Rate
Fs = 500;

% Position of the IMU with respect to the head center
r = [-0.1 0 -0.04];

% 1) Instantiate the transfer-function container
vest = VestibularTFs_DS;

% 2) Retrieve motion-to-afferent transfer functions
[H_Canal_Motion_Regular, H_Canal_Motion_Irregular, ~] = vest.canalMotionTFs();
[Hoto_Regular, Hoto_Irregular, ~] = vest.otolithMotionTFs();

% 3) Separate acceleration and angular-velocity channels
accel = data(:, 1:3)./9.81;     % linear acceleration [g]
angvel = data(:, 4:6);          % angular velocity [deg s⁻¹]

% 4) Calculate angular acceleration
angacc = [0 0 0;diff(deg2rad(angvel))*Fs];

% 5) Transform IMU data to the three canals 
% LA = [-0.58954; -0.78750; 0.17971];
% LH = [-0.31645; 0.04108; -0.94772];
% LP = [-0.69611; 0.66820; 0.26257]; 

RA = [0.58905; -0.78928; -0.17334];
RH = [0.32891; -0.03566; 0.94369]; 
RP = [0.69248; 0.66561; -0.27826]; 

% 6) get rotation matrix for right and left side canal components
% L_matrix = [LA LH LP];
R_matrix = [RA RH RP];

% 7) project the angular velocities onto the canal components done only for right ear
% L_comp = L_matrix'*sigIn;
angvel = angvel*R_matrix;

% 8) Create a matrix of r's to avoid a loop
r = repmat(r,[length(angvel) 1]);

% 9) Transform linear acceleration to head center
accel = accel + cross(angacc, r) + cross(deg2rad(angvel),cross(deg2rad(angvel),r));

% 10) Define the low pass filter as the vestibular transfer functions are
% only valid up to 20Hz
[B,A] = butter(4,20/250);

% 11) Filter data
accel = filtfilt(B,A, accel);
angvel = filtfilt(B,A, angvel);

% 12) Build a time vector for simulation
N  = size(data, 1);
t  = (0:N-1)' / Fs;

% 13) Pre-allocate outputs
OtoReg   = zeros(N,3);
OtoIrreg = zeros(N,3);
CanReg   = zeros(N,3);
CanIrreg = zeros(N,3);

% 14) Apply each TF to its corresponding column
for k = 1:3
    OtoReg(:,k)   = lsim(Hoto_Regular,   accel(:,k), t);  % otolith regular
    OtoIrreg(:,k) = lsim(Hoto_Irregular, accel(:,k), t);  % otolith irregular
    CanReg(:,k)   = lsim(H_Canal_Motion_Regular, angvel(:,k), t); % canal regular
    CanIrreg(:,k) = lsim(H_Canal_Motion_Irregular, angvel(:,k), t); % canal irregular
end

% 15) Save OtoReg, OtoIrreg, CanReg, CanIrreg, spk, fr as CSV
Vestibular_Afferents = [OtoReg, OtoIrreg, CanReg, CanIrreg];

% Remove the first 50 points to remove artifacts
Vestibular_Afferents = Vestibular_Afferents(50:end, :);

writematrix(Vestibular_Afferents, 'C:\Users\Data_Science_Team\Vestibular_Afferents.csv'); 
