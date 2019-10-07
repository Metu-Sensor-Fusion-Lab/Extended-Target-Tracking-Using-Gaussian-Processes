% This file demonstrates an example utilization of the algorithm proposed
...in the following paper:
    
% "Extended target tracking using Gaussian processes" (IEEE Trans. on Signal Process.,
...vol. 63, no. 16, pp. 4165–4178, Apr.  2015) by N. Wahlström and E. Özkan

% Dependencies: None


% Author:   Murat Kumru <kumru@metu.edu.tr>
% Date:     07.10.2019


clc;
clear all;
close all;

%% Paramaters
% Simulation parameters
expDuration = 30;           % Experiment duration (in seconds)
T = 1;                      % Sampling time (in seconds)
numMeasPerInstant = 20;     % Number of point measurements acquired in an instant
stdMeas = 0.1;              % Standart deviation of the measurements
isMeasurementPartial = 0;	% Flag for partial measurement
motionType = 2;             % 1: Still, 2: Constant velocity, 3: Sinusoidal motion
eps = 1e-6;                 % A small scalar
% Object parameters
objType = 2;                % objType = 1:Circle, 2:Square, 3:Triangle
switch objType
    case 1
        radius = 3;
        objParameters = radius;
    case 2
        edgeLength = 3;
        objParameters = edgeLength;
    case 3
        sideEdgeLength = 5;
        bottomEdgeLength = 4;
        objParameters = [sideEdgeLength bottomEdgeLength];
    case 4
        objParameters = [];
end
paramMeas = {expDuration, T, numMeasPerInstant, stdMeas, objType, objParameters...
    , motionType, isMeasurementPartial};

% GP parameters
numBasisAngles = 50;        % The number of basis angles (on which the extent is to be maintained)
meanGP = 0;                 % Mean of the GP model
stdPriorGP = 2;             % Prior standart deviation
stdRadiusGP = 0.8;          % Radius standart deviation
lengthScaleGP = pi/8;       % Length scale (smaller lengthscale implies edgier/ more spiky extents)
stdMeasGP = stdMeas;        % Measurement std in the GP model (it may be set different from the original value)
paramGP = {meanGP, stdPriorGP, stdRadiusGP, lengthScaleGP, stdMeasGP};

% Process paramaters
stdCenter = 1e-2;       % Std dev of the position process noise
stdPsi = 1e-4;          % Std dev of the orientation angle (psi) process noise
alpha = 1e-4;           % Forgetting factor (to used in the process model of the extent)
paramEKF = [stdCenter stdPsi alpha];

% Determine the basis angles
basisAngleArray = transpose(linspace(0, 2*pi, numBasisAngles+1));
basisAngleArray(end) = [];  % Trash the end point to eleminate the repetition at 2*pi


%% Generate Measurements and Ground Truth
[measComplete, groundTruth] = generate_measurements_2D(paramMeas); 


%% Determine the initial distribution
% Initial means
c0 = transpose(mean(measComplete(1:numMeasPerInstant, 2:3), 1));     % Compute Initial center as the average of the first set of measurements
v0 =  zeros(2,1);   % Initial velocity
psi0 = 0;           % Initial orientation (yaw) angle
psiDot0 = 0;        % Initial yaw rate
f0 = meanGP * ones(numBasisAngles, 1);  % Initial extent
% Initial covariances
P0_center = 10 * eye(2);
P0_velocity = 1 * eye(2);
P0_psi = 1e-5;
P0_psiDot = 1e-5;
P0_extent = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP);
P0 = blkdiag(P0_center, P0_psi, P0_velocity, P0_psiDot, P0_extent);

% Initialize the filter estimate
estState = [c0; psi0; v0; psiDot0; f0];
estStateCov = P0;


%% Process Model
F_tracking = kron([1 T; 0 1], eye(3));          	% Constant velocity model
F_extent = exp(-alpha*T) * eye(numBasisAngles);     % A forgetting factor is included
F = blkdiag(F_tracking, F_extent);

Q_tracking = kron([T^3/3 T^2/2; T^2/2 T], diag([stdCenter^2 stdCenter^2 stdPsi^2]));
Q_extent = (1-exp(-2*alpha*T)) * P0_extent;
Q = blkdiag(Q_tracking, Q_extent);

processModel = {F, Q};


%% Open a Figure (to illustrate tracking outputs)
figure;
grid on; hold on;
ylim([min(measComplete(:,3))-2 max(measComplete(:,3))+2]);
xlim([min(measComplete(:,2))-2 max(measComplete(:,2))+2]);
ylabel('Y axis');
xlabel('X axis');
title('Extended Target Tracking using Gaussian Processes');

%% Initialize Simulation Log
numInstants = ceil(expDuration/ T);
simLog.Parameters = {paramMeas, paramGP, paramEKF, processModel, basisAngleArray};
simLog.GroundTruth = groundTruth;
%simLog.GroundTruth = groundTruth;
simLog.TrackingData(numInstants) = struct('time', [], 'measurement', []...
    , 'stateEstimated', [], 'stateCovariance', []);


%% GPETT2D LOOP
time = 0;
for k = 1:numInstants
    
    % Extract current measurements
    curMeasArray = measComplete(abs(measComplete(:, 1)-time)<eps, 2:3);
    
    % Call the GPETT2D filter
    [estState, estStateCov] = filter_GPETT2D(processModel, estState, estStateCov...
        , curMeasArray, paramGP, basisAngleArray, eps);
    
    % Visualize the results
    visualize_tracking_outputs(curMeasArray, estState, estStateCov, groundTruth...
        , basisAngleArray, time);
    
    % Log data
    simLog.TrackingData(k) = struct('time', time, 'measurement', curMeasArray...
        , 'stateEstimated', estState, 'stateCovariance', estStateCov);
    
    % Update time
    time = time + T;    
end
