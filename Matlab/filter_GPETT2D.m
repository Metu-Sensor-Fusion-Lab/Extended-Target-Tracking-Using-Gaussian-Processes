function [ estState, estStateCov ] = filter_GPETT2D( processModel, prevEstState...
    , prevEstStateCov, measArray, paramGP, basisAngleArray, eps)
% FILTER_ETTGP2D basically implements the algorithm proposed in the paper

% Input:
%           processModel:           cell array of 2: {F, Q}
%           prevEstState:           previous state estimate (column vector)
%           prevEstStateCov:        previous covariance matrix
%           measArray:              each row is in the form [xMeas yMeas]
%           paramGP:                parameters of the underlying GP model
%                                   paramGP = {meanGP, stdPriorGP, stdRadiusGP, lengthScaleGP, stdMeasGP}
%           basisAngleArray:        angle array keeping the angles at which
%                                   the extent is maintained (a column vector)
%           eps:                    a small scalar, i.e., 1e-6
% Output:
%               estState:           state estimate (column vector)
%               estStateCov:        covariance matrix

% State space model of the system is as follows:
%           X_k = [c_k; psi_k; v_k; psiDot_k; f_k]      	state vector
%           c_k = [x_k; y_k]                                center of the object
%           v_k = [xDot_k; yDot_k]                        	velocity of the object
%           f_k = [f0_k; f1_k; ...; fN_k]               	extent of the object
%
% Process model:
%           X_(k+1) = F * X_k + w,  w~ N(0,Q)
% Measurement model:
%           y_k = h(X_k) + ek,      where e_k ~ N(0, R_k)


% Author:   Murat Kumru <kumru@metu.edu.tr>


% Extract matrices defining the proces model
F = processModel{1};
Q = processModel{2};
% Extract the necessary parameters of the GP model
lengthScaleGP = paramGP{4};
stdMeasGP = paramGP{5};

% Compute the inverse of P0_extent since it will repeatedly be used in the system dynamic model
persistent inv_P0_extent;
if isempty(inv_P0_extent)
    P0_extent = compute_GP_covariance(basisAngleArray, basisAngleArray, paramGP);   % Covariance of the target extent
    P0_extent = P0_extent + eps*eye(size(basisAngleArray,1));                       % To prevent numerical errors thrown during matrix inversion
    chol_P0_extent = chol(P0_extent);
    inv_chol_P0_extent = inv(chol_P0_extent);
    inv_P0_extent = inv_chol_P0_extent * inv_chol_P0_extent';
end


%% Process Update
predState = F * prevEstState;                	% Predicted state
predStateCov = F * prevEstStateCov * F' + Q;	% Covariance of the predicted state


%% Measurement Update
curNumMeas = size(measArray, 1);
numStates = size(F,1);

% Extract relevant information from the predicted state
predPos = predState(1:2);
predPsi = predState(3);
predExtent = predState(7:end);

% In the below loop, the following operations are performed for each 2D point measurement:
% 1. Compute measurement predictions from the predicted state by relying on
% the nonlinear measurement model.
% 2. Obtain the corresponding measurement covariance matrix.
% 3. Linearize the measurement model to employ in EKF equations.
predMeas = zeros(curNumMeas * 2, 1);                    % of the form [x1; y1; z1;...; xn; yn; zn]
measCov = zeros(curNumMeas*2);                          % measurement noise covariance matrix
linMeasMatrix = zeros(curNumMeas * 2, numStates);       % linearized measurement model
for i = 1:curNumMeas
    iMeas = transpose(measArray(i, :));    % Select one measurement of the form [xMeas; yMeas]
    
    % Compute the following variables to exploit in measurement model
    diffVector = iMeas - predPos;                           % The vector from the center to the measurement
    diffVectorMag = norm(diffVector);
    orientVector = diffVector / diffVectorMag;              % The orientation vector (of magnitude 1)
    angleGlobal = atan2(diffVector(2), diffVector(1));      % Angle in the global coordinate frame
    angleLocal = angleGlobal - predPsi;                  	% Angle in the local coordinate frame
    angleLocal = mod(angleLocal, 2*pi);
    
    covMeasBasis = compute_GP_covariance(angleLocal, basisAngleArray, paramGP);
    
    H_f = covMeasBasis * inv_P0_extent; % GP model relating the extent and the radial
    ... function value evaluated at the current measurement angle
        
    % Obtain the measurement prediction by the original nonlinear meas model
    iPredMeas = predPos + orientVector * H_f * predExtent; % [xPredicted; yPredicted]
    predMeas(1+(i-1)*2 : i*2) = iPredMeas;
    
    % Obtain the covariance of the current measurement set
    % Definition: iCovMeas = k(uk,uk), uk: argument of the current measurement
    iMeasCov = compute_GP_covariance(angleLocal, angleLocal, paramGP);
    % Definition: iR = k(uk,uk) - k(uk, uf) * inv(K(uf,uf)) * k(uf, uk)
    iR_f = iMeasCov - covMeasBasis * inv_P0_extent * covMeasBasis';
    measCov(((i-1)*2+1):i*2, ((i-1)*2+1):i*2) = stdMeasGP^2 * eye(2) + ...
        orientVector * iR_f * orientVector';
    
    % Obtain the linearized model (Details in the Appendix B of the mentioned paper)
    dP_dW = ((diffVector*diffVector')./ diffVectorMag^3 - eye(2)./diffVectorMag);
    dHf_du = -1/ lengthScaleGP^2 * sin(angleLocal - basisAngleArray') .* covMeasBasis...
        * inv_P0_extent;
    dThetaG_dw = 1/ diffVectorMag^2 * [diffVector(2) -diffVector(1)];
    
    H_lin_pos = eye(2) + dP_dW * (H_f * predExtent)...
        + orientVector * dThetaG_dw * (dHf_du * predExtent);
    H_lin_psi = -orientVector * dHf_du * predExtent;
    H_lin_extent = orientVector * H_f;
    
    % Linearized measurement model for the current measurement
    % = [dh/dxc, dh/dpsi, zeros(due to diff wrt velocities), dh/dextent]
    linMeasMatrix(1+(i-1)*2:i*2, :) = [H_lin_pos H_lin_psi zeros(2,3) H_lin_extent];
end

% Put the measurements in a column the form: [x1; y1;...; xn; yn]
tempArray = transpose(measArray);
curMeasArrayColumn = tempArray(:);

% Realize measurement update
kalmanGain = predStateCov * linMeasMatrix' ...
    / (linMeasMatrix*predStateCov*linMeasMatrix' + measCov);
estState = predState + kalmanGain * (curMeasArrayColumn - predMeas);
estStateCov = (eye(numStates) - kalmanGain*linMeasMatrix) * predStateCov;
estStateCov = (estStateCov + estStateCov')/2; % To make the covariance matrix symmetric (needed due to numeric errors)

end
