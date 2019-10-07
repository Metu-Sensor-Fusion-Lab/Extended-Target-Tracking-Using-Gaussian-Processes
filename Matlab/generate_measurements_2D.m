function [measMatrix, groundTruth] = generate_measurements_2D(measParamaters)
% This function produces point measurements originated from the contour of
% a specified object which is moving according to the specied patter type.

% Input:
%           measParamaters:     the parameters utilized

% Output:
%           measMatrix:         the measurement matrix comprised of rows
%                               in the form: [timeStamp xPosition yPosition]
%           groundTruth:       the ground truth information of the corresponding scenario.
%                               It is a structure with the following fields: 'objectDescription', 'dataRecord'.
%                               'objectDescription' includes object type (e.g., sphere, box...) and
%                               the parameters of the geometry.
%                               Each row of the 'dataRecord' is in the form of [timeStamp center' psi]
%                               where centerPosition = [cX cY]' and psi describes the orientation angle.


% Author:   Murat Kumru <kumru@metu.edu.tr>


% Extract measurement parameters
experimentDuration = measParamaters{1};
timeStep = measParamaters{2};
numMeasPerInstant = measParamaters{3};
stdMeas = measParamaters{4};
objType = measParamaters{5};
objParameters = measParamaters{6};
motionType = measParamaters{7};
isPartiallyObserved = measParamaters{8};

numIterations = ceil(experimentDuration/ timeStep);
measMatrix = zeros(numIterations*numMeasPerInstant, 3);
gtDataLog = zeros(numIterations, 4); % : [timeStamp centerPosition' psi]

positionObj = [0 0];    % Initialize the center position of the object
orientationObj = -pi + 2* pi * rand; % Randomly select a constant orientation
for i = 1:numIterations
    
    curTime = (i-1)*timeStep;
    % Simulate the motion of the object
    switch motionType
        case 1 % No movement
            velocityObj = [0 0];
            
        case 2 % Constant velocity
            velocityObj = [0.3 0.2];
            
        case 3 % Sinusoidal motion
            velocityObj = [0.05 sin(0.01*2*pi*curTime)];
    end
    positionObj = positionObj + timeStep * velocityObj;       
        
    % Produce contour measurements regarding the center
    curMeasArray = produceContourMeas(positionObj, orientationObj, numMeasPerInstant, stdMeas, objType,...
        objParameters, isPartiallyObserved);
    
    % Log the groundtruth
    gtDataLog(i, :) = [curTime positionObj orientationObj];
    
    % Insert the measurements to the output matrix with relevant time stamp
    measMatrix((i-1)*numMeasPerInstant + 1 : (i*numMeasPerInstant), :) = ...
        [ones(numMeasPerInstant, 1)*curTime curMeasArray];
end

groundTruth = struct('objectDescription', [objType, objParameters]...
            , 'dataLog', gtDataLog);
end


function [posArray] = produceContourMeas(positionObj, orientationObj, numMeasurements, stdMeas, objType...
    , objParameters, isPartiallyObserved)
switch objType
    case 1 % Circle
        radius = objParameters;
        
        if isPartiallyObserved
            thetaArray = pi+ pi * rand(numMeasurements, 1);
            rhoArray = radius + stdMeas * randn(numMeasurements, 1);
        else
            thetaArray = 2*pi * rand(numMeasurements, 1);
            rhoArray = radius + stdMeas * randn(numMeasurements, 1);
        end
        
        [xArray, yArray] = pol2cart(thetaArray, rhoArray);
        
    case 2 % Square
        edgeLen = objParameters;
        
        % First, produce a random array for the selection of the edges
        if isPartiallyObserved
            edgeArray = randi([1 3], numMeasurements, 1); % Upper vertex is not sampled
        else
            edgeArray = randi([0 3], numMeasurements, 1);
        end
        
        % Initialize the outputs
        xArray = zeros(numMeasurements, 1);
        yArray = zeros(numMeasurements, 1);
        pointer = 1;
        for iEdge = 0:3
            num = sum(edgeArray == iEdge);
            switch iEdge
                case 0 % variable X; Y is set to positive edge
                    xArray(pointer:(pointer+num-1)) = -edgeLen/2 + edgeLen*rand(num, 1);
                    yArray(pointer:(pointer+num-1)) = edgeLen/2*ones(num,1) + stdMeas * randn(num, 1);
                    
                case 1 % variable X; Y is set to negative edge
                    xArray(pointer:(pointer+num-1)) = -edgeLen/2 + edgeLen*rand(num, 1);
                    yArray(pointer:(pointer+num-1)) = -edgeLen/2*ones(num,1) + stdMeas * randn(num, 1);
                    
                case 2 % variable Y; X is set to positive edge
                    xArray(pointer:(pointer+num-1)) = edgeLen/2*ones(num,1) + stdMeas * randn(num, 1);
                    if isPartiallyObserved
                        yArray(pointer:(pointer+num-1)) = -edgeLen/2 + edgeLen/2*rand(num, 1);
                    else
                        yArray(pointer:(pointer+num-1)) = -edgeLen/2 + edgeLen*rand(num, 1);
                    end
                    
                case 3 % variable Y; X is set to negative edge
                    xArray(pointer:(pointer+num-1)) = -edgeLen/2*ones(num,1) + stdMeas * randn(num, 1);
                    if isPartiallyObserved
                        yArray(pointer:(pointer+num-1)) = -edgeLen/2 + edgeLen/2*rand(num, 1);
                    else
                        yArray(pointer:(pointer+num-1)) = -edgeLen/2 + edgeLen*rand(num, 1);
                    end
            end
            pointer = pointer + num;
        end
        
    case 3 % Triangle
        % It is assumed to be a isosceles triangle
        sideEdgeLength = objParameters(1);
        bottomEdgeLength = objParameters(2);
        height = sqrt(sideEdgeLength^2-(bottomEdgeLength/2)^2);
        
        % Initialize the outputs
        xArray = zeros(numMeasurements, 1);
        yArray = zeros(numMeasurements, 1);
        for j = 1:numMeasurements
            if rand < bottomEdgeLength/ (bottomEdgeLength+2*sideEdgeLength)
                % Measurement from bottom edge
                if isPartiallyObserved
                    xArray(j) = -1/3*height + stdMeas*randn;
                    yArray(j) = (bottomEdgeLength*rand - bottomEdgeLength/2);
                else
                    xArray(j) = -1/3*height + stdMeas*randn;
                    yArray(j) = (bottomEdgeLength*rand - bottomEdgeLength/2);
                end
            else
                % Measurement from one of the side edges
                if isPartiallyObserved
                    xArray(j) = (height*rand -1/3*height);
                    yArray(j) = -bottomEdgeLength/2*(1-(xArray(j)/height + 1/3)) + stdMeas*randn;
                else
                    xArray(j) = (height*rand -1/3*height);
                    randomSign = (rand > 0.5)*2 - 1;
                    yArray(j) = randomSign * bottomEdgeLength/2*(1-(xArray(j)/height + 1/3)) + stdMeas*randn;
                end
            end
        end
end

% Compute the rotation matrix
rotMatrix = [cos(orientationObj) -sin(orientationObj); sin(orientationObj) cos(orientationObj)];
% Rotate measurements
rotXYArray = rotMatrix * [xArray'; yArray'];

xArray = rotXYArray(1, :)';
yArray = rotXYArray(2, :)';


posArray = [xArray+positionObj(1) yArray+positionObj(2)];
end

