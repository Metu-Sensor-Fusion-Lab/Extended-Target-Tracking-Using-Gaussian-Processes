function [] = visualize_tracking_outputs(measurements, estState, estStateCov, groundTruth, extentAnglesLocal, time)

% Extract relevant information from the state
estPos = estState(1:2);
estPsi = estState(3);
estExtent = estState(7:end);

% Extract the variance of the extent
stdState = sqrt(diag(estStateCov));
stdExtent = stdState(7:end);

% Extract ground truth values of the position and the orientation angle
objType = groundTruth.objectDescription(1);
objParam = groundTruth.objectDescription(2:end);
gtKinematics = groundTruth.dataLog(abs(groundTruth.dataLog(:,1)-time)<1e-10, 2:end);
gtCenter =gtKinematics(1:2)';
gtPsi = gtKinematics(3);

extentAnglesGlobal = extentAnglesLocal + estPsi; % Calculate angles of the extent in global frame

cla(); % Clear the current axis

%% Plot the ground truth
switch objType
    case 1 % Circle
        plot_circle(gtCenter, objParam, [0.4 0.7 0.3]);
    case 2 % Square
        plot_square(gtCenter, gtPsi, objParam, [0.4 0.7 0.3]);
    case 3 % Triangle
        plot_triangle(gtCenter, gtPsi, objParam, [0.4 0.7 0.3]);
end


%% Plot measurements
plot(measurements(:,1), measurements(:,2), 'rx', 'LineWidth', 2);

%% Plot the center position
plot(estPos(1), estPos(2), 'k+', 'LineWidth', 2);

%% Plot the heading by a line
plot([estPos(1) estPos(1)+3*cos(estPsi)], [estPos(2) estPos(2)+3*sin(estPsi)], 'k', 'LineWidth', 2);

%% Plot the estimated circumference
[xEstimated, yEstimated] = pol2cart(extentAnglesGlobal, estExtent);
xEstimated = xEstimated + estPos(1); % Shift to the estimated center position 
yEstimated = yEstimated + estPos(2); % Shift to the estimated center position 

xEstimatedPlot = [xEstimated; xEstimated(1)]; % To make it appear as a enclosed polygon
yEstimatedPlot = [yEstimated; yEstimated(1)]; % To make it appear as a enclosed polygon
plot(xEstimatedPlot, yEstimatedPlot, 'Color', [0 0.4 0.9], 'LineWidth', 2);

%% Plot the confidence region of 1-std of the extent
% The uncertainity will be represented with a filled polygon
% Compute the vertices of the inner polygon
[xInner, yInner] = pol2cart(extentAnglesGlobal, max(estExtent - stdExtent, 0)); % Extent can not be smaller than 0
xInner = xInner + estPos(1); % Shift to the estimated center position 
yInner = yInner + estPos(2); % Shift to the estimated center position 

% Compute the vertices of the outer polygon
[xOuter, yOuter] = pol2cart(extentAnglesGlobal, (estExtent + stdExtent));
xOuter = xOuter + estPos(1); % Shift to the estimated center position 
yOuter = yOuter + estPos(2); % Shift to the estimated center position 

% Prepare the arrays for plotting 
xInnerPlot = [xInner; xInner(1)]; % To make it appear as a enclosed polygon
yInnerPlot = [yInner; yInner(1)]; 
xOuterPlot = [xOuter; xOuter(1)]; 
yOuterPlot = [yOuter; yOuter(1)]; 

hFill = fill([xOuterPlot; xInnerPlot] , [yOuterPlot; yInnerPlot], 'b');
hFill.FaceColor = [0.6 0.7 1];
hFill.EdgeColor = 'none';
hFill.FaceAlpha = 0.5;
uistack(hFill,'bottom');

drawnow();

end


function plot_circle(center, objectParameters, color)
% Plots a circle

numControlPoints = 100;
controlAngleArray = transpose(linspace(0, 2*pi, numControlPoints));

% Extract object parameters
radius = objectParameters;% Ellipsoid

% Produce cartesian points on the unit sphere
[xCircle, yCircle] = pol2cart(controlAngleArray, radius*ones(numControlPoints));

% Translate these points
xCircle = xCircle+ center(1);
yCircle = yCircle + center(2);

plot(xCircle, yCircle, 'Color', color, 'LineWidth', 2);
end


function plot_square(center, psi, objectParameters, color)
% Plots a square

% Extract object parameters
edgeLen = objectParameters(1);

% Define the vertices of the box
vertices_L = [-edgeLen  -edgeLen;
    edgeLen  -edgeLen;
    edgeLen  edgeLen;
    -edgeLen  edgeLen;
    -edgeLen  -edgeLen] * 0.5;

% Compute the rotation matrix regarding the orientation angle
rotMatrix = compute_rotation_matrix(psi); 

% Transform the local points into the global frame
vertices_G = transpose(rotMatrix * vertices_L');    % Rotation
vertices_G = [vertices_G(:,1)+center(1), vertices_G(:,2)+center(2)]; % Translation

plot(vertices_G(:,1), vertices_G(:,2), 'Color', color, 'LineWidth', 2);
end


function plot_triangle(center, psi, objectParameters, color)
% Plots a triangle

% Extract object parameters
sideEdgeLength = objectParameters(1);
bottomEdgeLength =objectParameters(2);
height = sqrt(sideEdgeLength^2-(bottomEdgeLength/2)^2);

% Define the vertices of the box
vertices_L = [-height/3  -bottomEdgeLength/2;
    2/3*height  0;    
    -height/3 bottomEdgeLength/2;
    -height/3  -bottomEdgeLength/2];

% Compute the rotation matrix regarding the orientation angle
rotMatrix = compute_rotation_matrix(psi); 

% Transform the local points into the global frame
vertices_G = transpose(rotMatrix * vertices_L');    % Rotation
vertices_G = [vertices_G(:,1)+center(1), vertices_G(:,2)+center(2)]; % Translation

plot(vertices_G(:,1), vertices_G(:,2), 'Color', color, 'LineWidth', 2);
end


function rotationMatrix =  compute_rotation_matrix(psi)
rotationMatrix = [cos(psi) -sin(psi); sin(psi) cos(psi)];
end
