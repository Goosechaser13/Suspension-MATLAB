%Suspension analysis

%% Mounting Positions (ENTER COORDS RELATIVE TO SINGLE POINTS)
% Must be Coplanar suspension
%% IGNORE Z, GET HORI,VERT DISTANCE FROM A REFERENCE POINT
% in X ; Y ; Z format
MSpring = [-9.359109; 15.858556 ;0]; % Spring mount to chassis X Y Z in.
MBell = [-11;7.605988;0]; %Bellcrank pivot point
MPush = [-21.72;4.21;0]; % Where pusharm attaches to lower arm

MBellSpr = [-14.404428; 11.010416 ;0]; %Point where bellcrank connects to spring
MBellPush = [-15.63532; 9.434057 ;0]; % bellcrank connects to rod
MBellARB = [ -13.774781 ; 8.700301 ; 0];

MArmL = [-7.98395 ; 3 ; 6.35]; % Where lower arm connects to chassis
MArmU = [-9.633238 ; 7.305803 ; 6.35]; % upper arm to chassis

MArmL_End = [-20.567206; 3.142493; 0]; % Where knuckle connects to lower arm
MArmU_End = [-21.776179 ; 10.6 ; -0]; %Where knuckle connects to upper arm

minWheelTravel = 1.9685; %minimum travel in inches set by rules
v_knuckle = (MArmU_End-MArmL_End)/norm(MArmU_End-MArmL_End)
MWheel = v_knuckle*(0.97) + MArmL_End;

MARB = [-7.98395 ; 3 ; 4.6];

%% 
points = [MSpring,MBell,MPush,MBellPush,MBellSpr, ... 
    MArmL, MArmU, MArmL_End , MArmU_End,MWheel,MBellARB,MARB];
% Points dictionary:
%  1 : Spring Mount   2 : Bellcrank Mount
%  3 : Pushrod Mount  4 : Pushrod Bell M
%  5 : Spring Bell M  6 : Lower Arm M
%  7 : Upper Arm M    8 : Lower Arm end
%  9 : Upper Arm end  10 :Wheel Mount
% 11 : ARB Bell M     12: ARB Mount
% Damper settings
x_rest = 7; %rest length for damper
x_min = 5.94;
x_max = 8.19;
%% Analysis

%Create motion ratio graph
%  This is done by change damper size and evaluting ratio of the wheel change
% MR = deltaSpring / deltaWheel


n = 15 ; % num of iterations
totalPoints = zeros([3,length(points),n]);
x = linspace(x_min,x_max,n); %damper length
for i = 1:length(x)
    totalPoints(:,:,i) = DamperX_to_Points(points,x(i));
end
wheelHeight = squeeze(totalPoints(2,10,:));
% set wheelHeight relative to the default point:
wheelHeight = wheelHeight - MWheel(2);
deltaX = diff(x);
deltaWheel = diff(wheelHeight);
MR = gradient(x,wheelHeight);
 
figure(1);
subplot(1,2,1);
plot(wheelHeight,x)
title("Motion Ratio");

if( (max(wheelHeight) - min(wheelHeight)) < minWheelTravel )
    text(0.5,0.5,0,'BELOW MIN WHEEL TRAVEL','Units','normalized',...
        'HorizontalAlignment','center','VerticalAlignment','middle' ...
        ,Color=[1,0,0],FontWeight='bold') %% If wheel doesn't travel enough, warn
end
xlabel('Wheel Distance from rest'); ylabel('Damper Extension');

subplot(1,2,2)
plot(wheelHeight,MR);
title("Instantaneous Motion Ratio");
ytickformat('%.2f');
xlabel('Wheel Height distance from rest'); ylabel('Motion Ratio');

visualizeSuspension(totalPoints(:,:,1));

%% Move Geometry by Damper Length


%% Visualize function
function visualizeSuspension(coords)
    figure(2)
    X = coords(1,:);
    Y = coords(2,:);
    Z = coords(3,:);
    
    scatter3(X,Z,Y,36,'filled')
    grid on;
    xlabel('X'); ylabel('Z'); zlabel('Y');
    %xlim([-100,50]);
    %zlim([-100,50]);
    %zlim([-50,100]);
    axis equal;
    view(0,0);
    title('3D visualization');
    rotate3d on;
    hold on 
    fill3([X(2),X(4),X(5)], [Z(2)+0.01,Z(4)+0.01,Z(5)+0.01],[Y(2),Y(4),Y(5)],'g')
    drawLine(coords(:,1),coords(:,5),'-ob') %spring line
    drawLine(coords(:,3),coords(:,4),'-ob')
    drawLine( coords(:,6),coords(:,8),'-or')
    drawLine( coords(:,7),coords(:,9),'-or')
    
    hold off
end


function drawLine(p1,p2,parameters)
    plot3( [p1(1),p2(1)] ,[p1(3),p2(3)],[p1(2),p2(2)],parameters,'LineWidth',2,'MarkerFaceColor','b' );
end
