%Suspension analysis

%% Mounting Positions (ENTER COORDS RELATIVE TO SINGLE POINTS)
% Must be Coplanar suspension
%% IGNORE Z, GET HORI,VERT DISTANCE FROM A REFERENCE POINT
% in X ; Y ; Z format
MSpring = [-9.359109; 15.858556 ;0]; % Spring mount to chassis X Y Z in.
MBell = [-9.3;7;0]; %Bellcrank pivot point
MPush = [-18.89;3;0]; % Where pusharm attaches to lower arm

MBellSpr = [-12.403644; 10.67471 ;0]; %Point where bellcrank connects to spring
MBellPush = [-13.763625; 9.208269 ;0]; % bellcrank connects to rod
MBellARB = [ -11.833428 ; 8.25485 ; 0];

MArmL = [-7.98395 ; 3 ; 6.93]; % Where lower arm connects to chassis
MArmU = [-9.633238 ; 7.305803 ; 6.93]; % upper arm to chassis

MArmL_End = [-20.567206; 3.142493; 0]; % Where knuckle connects to lower arm
MArmU_End = [-21.776179 ; 10.6 ; -0]; %Where knuckle connects to upper arm

minWheelTravel = 1.9685; %minimum travel in inches set by rules
v_knuckle = (MArmU_End-MArmL_End)/norm(MArmU_End-MArmL_End);
MWheel = v_knuckle*(3.12) + MArmL_End;

MARB = [-10.6 ; 4 ; 4];
L_ARB = 3;
MARB_End = [-10.6 ; 5.02 ; MARB(3)-L_ARB];

MBellARB = [-13.201524 ; 7.998113 ; 0.000000];
MARB     = [-12 ; 2.17 ; 2.9];
MARB_End = [-12 ; 2.23 ; .1];

v_arm =  (MARB_End - MARB) / norm(MARB_End - MARB);
L2 = 3;
v_arbbell =(MBellARB - MBell) / norm(MBellARB - MBell);
r = 4;
%MBellARB = (v_arbbell * r) + MBell
%MARB_End = (v_arm * L2) +MARB
%% 
points = [MSpring,MBell,MPush,MBellPush,MBellSpr, ... 
    MArmL, MArmU, MArmL_End , MArmU_End,MWheel,MBellARB,MARB,MARB_End];
% Points dictionary:
%  1 : Spring Mount   2 : Bellcrank Mount
%  3 : Pushrod Mount  4 : Pushrod Bell M
%  5 : Spring Bell M  6 : Lower Arm M
%  7 : Upper Arm M    8 : Lower Arm end
%  9 : Upper Arm end  10 :Wheel Mount
% 11 : ARB Bell M     12: ARB Mount
% 13 : ARB End
% Damper settings
x_rest = 7; %rest length for damper
x_min = 5.94;
x_max = 8.19;
%% Analysis
fprintf("DropLink: %.4f \n",norm(points(:,11)-points(:,13)))
fprintf("Bar Length: %.4f \n",MARB(1)*-2)
fprintf("Arm Length: %.4f \n",norm(points(:,12)-points(:,13)))
%Create motion ratio graph
%  This is done by change damper size and evaluating ratio of the wheel change
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

%MR ARB 
ARBEnd = squeeze(totalPoints(2,13,:));
MR_ARB = gradient(ARBEnd,wheelHeight);
 
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

%visualizeSuspension(totalPoints(:,:,7));
visualizeSuspension(points());

% PLOT ARB MR
figure(3);
subplot(1,2,1);
plot(wheelHeight,ARBEnd)
title("Motion Ratio of ARB");

if( (max(wheelHeight) - min(wheelHeight)) < minWheelTravel )
    text(0.5,0.5,0,'BELOW MIN WHEEL TRAVEL','Units','normalized',...
        'HorizontalAlignment','center','VerticalAlignment','middle' ...
        ,Color=[1,0,0],FontWeight='bold') %% If wheel doesn't travel enough, warn
end
xlabel('Wheel Distance from rest'); ylabel('ARB Travel');

subplot(1,2,2)
plot(wheelHeight,MR_ARB);
title("Instantaneous Motion Ratio");
ytickformat('%.2f');
xlabel('Wheel Height distance from rest'); ylabel('Motion Ratio ARB');

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
    fill3([X(2),X(4),X(11)], [Z(2)+0.01,Z(11)+0.01,Z(5)+0.01],[Y(2),Y(4),Y(11)],'c')
    drawLine(coords(:,1),coords(:,5),'-ob') %spring line
    drawLine(coords(:,3),coords(:,4),'-ob')
    drawLine( coords(:,6),coords(:,8),'-or')
    drawLine( coords(:,7),coords(:,9),'-or')
    drawLine( coords(:,12),coords(:,13),'-ob')
    drawLine( coords(:,11),coords(:,13),'-or')
    
    hold off
end


function drawLine(p1,p2,parameters)
    plot3( [p1(1),p2(1)] ,[p1(3),p2(3)],[p1(2),p2(2)],parameters,'LineWidth',2,'MarkerFaceColor','b' );
end
