clc; clear

%% Suspension Values *EDIT M VARIABLES*
MSpring = [10.1;15.11;0]; % Spring mount to chassis X Y Z in.
MBell = [11.5;7.92;0]; %Bellcrank pivot point
MPush = [21.72;4.21;0];

MBellSpr = [15.151; 10.263;0]; %Point where bellcrank connects to spring
MBellPush = [14.566; 9.294;0]; % bellcrank connects to rod

MArmL = [10.35 ; 2.75 ; 0]; % Where lower arm connects to chassis
MArmU = [12 ; 7.67 ; -0]; %upper arm to chassis

MArmLEnd = MPush; % Where knuckle connects to lower arm
MArmUEnd = [21.72 ; 6.62 ; -0]; %Where knuckle connects to upper arm
MknuckleDistance = 0.97;


v_knuckle = (MArmUEnd-MArmLEnd)/norm(MArmUEnd-MArmLEnd);
MWheel = v_knuckle*(MknuckleDistance) + MArmLEnd;
points = [MSpring,MBell,MPush,MBellPush,MBellSpr, ... 
    MArmL, MArmU, MArmLEnd , MArmUEnd,MWheel];

%% fixed Parameters
springLength = norm(MSpring - MBellSpr);
pushrodLength = norm(MPush - MBellPush);
lowerArmLength = norm(MArmL - MArmLEnd);
upperArmLength = norm(MArmU - MArmUEnd);
bellSprLength = norm(MBellSpr - MBell);
bellPushLength = norm(MBellPush - MPush);

fixedParam = [springLength,pushrodLength,lowerArmLength,upperArmLength]

%% Solver function

function F = VectorEquations(points)
    %unpack variables
    [MSpring,MBell,MPush,MBellPush,MBellSpr, ... 
        MArmL, MArmU, MArmLEnd , MArmUEnd,MWheel] = points;

    % Calculate fixed link lengths and enforce geometric constraints

    F = [norm(MBell - MBellSpr) - springLength;
        norm(MBell - MBellPush) - pushrodLength;
        norm(MArmL - MArmLEnd) - lowerArmLength;
        norm(MArmU - MArmUEnd) - upperArmLength];
end