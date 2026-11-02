clc; clear;

data = readtable("CarData.xlsx",VariableNamingRule="preserve");

%extract data
W = data{1,4}; % car weight in lbs
W_distribution = data{1,5};
WUC = data{2,4}; %unsprung weight
WF = data{1,2}; WR = data{1,3}; %weight on front/rear
WUF = WUC*W_distribution; WUR = WUC*(1-W_distribution); %unsprung weight
WSF = WF-WUF; WSR = WR-WUR; WS = WSF+WSR; %sprung weight

FLTD = data{14,2}; %Front Load Transfer distribution
RLTD = data{14,3}; %Rear Load Transfer distribution

wheelBase = data{6,4};
trackF = data{7,2}; trackR = data{7,3}; %track width in.

kSpring = data{8,4}; %lbs / in
MRF = data{5,2}; MRR = data{5,3}; % motion ratios
KRF = kSpring*(MRF^2); KRR = kSpring*(MRR^2); %lbs/in wheel rates

ZF = data{10,2} ; ZR = data{10,3};

TireRate = data{12,4}; % tire rate
TireRadius = data{13,4};
SLR = TireRadius*0.9; % tire static loaded radius (assuming 90% of radius)
CG_height = data{9,4}; %in

g = 32.174; % gravity in feet
g_in = 386.4; % gravity in inch/s^2

cachedData = {W,W_distribution,WUC,WF,WR,WUF,WUR,WSF,WSR,WS,FLTD,RLTD,wheelBase,trackF,trackR,kSpring...
    MRF,KRF,KRR,ZF,ZR,TireRate,TireRadius,SLR,CG_height,g,g_in}

%% Calculating Bar stiffness
function [KF_bar,KR_bar] = RollToBarK(RG,IF,IR,LF,LR,cachedData)
    trackF = cachedData{14}; trackR = cachedData{15};
    [KARBF,KARBR] = RollCalculator(RG,cachedData);
    %IF bar/wheel_center movement, MR for roll bar
    %IR
    %LF %ARB lever arm length in
    %LR %in
    KF_bar = (KARBF / (IF^2))* ( (LF/12) ^2/ (trackF/12) ^2); % torsion bar rate lb*in/deg
    KR_bar = (KARBR / (IR^2))* ( (LR/12) ^2/ (trackR/12) ^2);
end
RG = 0.5
[KF_bar,KR_bar] = RollToBarK(RG,1,1,3,3,cachedData);
fprintf("Assuming LLTD: F: %.1f%% R: %.1f%% , Roll gradient: %.2f [deg/g] \n",FLTD*100,RLTD*100,RG)
fprintf("Front Bar stiffness = %.3f [lb*ft/deg] \n",KF_bar);
fprintf("Rear Bar stiffness = %.3f [lb*ft/deg] \n",KR_bar);


%% Plot Data
RG = 0.2:0.05:2; % deg/g , Desired Roll Gradient
stiffnessForRG = zeros([2,length(RG)]);
for i = 1:length(RG)
    [KARBF,KARBR] = RollCalculator(RG(i),cachedData);
    stiffnessForRG(1,i) = KARBF;
    stiffnessForRG(2,i) = KARBR;
end

figure(1);
hold on 
plot(RG,stiffnessForRG(1,:),'DisplayName','Front Stiffness');
plot(RG,stiffnessForRG(2,:),'DisplayName','Rear Stiffness');
title("Stiffness for desired Roll Gradient");
legend();
xlabel('Roll Gradient [deg/g]');
ylabel('Stiffness [lb*ft/deg]');
theme(gcf,'light');

function [KARB,KARBR] = RollCalculator(RG,cachedData)
    %% Calculations (pg 603 in millitch)
    [W,W_distribution,WUC,WF,WR,WUF,WUR,WSF,WSR,WS,FLTD,RLTD,wheelBase,trackF,trackR,kSpring...
        MRF,KRF,KRR,ZF,ZR,TireRate,TireRadius,SLR,CG_height,g,g_in] = deal(cachedData{:});

    F_freq = ( 1/ (2*pi)) * sqrt( (KRF * g_in) / (WSF/2) ); % front ride freq, not needed
    R_freq = ( 1/ (2*pi)) * sqrt( (KRR * g_in) / (WSR/2) ); % rear ride freq
    
    
    % calc sprung mass CG height
    hs = ( W*CG_height - WUF*SLR - WUR*SLR ) / (WS);
    
    Sprung_distribution = WSF / (WS); % sprung mass weight distribution
    
    % rolling moment lever arm, the distance from sprung CG and roll axis
    hRM = hs - (ZF + (ZR - ZF)*(1-Sprung_distribution));
    
    % now calculate the roll moment (lbs*ft per lateral g) this creates
    MRG = (hRM*WS) / 12; %divide by 12 to convert to feet for gs
    
    % Calculated Roll rate required, using roll gradient:
    KRollReq = MRG / RG; % lb*ft / deg
    
    % Roll rate provided by springs:
    KSF = (KRF * trackF^2) / (4320/pi); %units: lb*ft/deg
    KSR = (KRR * trackR^2) / (4320/pi); %units: lb*ft/deg
    KS = KSF + KSR; % roll rate available from springs
    
    % ROLL RATE REQUIRED from ARBs:
    KARB = KRollReq - (KS); % lb*ft/deg
    
    
    totalLoadTrans = W*CG_height / ((trackF+trackR) /2 ); % per lateral g
    FLT = totalLoadTrans*FLTD; %lb/g
    RLT = totalLoadTrans*RLTD; %lb/g
    
    %Calculate roll stiffness from pg 604
    KF = (FLT - (WSF*ZF/trackF) - (WUF*SLR/trackF) ) / (12*RG/trackF ); %lb*ft/deg
    KR = KRollReq - KF; %lb*ft/deg rear stiffness
    
    KARBF = KF - KSF; % stiffness contributiuon to total car roll rate
    KARBR = KRollReq - KF - KSR;

end

