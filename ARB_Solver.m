clc; clear;

data = readtable("CarData.xlsx",VariableNamingRule="preserve");

%% Sort Data

RG = 0.5 % deg/g , Desired Roll Gradient

W = data{1,4}; % car weight in lbs
W_distribution = data{1,5};
WUC = data{2,4}; %unsprung weight


%extract data
WF = data{1,2}; WR = data{1,3}; %weight on front/rear
WUF = WUC*W_distribution; WUR = WUC*(1-W_distribution); %unsprung weight
WSF = WF-WUF; WSR = WR-WUR; WS = WSF+WSR; %sprung weight


wheelBase = data{6,4};
trackF = data{7,2}; trackR = data{7,3}; %track width

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


%% Calculations (pg 603 in millitch)
F_freq = ( 1/ (2*pi)) * sqrt( (KRF * g_in) / (WSF/2) ); % front ride freq, not needed
R_freq = ( 1/ (2*pi)) * sqrt( (KRR * g_in) / (WSR/2) ); % rear ride freq

Sprung_distribution = WSF / (WS); % sprung mass weight distribution
% calc sprung mass CG height
hs = ( W*CG_height - WUF*SLR - WUR*SLR ) / (WS);

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
KARB = KRollReq - (KS) % lb*ft/deg