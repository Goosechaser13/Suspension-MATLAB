
W = 580 % car weight in lbs


FR_distribution = .45 %distribution of front/total weight

WF = 580 * FR_distribution; %weight on front
WR = W - WF; % weight on rear


kSpring = 150 %lbs / in
MR = 1 %motion ratio
g = 386.4 % gravity in inch/s^2