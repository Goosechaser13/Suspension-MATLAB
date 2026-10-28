
w_axle = 580 %lbs
kSpring = 150 %lbs / in
MR = 1 %motion ratio
g = 386.4 % gravity in inch/s^2

m_axle = w_axle / g %turn weight to mass 
kWheel = kSpring * (MR)^2

rideFreq = sqrt(kSpring/m) %rads/s
Kw = (k)