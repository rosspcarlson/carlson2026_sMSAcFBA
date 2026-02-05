function[mF, sF] = functionEc4800Paper1(i, zmax, rowtrack) 

%Ec4800AeroCost/NoCost
%December 2022-December 2023, March 2024, June/July 2024, March/June 2025  
%Ross P. Carlson
%cost analysis for aerobic cultures, glucose and iron limited (nitrogen
%separate due to saturation of transporter)
%analyzes median experimental chemostat data using surface area constraints
%using temperature corrected TO numbers
%corrects the transporter saturation term for the chemostat conditions

%SA:V Si et al 2017 MG1655 (fMSA 1-20%)
MREc = [1.4089;	2.8179;	4.2268;	5.6358; 7.0447;	8.4537;	9.8627;	11.2716;	12.6806;	14.08953;	15.4985;	16.9075;	18.3164;	19.7254;	21.1343;	22.5433;	23.9522;	25.3612;	26.7702;	28.1791]; %non growth parameters for (1-20%) MRE 
MREmu = [0.2883;	0.5766;	0.865;	1.1533; 1.4416;	1.73;	2.0183;	2.3066;	2.595;	2.8833;	3.1716;	3.46;	3.7483;	4.0366;	4.325;	4.6133;	4.9016;	5.19;	5.4783;	5.7666] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 MG1655 (fMSA 2-20%, 6.5%)
%MREc = [2.817; 4.226; 5.635; 7.044; 8.453; 9.158; 9.862; 11.271; 12.680; 14.089; 15.498; 16.907; 18.316; 19.725; 21.134; 22.543; 23.952; 25.361; 26.770; 28.179];
%MREmu = [0.576; 0.865; 1.153; 1.441; 1.73; 1.874; 2.0183; 2.306; 2.595; 2.883; 3.1716; 3.46; 3.748; 4.0366; 4.325; 4.613; 4.901; 5.19; 5.478; 5.766] ; %growth associated parameters for (2-20, 7.5%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 1-20%)
%MREc = [1.9370;	3.8741;	5.8112;	7.7483; 9.6854;	11.6225;	13.5595;	15.4966;	17.4337;	19.3708;	21.3079;	23.245;	25.1820;	27.1191;	29.0562;	30.9933;	32.9304;	34.8675;	36.8045;	38.7416]; %non growth parameters for (1-20%) MRE 
%MREmu = [0.5508;	1.1016;	1.6525;	2.2033; 2.7541;	3.305;	3.8558;	4.4066;	4.9575;	5.5083;	6.0591;	6.61;	7.1608;	7.7116;	8.2625;	8.8133;	9.3641;	9.915;	10.4658;	11.0166] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 2-20%, 6.5%)
%MREc = [3.874;	5.811;	7.7483;	9.685;	11.622;	12.591;	13.559;	15.496;	17.433;	19.370;	21.307; 23.245;	25.182; 27.119; 29.056; 30.993; 32.930; 34.867; 36.804; 38.741]; %non growth parameters for (1-20%) MRE 
%MREmu = [1.101;    1.652; 2.203; 2.754; 3.305; 3.580; 3.855; 4.406; 4.957;	5.508; 6.059; 6.61; 7.160; 7.711; 8.262; 8.813; 9.364; 9.915; 10.465; 11.0166] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 1-20%)SA:V +10%
%MREc = [2.130;	4.261;	6.392;	8.523;	10.654;	12.785;	14.915;	17.046;	19.177;	21.308;	23.439;	25.57;	27.700;	29.831;	31.962;	34.093;	36.224;	38.355;	40.485;	42.616]; %non growth parameters for (1-20%) MRE 
%MREmu = [0.605;	1.211;	1.817;	2.423;	3.029;	3.635;	4.241;	4.847;	5.453;	6.059;	6.665;	7.271;	7.876;	8.482;	9.088;	9.694;	10.300;	10.906;	11.512;	12.118] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 1-20%)SA:V +25%
%MREc = [2.4214;	4.842;	7.264;	9.6858;	12.107;	14.528;	16.950;	19.371;	21.793;	24.214;	26.636;	29.057;	31.478;	33.900;	36.321;	38.743;	41.164;	43.586;	46.007;	48.429]; %non growth parameters for (1-20%) MRE 
%MREmu = [0.688;	1.377;	2.065;	2.754;	3.442;	4.131;	4.819;	5.5083;	6.1968;	6.885;	7.5739;	8.262;	8.951;	9.639;	10.328;	11.0166;	11.705;	12.393;	13.0822;	13.7708] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 1-20%)SA:V +50%
%MREc = [2.905625;	5.81125;	8.716875;	11.6225;	14.528125;	17.43375;	20.339375;	23.245;	26.150625;	29.05625;	31.961875;	34.8675;	37.773125;	40.67875;	43.584375;	46.49;	49.395625;	52.30125;	55.206875;	58.1125]; %non growth parameters for (1-20%) MRE 
%MREmu = [0.82625; 1.6525; 2.47875; 3.305; 4.13125; 4.9575; 5.78375; 6.61; 7.43625; 8.2625; 9.08875; 9.915; 10.74125; 11.5675; 12.39375; 13.22; 14.04625; 14.8725; 15.69875; 16.525] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 1-20%)SA:V -10%
%MREc = [1.743;	3.486;	5.230;	6.973;	8.716;	10.460;	12.203;	13.947;	15.690;	17.433;	19.177;	20.920;	22.663;	24.407;	26.150;	27.894;	29.637;	31.380;	33.124;	34.867]; %non growth parameters for (1-20%) MRE 
%MREmu = [0.495;	0.991;	1.487;	1.983;	2.478;	2.974;	3.470;	3.966;	4.461;	4.957;	5.453;	5.949;	6.444;	6.940;	7.436;	7.932;	8.427;	8.923;	9.419;	9.915] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 1-20%)SA:V -25%
%MREc = [0.968;	1.937;	2.905;	3.874;	4.842;	5.811;	6.779;	7.748;	8.716;	9.685;	10.653;	11.622;	12.591;	13.559;	14.528;	15.496;	16.465;	17.433;	18.402;	19.370]; %non growth parameters for (1-20%) MRE 
%MREmu = [0.275; 0.550; 0.826; 1.101; 1.377; 1.652; 1.927; 2.203; 2.478; 2.754; 3.029; 3.305; 3.580; 3.855; 4.13125; 4.406; 4.682; 4.957; 5.232; 5.508] ; %growth associated parameters for (1-20%) MRE

%SA:V Si et al 2017 NCM3722 (fMSA 1-20%)SA:V -50%
%MREc = [1.452;	2.905;	4.358;	5.811;	7.264;	8.716;	10.169;	11.622;	13.0753;	14.528;	15.980;	17.433;	18.886;	20.339;	21.792;	23.245;	24.697;	26.150;	27.603;	29.056]; %non growth parameters for (1-20%) MRE 
%MREmu = [0.413;	0.826;	1.239;	1.652;	2.065;	2.478;	2.891;	3.305;	3.718;	4.131;	4.544;	4.957;	5.370;	5.783;	6.196;	6.61;	7.0231;	7.436;	7.849;	8.262] ; %growth associated parameters for (1-20%) MRE

% limiting substrate transporter saturation for growth rates in data set,
% PtsG glucose transporter saturation, wild type phenotypes both strains
%satnGlc = [ 0.149	0.223	0.298	0.373	0.447	0.497	0.546	0.597	0.646	0.695	0.746	0.795	0.844	0.895	0.944	0.994	1	1	1	1	1	1	1   0.1	0.15	0.2	0.25	0.3	0.35	0.4	0.45	0.5	0.55	0.6	0.65	0.7	0.75	0.8	0.85	0.9	0.95	1	1	1	1	1];
%PtsG glucose transporter saturation,MG1655 wildtype, NCM3722 saturated ptsG, recombinant expression
satnGlc = [ 0.149	0.223	0.298	0.373	0.447	0.497	0.546	0.597	0.646	0.695	0.746	0.795	0.844	0.895	0.944	0.994	1	1	1	1	1	1	1   1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1];
%Amt NH4+ transporter saturation
satnNH3 = [ 1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1.00	1.00	1.00	1.00	1.00	1.00   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 1 1 1 1 1 1]; 

%experimental data (order of columns: growth
%rate, glucose, pyruvate, succinate, lactate out, lactate in, formate,
%acetate out, acetate in, ethanol out)

%rows: 1-23, cstr/batch, glucose limited, consensus data MG1655
%rows: 24-46, cstr/batch, glucose limited, consensus data NCM3722

RM = [0.100	1.566	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
0.150	2.223	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
0.200	2.880	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
0.250	3.537	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
0.300	4.194	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
0.333	4.628	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000
0.366	5.061	0.000	0.000	0.000	0.000	0.000	0.202	0.000	0.000
0.400	5.508	0.000	0.000	0.000	0.000	0.000	0.527	0.000	0.000
0.433	5.942	0.000	0.000	0.000	0.000	0.000	0.843	0.000	0.000
0.466	6.375	0.000	0.000	0.000	0.000	0.000	1.159	0.000	0.000
0.500	6.822	0.000	0.000	0.000	0.000	0.000	1.484	0.000	0.000
0.533	7.256	0.000	0.000	0.000	0.000	0.000	1.800	0.000	0.000
0.566	7.689	0.000	0.000	0.000	0.000	0.000	2.115	0.000	0.000
0.600	8.136	0.000	0.000	0.000	0.000	0.000	2.441	0.000	0.000
0.633	8.570	0.000	0.000	0.000	0.000	0.000	2.757	0.000	0.000
0.666	9.003	0.000	0.000	0.000	0.000	0.000	3.072	0.000	0.000
0.700	9.450	0.000	0.000	0.000	0.000	0.000	3.398	0.000	0.000
0.733	9.884	0.000	0.000	0.000	0.000	0.000	3.713	0.000	0.000
0.766	10.317	0.000	0.000	0.000	0.000	0.000	4.029	0.000	0.000
0.800	10.764	0.000	0.000	0.000	0.000	0.000	4.354	0.000	0.000
0.833	11.198	0.000	0.000	0.000	0.000	0.000	4.670	0.000	0.000
0.866	11.631	0.000	0.000	0.000	0.000	0.000	4.986	0.000	0.000
0.900	12.078	0.000	0.000	0.000	0.000	0.000	5.311	0.000	0.000
0.1	1.50	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.15	2.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.2	2.50	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.25	3.02	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.3	3.55	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.35	4.09	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.4	4.63	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.45	5.18	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.5	5.73	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.55	6.26	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.60	6.80	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0
0.65	7.40	0.00	0.00	0.00	0.00	0.00	0.13	0.00	0
0.70	8.00	0.00	0.00	0.00	0.00	0.00	0.25	0.00	0
0.75	8.77	0.00	0.00	0.00	0.00	0.00	0.80	0.00	0
0.80	9.35	0.00	0.00	0.00	0.00	0.00	1.40	0.00	0
0.85	10.40	0.00	0.00	0.00	0.00	0.00	3.00	0.00	0
0.90	11.40	0.00	0.00	0.00	0.00	0.00	5.00	0.00	0
0.95	12.70	0.00	0.00	0.00	0.00	0.00	7.30	0.00	0
1.00	14.00	0.00	0.00	0.00	0.00	0.00	9.30	0.00	0
1.05	15.43	0.00	0.00	0.00	0.00	0.00	11.78	0.00	0
1.10	17.06	0.00	0.00	0.00	0.00	0.00	14.34	0.00	0
1.15	18.87	0.00	0.00	0.00	0.00	0.00	17.06	0.00	0
1.20	20.86	0.00	0.00	0.00	0.00	0.00	19.94	0.00	0];

%error in experimental data, same order as original rates  
RMEr = [0.01	0.11	0	0	0	0	0	0.00	0	0.000
0.01	0.16	0	0	0	0	0	0.00	0	0.000
0.01	0.20	0	0	0	0	0	0.00	0	0.000
0.01	0.25	0	0	0	0	0	0.00	0	0.000
0.02	0.29	0	0	0	0	0	0.00	0	0.000
0.02	0.32	0	0	0	0	0	0.00	0	0.000
0.02	0.35	0	0	0	0	0	0.00	0	0.000
0.02	0.39	0	0	0	0	0	0.00	0	0.000
0.02	0.42	0	0	0	0	0	0.06	0	0.000
0.02	0.45	0	0	0	0	0	0.08	0	0.000
0.02	0.48	0	0	0	0	0	0.10	0	0.000
0.02	0.51	0	0	0	0	0	0.13	0	0.000
0.02	0.54	0	0	0	0	0	0.15	0	0.000
0.02	0.57	0	0	0	0	0	0.17	0	0.000
0.02	0.60	0	0	0	0	0	0.19	0	0.000
0.02	0.63	0	0	0	0	0	0.22	0	0.000
0.02	0.66	0	0	0	0	0	0.24	0	0.000
0.02	0.69	0	0	0	0	0	0.26	0	0.000
0.02	0.72	0	0	0	0	0	0.28	0	0.000
0.02	0.75	0	0	0	0	0	0.30	0	0.000
0.02	0.78	0	0	0	0	0	0.33	0	0.000
0.02	0.81	0	0	0	0	0	0.35	0	0.000
0.02	0.85	0	0	0	0	0	0.37	0	0.000
0.01	0.11	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.14	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.18	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.21	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.25	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.29	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.32	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.36	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.40	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.44	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.48	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.52	0	0	0	0.00	0.00	0.00	0.00	0
0.01	0.56	0	0	0	0.00	0.00	0.00	0.00	0
0.02	0.61	0	0	0	0.00	0.00	0.00	0.00	0
0.02	0.65	0	0	0	0.00	0.00	0.10	0.00	0
0.02	0.73	0	0	0	0.00	0.00	0.21	0.00	0
0.02	0.80	0	0	0	0.00	0.00	0.35	0.00	0
0.02	0.89	0	0	0	0.00	0.00	0.51	0.00	0
0.02	0.98	0	0	0	0.00	0.00	0.65	0.00	0
0.02	1.08	0	0	0	0.00	0.00	0.82	0.00	0
0.02	1.19	0	0	0	0.00	0.00	1.00	0.00	0
0.02	1.32	0	0	0	0.00	0.00	1.19	0.00	0
0.02	1.46	0	0	0	0.00	0.00	1.40	0.00	0];

totno = 100 ; %100 perturbations
rxnno = 91; %reactions in model

%vector to be pulled from function

mF = zeros(rowtrack, zmax); %median flux for each reaction or function
sF = zeros(rowtrack, zmax); %standard deviation for each reaction or function

%number of experimental conditions
LR = size(RM,1);
%number of experimental rates considered
LC = size(RM,2);

z = 1 ; %starting MSA % 

while z < zmax +1  %cycles through MSA %, set with z and zmax

%output
MREout = zeros(rxnno + 6,totno);
 
loop = 1; 

%temporary RM with perturbations
RMtemp = zeros(LR, LC);

%one hundred perturbations loop
while loop <totno + 1  

%module for perturbing the experimental values within their limits

x = 1; % row counter
y = 1;  % column counter

while x < LR+1   
    
    y = 1 ;
    while y < LC+1 
        
        RMtemp(x,y) = RM(x,y) + (randi([-1,1],1,1))*RMEr(x,y)*rand; % perturbation routine, normal 
        %RMtemp(x,y) = RM(x,y) +  RMEr(x,y)*(1 + (-1-1).*rand(1,1)); % perturbation routine, uniform
        y = y +1 ;
         
    end
    
    x = x + 1 ;
    
end

FBAfunction4800Paper1('Models/EcMREpaper1.mat', 'ATPm',RMtemp(i,1),RMtemp(i,2),RMtemp(i,3), RMtemp(i,4), RMtemp(i,5), RMtemp(i,6), RMtemp(i,7), RMtemp(i,8), RMtemp(i,9),RMtemp(i,10),MREc(z), MREmu(z), satnGlc(i),satnNH3(i));

if isempty(ans)
    MREout(1:rxnno,loop) = 0; %fluxes output
    MREout(rxnno+1, loop) = 0; %PO number
    MREout(rxnno+2, loop) = 0; %ATP synthase
    MREout(rxnno+3,loop) = 0 ;% glycolysis split
    MREout(rxnno+4,loop) = 0; %NADPH from PPP
    MREout(rxnno+5,loop) = 0 ; %NADPH from PntAB
    MREout(rxnno+6,loop) = 0; %NADPH from Icd
       
else
    
    MREout(1:rxnno,loop) = ans ; %fluxes output
    MREout(rxnno+1,loop) = ans(45)/(ans(29) + ans(39) + ans(42) + ans(43) + ans(45) + ans(46)) + (0.5*ans(47))/(ans(47) + ans(48) + ans(49)) + ans(49)/(ans(47) + ans(48) +  ans(49)); %PO
    MREout(rxnno+2,loop) = (ans(50)+abs(ans(51)))*118.56/(3.7864*exp(ans(64)*0.7946)); %ATP synthase enzymes per nm^2
    MREout(rxnno+3,loop) = ans(4)/(ans(4) + ans(14)) ;% glycolysis split
    MREout(rxnno+4,loop) = ans(14)/(ans(14)+ans(26)+ans(53)); %NADPH from PPP
    MREout(rxnno+5,loop) = ans(53)/(ans(14)+ans(26)+ans(53)) ; %NADPH from PntAB
    MREout(rxnno+6,loop) = ans(26)/(ans(14)+ans(26)+ans(53)); %NADPH from Icd
      
end

loop = loop + 1 ;

end
 
%Excel output file information

Descr = ["250602Ec4800NCM71_" + i + ".xlsx"];    %***

file = Descr; 
%rang = 'A2:CV93';                        
%coln = 'A1:BF1';                         

if z== 1  
         
        tab = 'MSA1%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA1_empty ', i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
            
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end
if z== 2  
         
        tab = 'MSA2%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA2_empty ', i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end

if z== 3  
         
        tab = 'MSA3%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA3_empty ', i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end

 if z== 4  
         
        tab = 'MSA4%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA4_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
           mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end
if z== 5  
         
        tab = 'MSA5%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA5_empty ', i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
            
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end
if z== 6  
         
        tab = 'MSA6%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA6_empty ', i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end

if z== 7  
         
        tab = 'MSA7%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA7_empty ', i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end

 if z== 8  
         
        tab = 'MSA8%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA8_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
           mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
 end
   
 if z== 9  
         
        tab = 'MSA9%'; 
        MREout(:,all(MREout ==0)) = [];
        if isempty(MREout)
            'MSA9_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end   
 end
if z== 10  
         
       tab = 'MSA10%'; 
       MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
       if isempty(MREout)
            'MSA10_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
       else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
end
if z== 11  
         
      tab = 'MSA11%'; 
      MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
      if isempty(MREout)
            'MSA11_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
      else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
end
if z== 12  
         
        tab = 'MSA12%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA12_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
end
if z== 13  
         
        tab = 'MSA13%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA13_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
end
if z== 14  
         
        tab = 'MSA14%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA14_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            xlswrite(file,MREout,tab);
        end 
end
if z== 15  
        
        tab = 'MSA15%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
        if isempty(MREout)
            'MSA15_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
        else
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            xlswrite(file,MREout,tab);
        end 
end     
if z== 16  
         
        tab = 'MSA16%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
       if isempty(MREout)
            'MSA16_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
       else
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            xlswrite(file,MREout,tab);
       end 
       end      

if z== 17  
         
        tab = 'MSA17%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
       if isempty(MREout)
            'MSA17_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
       else
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            xlswrite(file,MREout,tab);
       end 
       end      

if z== 18  
         
        tab = 'MSA18%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
       if isempty(MREout)
            'MSA18_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
       else
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            xlswrite(file,MREout,tab);
       end 
       end      

if z== 19  
         
        tab = 'MSA19%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
       if isempty(MREout)
            'MSA19_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
       else
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            xlswrite(file,MREout,tab);
       end 
end      

if z== 20  
         
        tab = 'MSA20%'; 
        MREout(:,all(MREout ==0)) = [];  %sorts empty data sets
       if isempty(MREout)
            'MSA20_empty ',i
            mF(:,z) = zeros(rowtrack, 1);
            sF(:,z) = zeros(rowtrack, 1);
       else
            MREout = [MREout (mean(MREout, 2))  (std(MREout,0, 2))]; %adds columns with reaction flux mean and standard deviation
            mF(:,z) = (mean(MREout, 2));
            sF(:,z) = (std(MREout,0, 2));
            xlswrite(file,MREout,tab);
       end 
end      
     z = z + 1 ;
end   
end
 
