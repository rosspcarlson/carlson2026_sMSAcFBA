
%EcMREjobGLA2500
%October 2022 - May 2023 
%Modified to account for PtsG saturation, May 2024
%Ross P. Carlson
% Calculates phenotype constrained by membrane surface area, fixes two
% resource demands mu and qmATP and minimizes substrate (glucose) flux to achieve
% requirements, no membrane surface area constraint for secretion
%modified for TO analysis (temp corrected TO parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%initial values

k = 1; % overall counter per z  
y = 1; % counter for mu loop
z = 1; % counter for MRE values
zmax = 15 ; % counter for MRE percentages
MWglc = 0.18; %molecular weight glucose, mmol basis, units g     
MWlac = 0.09; %molecular weight lactate, mmol basis, units g     
MWace = 0.06; %molecular weight acetate, mmol basis, units g      
subglc = 1 ; % reaction glucose 
sublac = 69 ; % reaction lactate in  
subace = 66 ; % reaction acetate in  
biom = 63; % model row for biomass growth

%SA:V Si et al 2017
MREc = [1.408958333;	2.817916667;	4.226875;	5.635833333;	7.044791667;	8.45375;	9.862708333;	11.27166667;	12.680625;	14.08958333;	15.49854167;	16.9075;	18.31645833;	19.72541667;	21.134375]; %non growth parameters for (1-15%) MRE 
MREmu = [0.288333333	;0.576666667	;0.865	;1.153333333	;1.441666667	;1.73	;2.018333333	;2.306666667	;2.595	;2.883333333	;3.171666667	;3.46	;3.748333333	;4.036666667;	4.325] ; %growth associated parameters for (1-15%) MRE

ModRxnList = ["Rx1 G1 ";   "Rx2 G2a ";    "Rx3  G2b "; "Rx4  G3r "; "Rx5  G4";   "Rx6  G5";   "Rx7  G6r";  "Rx8  G7r";  "Rx9  G8r";  "Rx10  G9r"; "Rx11  G10r";    "Rx12  G11r";    "Rx13  G12r";    "Rx14  PPP1r";   "Rx15  PPP2r";   "Rx16  PPP3";    "Rx17  PPP4";    "Rx18  PPP5";    "Rx19  PPP6r";   "Rx20  PPP7r";   "Rx21  PPP8r";...
    "Rx22  PPP9r ";  "Rx23  PPP10r "; "Rx24  TCA1";    "Rx25  TCA2r";   "Rx26:  TCA3";   "Rx27  TCA4";    "Rx28  TCA5r ";  "Rx29  TCA6";    "Rx30  TCA7";    "Rx31  TCA8r";   "Rx32  TCA9r";   "Rx33  AP1"; "Rx34  AP2"; "Rx35  AP3"; "Rx36  IM1"; "Rx37  IM2"; "Rx38  IM3"; "Rx39  IM3b";    "Rx40  IM4"; "Rx41  IM5r";    "Rx42  IM6 ";...
    "Rx43  IM7"; "R44  IM8"; "Rx5  R1";  "Rx46  R2";  "Rx47  R3";  "Rx48  R4 "; "Rx49  R5";  "Rx50  R6r (ATPase)";    "Rx51  R7";  "Rx52  R8";  "Rx53  GS1r ";   "Rx54  GS2"; "Rx55  AA1"; "Rx56  AA2"; "Rx57  AA3r ";   "Rx58  ATPm";   "Rx59  ADK ";    "Rx60  PYRO ";   "Rx61  PH";  "Rx62  Phin";    "Rx63  BIOM ";   "Rx64  T1 (Eto)";...
    "Rx65  T2 (Ace)";    "Rx66  T2in";    "Rx67  T3 (NH3)";    "Rx68 T4 (Lac)"; "Rx69  T4in";    "Rx70  T5 (For)";    "Rx71  T6 (CO2)";    "Rx72  T6in ";   "Rx73  T7 (O2)"; "Rx74  T8";  "Rx75  T9 (H2O)";    "Rx76  T9in";    "Rx77  T10"; "Rx78  T11 (Pyr) ";  "Rx79  T12 (Mns)";   "Rx80  T13 (Xyl)";   "Rx81  T14 (Gly)";   "Rx82  T15 (Suc)";   "Rx83  T16 ";    "Rx84  T17 ";...
    "Rx85  T18"; "Rx86  T19"; "Rx87  T20"; "Rx88  T21"; "Rx89  T22"; "Rx90  T23 ";    "Rx91  T24"; "Rx92  T25 (MRE)"; "Rx93 T26 (Emol)"; "Rx94 T27 (Cmol)"; "Rx95 T28 (Mol)"; "Yx/s"; "P/O number"; "1/Yx/s"; "MRE/qATP"];

s = 95; %number of reactions in EcMREModel    

%specific rates for constraints

vATPm = [0 1 3 5 7 9	11	13	15 17 19 21	23 25 27	29 31 33	35 37 39	41	43 45 47	49 51 53 55	57 59 61 63	65	67 69	71 73	75	77	79	81 83 85	87 89 91 93 95 97 99 101 103 105 107 109 111 113 115 117 119 121];

n = size(vATPm, 2) ;% counter for number of increases in mATP

Output = zeros(s+2,(n)*(n+35)); %preallocate output matrix

while z < zmax +1 %cycles through MRE percentages 

  Output = zeros(s+2,(n)*(n+35)); %preallocate output matrix  
  mu = 0.1 ;%  initial growth rate
  k = 1 ;
  y = 1 ;
  
while y < n + 73 %cycles through mu 

i = 1;

while i < n %cycles through mATP

    qatp = vATPm(i) ;

FBAfunctionMREGlucose2500Q10('Models/EcMREnocostv2.mat','G1',qatp, mu, MREc(z), MREmu(z));  %   ***

if isempty(ans)
    Output(1:s,k) = 0;
    Output(s+1,k) = 0 ;
    Output(s+2,k) = 0 ;
else
    Output(1:s,k) = ans ;
    %biomass yield on substrate
    Output(s+1, k) = Output(biom,k)/((Output(subglc,k))*MWglc + (Output(sublac,k))*MWlac + (Output(subace,k))*MWace); %biomass yield on substrate
    %P/O number
    Output(s+2, k) = Output(45,k)/(Output(29,k) + Output(39,k) + Output(42,k) + Output(43,k) + Output(45,k) + Output(46,k) + Output(81,k)) + (0.5*Output(47,k))/(Output(47,k) + Output(48,k) + Output(49,k)) + Output(49,k)/(Output(47,k) + Output(48,k) +  Output(49,k)); %PO
    %gram substrate per gram biomass
    Output(s+3, k) = 1/(Output(s+1,k));
    %MRE per qATP
    Output(s+4, k) = (Output(92,k)-Output(63,k)*MREmu(z))/(Output(58,k));

end

k = k + 1;
i = i+1 ;

end

mu = mu + 0.01 ;

y = y + 1 ;
end

Descr = ['240602EcMRE_2500c', '.xlsx'];    %***
file = Descr; 
rang = 'B1:JIF99';                        
coln = 'A1:A99';                         

if z== 1  
        %Output1 = Output;
        tab = 'MRE1'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file,ModRxnList(:),tab,coln) ;   
        EucD1 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
   
 if z== 2  
        %Output2 = Output;
        tab = 'MRE2'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file,ModRxnList(:),tab,coln) ; 
        EucD2 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
 end
if z== 3  
        %Output3 = Output;
        tab = 'MRE3'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD3 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 4  
        %Output4 = Output;
        tab = 'MRE4'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD4 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 5  
        %Output5 = Output;
        tab = 'MRE5'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD5 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 6  
        %Output6 = Output;
        tab = 'MRE6'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ;
        EucD6 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 7  
        %Output7 = Output;
        tab = 'MRE7'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD7 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 8  
        %Output8 = Output;
        tab = 'MRE8';
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD8 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 9  
        %Output9 = Output;
        tab = 'MRE9'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD9 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 10  
        %Output10 = Output;
        tab = 'MRE10'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD10 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end

if z== 11  
        %Output11 = Output;
        tab = 'MRE11'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD11 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 12  
        %Output12 = Output;
        tab = 'MRE12'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD12 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 13  
        %Output13 = Output;
        tab = 'MRE13'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD13 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 14  
        %Output14 = Output;
        tab = 'MRE14'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD14 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
if z== 15  
        %Output15 = Output;
        tab = 'MRE15'; 
        Output(:,all(Output ==0)) = [];
        xlswrite(file,Output,tab,rang);
        xlswrite(file, ModRxnList(:),tab,coln) ; 
        EucD15 = [Output(1,:);Output(58,:); Output(63,:);Output(65,:);Output(96,:);Output(97,:)];
end
z = z + 1 ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%