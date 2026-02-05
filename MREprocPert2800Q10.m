%MREprocPert2800
%December 2022-May 2023
%Ross P. Carlson
%solves FBA problem maximizing ATP production given experimental flux
%constraints, experimental fluxes are perturbed within their confidence
%interval and the MRE enzyme turnover numbers are perturbed +/-5 to 95% to test
%robustness to assumed parameters

%no cost overflow metabolism
%MRE = 7% of surface area is protein
%examines the median data set for each condition, batch data only

%columns, specific rates (mmol/g cdw/h):  
%batch data (mu(1), gluc(2), pyr(3), suc(4), lactout(5), lactin(6) for(7), aceout(8), acein(9))

clear

pert = 0.95; % percent perturbation, change for each analysis
P = 7 ; %membrane surface area protein percentage

%SA:V Volkmer and Heinemann 2011 (3-29%,2%)
%MREmu = [0.391458333;	1.174375;	1.957291667;	2.740208333;	3.523125;	3.914583333;	4.306041667;	4.6975;	5.088958333;	5.871875;	7.04625;	8.22062;	9.78645833;	10.9608333;	-12.13520833]; %non growth parameters for MRE (3-29%)
%MREc = [1.020833333;	3.0625;	5.104166667;	7.145833333;	9.1875;	10.20833333;	11.22916667;	12.25;	13.27083333;	15.3125;	18.375;	21.4375;	25.52083333;	28.58333333;	31.64583333] ; %growth associated parameters for MRE (3-29%)
%SA:V Si et al 2017
MREc = [1.408958333;	2.817916667;	4.226875;	5.635833333;	7.044791667;	8.45375;	9.862708333;	11.27166667;	12.680625;	14.08958333;	15.49854167;	16.9075;	18.31645833;	19.72541667;	21.134375]; %non growth parameters for (1-15%) MRE 
MREmu = [0.288333333	;0.576666667	;0.865	;1.153333333	;1.441666667	;1.73	;2.018333333	;2.306666667	;2.595	;2.883333333	;3.171666667	;3.46	;3.748333333	;4.036666667;	4.325] ; %growth associated parameters for (1-15%) MRE

%median batch data sets from the MRE analysis, see output routines for
%condition order

RM = [0.65	9.10	0.00	0.00	0.00	0.00	0.00	2.83	0.00
0.71	7.31	0.00	0.00	0.00	4.98	0.00	4.02	0.00
0.68	9.29	0.00	0.00	0.00	0.00	0.00	2.72	0.00
0.62	8.19	0.00	0.00	0.00	3.45	0.00	5.34	0.00
0.38	0.00	0.00	0.00	0.00	17.69	0.00	7.76	0.00
0.36	0.00	0.00	0.00	0.00	17.12	0.00	7.13	0.00
0.22	0.00	0.00	0.00	0.00	0.00	0.00	0.00	9.92
];

%error in experimental data, same order as RM above
RMEr = [0.01	0.40	0.00	0.00	0.00	0.00	0.00	0.23	0.00
0.01	0.41	0.00	0.00	0.00	0.39	0.00	0.34	0.00
0.01	0.27	0.00	0.00	0.00	0.00	0.00	0.24	0.00
0.02	0.18	0.00	0.00	0.00	0.16	0.00	0.48	0.00
0.01	0.00	0.00	0.00	0.00	0.80	0.00	0.48	0.00
0.01	0.00	0.00	0.00	0.00	0.20	0.00	0.07	0.00
0.01	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.21
];

totno = 100 ; %100 perturbations
rxnno = 95; %reactions in model

%number of experimental conditions
LR = size(RM,1);
%number of experimental rates considered
LC = size(RM,2);

%output
MREout = zeros(LR,rxnno + 1,totno);
 
loop = 1; 

%temporary RM with perturbations
RMtemp = zeros(LR, LC);

%one hundred perturbations loop
while loop <totno + 1  

%%%%%module for perturbing the experimental values within their limits

x = 1; % row counter
y = 1;  % column counter

while x < LR+1   
    
    y = 1 ;
    while y < LC+1 
        
        %RMtemp(x,y) = RM(x,y) + (randi([-1,1],1,1))*RMEr(x,y)*rand  ;
        RMtemp(x,y) = RM(x,y) + RMEr(x,y)*(1 + (-1-1).*rand(1,1));
        y = y +1 ;
         
    end
    
    x = x + 1 ;
    
end

i = 1 ; %experimental condition counter
while i < LR+1 

FBAfunctionMREGLA2800Q10('Models/EcMREnocostv2.mat', 'ATPm',RMtemp(i,1),RMtemp(i,2),RMtemp(i,3), RMtemp(i,5), RMtemp(i,6), RMtemp(i,7), RMtemp(i,8), RMtemp(i,9),MREc(P), MREmu(P), pert);

  if isempty(ans)
    MREout(i,1:rxnno,loop) = 0; %fluxes
    MREout(i,rxnno+1, loop) = 0; %PO number
     
else
    MREout(i,1:rxnno,loop) = ans ;
    MREout(i,rxnno+1,loop) = ans(45)/(ans(29) + ans(39) + ans(42) + ans(43) + ans(45) + ans(46) + ans(81)) + (0.5*ans(47))/(ans(47) + ans(48) + ans(49)) + ans(49)/(ans(47) + ans(48) +  ans(49)); %PO
     
  end

i = i+1 ;
  
end

loop = loop + 1 ;

end

Descr = ['230520EcMREnocost2800_expDataQ10_pert95', '.xlsx'];    %***
file = Descr; 
rang = 'A2:CV97';                       %*** s+2
coln = 'A1:BF1';                        %*** s+2

        %Output1 = Output;
        tab = 'BatchG'; 
        POut1 = MREout ;
        POut1([2:7],:,:) = [] ;
        POp1 = squeeze(POut1) ;
        POp1(:,all(POp1 ==0)) = [];
        xlswrite(file,POp1,tab,rang);
           

        %Output2 = Output;
        tab = 'BatchGL'; 
        POut2 = MREout ;
        POut2([1 3:7],:,:) = [] ;
        POp2 = squeeze(POut2) ;
        POp2(:,all(POp2 ==0)) = [];
        xlswrite(file,POp2,tab,rang);
        
        %Output3 = Output;
        tab = 'BatchGA'; 
        POut3 = MREout ;
        POut3([1 2 4:7],:,:) = [] ;
        POp3 = squeeze(POut3) ;
        POp3(:,all(POp3 ==0)) = [];
        xlswrite(file,POp3,tab,rang);
   
        %Output4 = Output;
        tab = 'BatchGLA'; 
        POut4 = MREout ;
        POut4([1 2 3 5:7],:,:) = [] ;
        POp4 = squeeze(POut4) ;
        POp4(:,all(POp4 ==0)) = [];
        xlswrite(file,POp4,tab,rang);
         

        %Output5 = Output;
        tab = 'BatchL'; 
        POut5 = MREout ;
        POut5([1 2 3 4 6:7],:,:) = [] ;
        POp5 = squeeze(POut5) ;
        POp5(:,all(POp5 ==0)) = [];
        xlswrite(file,POp5,tab,rang);
        
        %Output6 = Output;
        tab = 'BatchLA'; 
        POut6 = MREout ;
        POut6([1:5 7],:,:) = [] ;
        POp6 = squeeze(POut6) ;
        POp6(:,all(POp6 ==0)) = [];
        xlswrite(file,POp6,tab,rang);
        
         %Output7 = Output;
        tab = 'BatchA'; 
        POut7 = MREout ;
        POut7([1:6],:,:) = [] ;
        POp7 = squeeze(POut7) ;
        POp7(:,all(POp7 ==0)) = [];
        xlswrite(file,POp7,tab,rang);
        
       
        
        
        
%clear all
