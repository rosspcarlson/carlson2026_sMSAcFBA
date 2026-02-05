
%function to run FBA using experimental measurements as constraints

function [FBAflux] = FBAfunction4800Paper1(model_name,objective, mu, gluc_fx, pyr_fx, succ_fx, lactout_fx, lactin_fx, form_fx, aceout_fx, acein_fx, etoh_fx, MREc, MREmu, satnGlc, satnNH3)
  
%enzyme corrections based on experimental proteomics data

ptsGCorr = 1; %omega parameter from Box 1
ATPCorr = 1; %omega parameter from Box 1
ETCCorr = 1; %omega parameter from Box 1

% Open model
model = readCbModel(model_name);
 
% Set model objective
model = changeObjective(model,objective);
 
% Set substrates and constraints 

model = changeRxnBounds(model,'G1',gluc_fx,'b'); % Glucose uptake rate
model = changeRxnBounds(model,'T1',etoh_fx,'b'); % Ethanol secretion rate

model = changeRxnBounds(model,'T2',aceout_fx,'b'); % Acetate secretion rate
model = changeRxnBounds(model,'T2cost',0,'b'); % Acetate secretion rate, cost
model = changeRxnBounds(model,'T2in',acein_fx,'b'); % Acetate uptake rate

model = changeRxnBounds(model,'T4',lactout_fx,'b'); % Lactate secretion  rate
model = changeRxnBounds(model,'T4cost',0,'b'); % Lactate secretion  rate, cost
model = changeRxnBounds(model,'IM3b',0,'b'); % with Lactate uptake rate, needed to avoid LDH cycling
model = changeRxnBounds(model,'IM3',0,'b'); % with Lactate secretion rate, needed to avoid LDH cycling
model = changeRxnBounds(model,'T4in',lactin_fx,'b'); % Lactate uptake rate  

model = changeRxnBounds(model,'T5',form_fx,'b'); % Formate secretion rate
model = changeRxnBounds(model,'T5cost',0,'b'); % Formate secretion rate, cost

model = changeRxnBounds(model,'T11',pyr_fx,'b'); %  pyruvate secretion
model = changeRxnBounds(model,'T11cost',0,'b'); %  pyruvate secretion, cost
model = changeRxnBounds(model,'T11in',0,'b'); %  pyruvate in

model = changeRxnBounds(model,'T12',succ_fx,'b'); % Succinate secretion rate
model = changeRxnBounds(model,'T12cost',0,'b'); % Succinate secretion rate, cost
model = changeRxnBounds(model,'T12in',0,'b'); % Succinate in

model = changeRxnBounds(model,'BIOM',mu,'b'); % Growth rate 
 
model = changeRxnBounds(model,'T25',MREc,'u'); % growth independent MRE constant
model = changeRxnBounds(model,'IM7',0,'b'); % Formate dehydrogenase, anaerobic enzyme
model = changeRxnBounds(model,'T1',0,'b'); % O2 sensitive ethanol dehydrogenase (very little effect on mu or Yx/S)
model = changeRxnBounds(model,'R4',0,'b'); % AppBC enzyme functionally similar to Cyd, higher affinity
model = changeRxnBounds(model,'R6r',0,'b'); % ATP synthase reverse direction, off aerobically to prevent futile cycling
%model = changeRxnBounds(model,'R6f',0,'b'); % ATP synthase forward direction, off anaerobically to prevent futile cycling
%model = changeRxnBounds(model,'T7r',0,'b'); % O2 uptake, anaerobic

% membrane surface area costs, parameters corrected for T, Q10 = 2
model.S(55,1) = (-0.097*ptsGCorr)/(satnGlc) ; %G1 PtsG ***adjusted for satn and experimental data
model.S(55,2) = -3.041 ; %G2a ABC glucose
model.S(55,29) = -0.028 ; %TCA6 Sdh
model.S(55,30) = -0.028 ; %TCA7 Frd
model.S(55,39) = -0.028 ; %IM3b Ldh
model.S(55,42) = -0.116 ; %IM6 Pox
model.S(55,43) = -0.185 ; %IM7 Fdh
model.S(55,45) = -0.378*(ETCCorr) ; %R1 Nuo ***adjusted for experimental data
model.S(55,46) = -0.01*(ETCCorr) ; %R2 Ndh ***adjusted for experimental data
model.S(55,47) = -0.077*(ETCCorr) ; %R3 Cyd ***adjusted for experimental data
model.S(55,49) = -0.049*(ETCCorr) ; %R5 Cyo ***adjusted for experimental data
model.S(55,50) = -0.087*(ATPCorr); %R6f ATP synthase, forward direction
model.S(55,51) = -0.087*(ATPCorr); %R6r ATP synthase, reverse direction

model.S(55,53) = -0.057 ; %R8 Pnt 
model.S(55,64) = -MREmu ; %growth rate dependendent MRE parameter in biomass reaction 
model.S(55,66) = 0 ; %T2 acetate out
model.S(55,67) = -0.065 ; %T2cost acetate
model.S(55,68) = -0.065 ; %T2in acetate

model.S(55,69) = -0.172/satnNH3 ; %T3 NH4+, adjusted for transporter saturation

model.S(55,70) = 0 ; %T4 lactate out
model.S(55,71) = -0.065 ; %T4cost lactate
model.S(55,72) = -0.065 ; %T4in lactate

model.S(55,73) = 0 ; %T5 formate out
model.S(55,74) = -0.065 ; %T5cost formate
 
model.S(55,75) = -0 ; %T6 CO2
model.S(55,76) = -0 ; %T6in CO2

model.S(55,78) = -0.151 ; %T8 SO4
model.S(55,79) = -0 ; %T9 H2O
model.S(55,80) = -0 ; %T9in H2O
model.S(55,81) = -0.151 ; %T10 PO4

model.S(55,82) = 0 ; %T11 pyruvate, no cost  
model.S(55,83) = -0.065 ; %T11cost pyruvate  
model.S(55,84) = -0.065 ; %T11in pyruvate  

model.S(55,85) = 0 ; %T12 succinate, no cost  
model.S(55,86) = -0.065 ; %T12cost succinate  
model.S(55,87) = -0.065 ; %T12in succinate  

% Solve model for optimization criterion
FBAsolution = optimizeCbModel(model,'max');
FBAflux = FBAsolution.x;

end