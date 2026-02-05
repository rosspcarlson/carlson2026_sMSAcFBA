
%FBAfunctionMREGLA2800
%December 2022 - April 2023
%Ross P. Carlson
%Function to run FBA using experimental measurements as constraints,
%modified to perturb the enzyme turnover numbers +/-'pert'% to test robustness of
% solutions



function [FBAflux] = FBAfunctionMREGLA2800Q10(model_name,objective, mu, gluc_fx, pyr_fx, lactout_fx, lactin_fx, form_fx, aceout_fx, acein_fx, MREc, MREmu, pert)
  
  % Open model
model = readCbModel(model_name);
 
% Set model objective
model = changeObjective(model,objective);
 
% Set substrates and constraints 

model = changeRxnBounds(model,'G1',gluc_fx,'b'); % Glucose uptake rate
model = changeRxnBounds(model,'T4',lactout_fx,'b'); % Lactate secretion  rate
model = changeRxnBounds(model,'IM3b',lactin_fx,'b'); % Lactate uptake rate, needed to avoid LDH cycling
model = changeRxnBounds(model,'T4in',lactin_fx,'b'); % Lactate uptake rate  
model = changeRxnBounds(model,'T2',aceout_fx,'b'); % Acetate secretion rate
model = changeRxnBounds(model,'T2in',acein_fx,'b'); % Acetate uptake rate
model = changeRxnBounds(model,'T5',form_fx,'b'); % Formate secretion rate
model = changeRxnBounds(model,'BIOM',mu,'b'); % Growth rate 
model = changeRxnBounds(model,'T11',-pyr_fx,'b'); %  pyruvate

model = changeRxnBounds(model,'T25',MREc,'u'); % growth independent MRE constant
model = changeRxnBounds(model,'IM7',0,'b'); % Formate dehydrogenase, anaerobic enzyme
model = changeRxnBounds(model,'T1',0,'b'); % O2 sensitive ethanol dehydrogenase (very little effect on mu or Yx/S)
model = changeRxnBounds(model,'R4',0,'b'); % AppBC enzyme functionally similar to Cyd, 0.93H+/e, unknown cellular role

model = changeRxnBounds(model,'T12',0,'b'); % mannose
model = changeRxnBounds(model,'T13',0,'b'); % xylose
model = changeRxnBounds(model,'T14',0,'b'); % glycerol
model = changeRxnBounds(model,'T15',0,'b'); % succinate
model = changeRxnBounds(model,'T16',0,'b'); % oxaloacetate
model = changeRxnBounds(model,'T17',0,'b'); % maltose
model = changeRxnBounds(model,'T18',0,'b'); % 2-oxoglutarate
model = changeRxnBounds(model,'T19',0,'b'); % phosphoenolpyruvate
model = changeRxnBounds(model,'T20',0,'b'); % fumarate
model = changeRxnBounds(model,'T21',0,'b'); % ribose
model = changeRxnBounds(model,'T22',0,'b'); % malate
model = changeRxnBounds(model,'T23',0,'b'); % mannitol
model = changeRxnBounds(model,'T24',0,'b'); % sorbitol

%temperature adjusted TO numbers
%membrane surface area costs 
model.S(55,1) = -0.097 + 0.097*pert*(1 + (-1-1).*rand(1,1));  %G1 PtsG 
model.S(55,2) = -3.041 + 3.041*pert*(1 + (-1-1).*rand(1,1));  %G2a ABC glucose
model.S(55,29) = -0.028 + 0.028*pert*(1 + (-1-1).*rand(1,1));  %TCA6 Sdh
model.S(55,30) = -0.028 + 0.028*pert*(1 + (-1-1).*rand(1,1));  %TCA7 Frd
model.S(55,39) = -0.028 + 0.028*pert*(1 + (-1-1).*rand(1,1));  %IM3b Ldh
model.S(55,42) = -0.116 + 0.116*pert*(1 + (-1-1).*rand(1,1));  %IM6 Pox
model.S(55,43) = -0.185 + 0.185*pert*(1 + (-1-1).*rand(1,1));  %IM7 Fdh
model.S(55,45) = -0.378 + 0.378*pert*(1 + (-1-1).*rand(1,1));  %R1 Nuo
model.S(55,46) = -0.010 + 0.010*pert*(1 + (-1-1).*rand(1,1));  %R2 Ndh
model.S(55,47) = -0.077 + 0.077*pert*(1 + (-1-1).*rand(1,1));  %R3 Cyd
model.S(55,49) = -0.049 + 0.049*pert*(1 + (-1-1).*rand(1,1));  %R5 Cyo
model.S(55,50) = -0.087 + 0.087*pert*(1 + (-1-1).*rand(1,1));  %R6 ATP synthase
model.S(55,52) = -0.057 + 0.057*pert*(1 + (-1-1).*rand(1,1));  %R8 Pnt 
model.S(55,63) = -MREmu ; %growth rate dependendent MRE parameter in biomass reaction 
%model.S(55,65) = -0.065 ; %T2 acetate out
model.S(55,66) = -0.065 + 0.065*pert*(1 + (-1-1).*rand(1,1));  %T2in acetate
model.S(55,67) = -0.172 + 0.172*pert*(1 + (-1-1).*rand(1,1));  %T3 NH4+
%model.S(55,68) = -0.065 ; %T4 lactate out
model.S(55,69) = -0.065 + 0.065*pert*(1 + (-1-1).*rand(1,1));  %T4in lactate
%model.S(55,70) = -0.065 ; %T5 formate out
model.S(55,71) = -0 ; %T6 CO2
model.S(55,72) = -0 ; %T6in CO2
model.S(55,74) = -0.151 + 0.151*pert*(1 + (-1-1).*rand(1,1)) ; %T8 SO4
model.S(55,75) = -0 ; %T9 H2O
model.S(55,76) = -0 ; %T9in H2O
model.S(55,77) = -0.151 + 0.151*pert*(1 + (-1-1).*rand(1,1)) ;  %T10 PO4
model.S(55,78) = 0 ; %pyruvate, no cost  

% Solve model for optimization criterion
FBAsolution = optimizeCbModel(model,'max');
FBAflux = FBAsolution.x;

end