
% Function to run FBA using experimental measurements as constraints

function [FBAflux] = FBAfunctionMREGlucose2500Q10(model_name,objective, qatp, mu, MREc, MREmu)
 
% Open model
model = readCbModel(model_name);
 
% Set model objective
model = changeObjective(model,objective);
 
% Set substrate         ***
%model = changeRxnBounds(model,'G1',0,'b'); % glucose
model = changeRxnBounds(model,'T4in',0,'b'); % lactate
model = changeRxnBounds(model,'IM3b',0,'b'); % lactate,(prevents IM3/IM3b LDH loop in place of ndh)
model = changeRxnBounds(model,'T2in',0,'b'); % acetate
model = changeRxnBounds(model,'T11',0,'b'); %  pyruvate
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

model = changeRxnBounds(model,'T25',MREc,'u'); % growth independent MRE constant
model = changeRxnBounds(model,'BIOM',mu,'b'); % specific growth rate (h^-1)
model = changeRxnBounds(model,'ATPm',qatp,'b'); % specific maintenance energy flux

model = changeRxnBounds(model,'IM1',0,'b'); % O2 sensitive Pfl (very little effect on mu or Yx/S)
model = changeRxnBounds(model,'T1',0,'b'); % O2 sensitive ethanol dehydrogenase (very little effect on mu or Yx/S)
model = changeRxnBounds(model,'R4',0,'b'); % AppBC enzyme functionally similar to Cyd, 0.93H+/e, high affinity for O2

%membrane surface area costs Q10 = 2.0
%growth rate dependent PtsG saturation, MG1655 (mu max = 0.67/h, Ks = 5 uM)
if mu < 0.67
    satn = 1.493*mu ;
else
    satn = 1 ;
end

model.S(55,1) = -0.097/satn ; %G1 PtsG 
model.S(55,2) = -3.041 ; %G2a ABC glucose
model.S(55,29) = -0.028 ; %TCA6 Sdh
model.S(55,30) = -0.028 ; %TCA7 Frd
model.S(55,39) = -0.028 ; %IM3b Ldh
model.S(55,42) = -0.116 ; %IM6 Pox
model.S(55,43) = -0.185 ; %IM7 Fdh
model.S(55,45) = -0.378 ; %R1 Nuo
model.S(55,46) = -0.01 ; %R2 Ndh
model.S(55,47) = -0.077 ; %R3 Cyd
model.S(55,49) = -0.049 ; %R5 Cyo
model.S(55,50) = -0.087 ; %R6 ATP synthase
model.S(55,52) = -0.057 ; %R8 Pnt 
model.S(55,63) = -MREmu ; %growth rate dependendent MRE parameter in biomass reaction 
%model.S(55,65) = -0.065 ; %T2 acetate out
model.S(55,66) = -0.065 ; %T2in acetate
model.S(55,67) = -0.172 ; %T3 NH4+
%model.S(55,68) = -0.065 ; %T4 lactate out
model.S(55,69) = -0.065 ; %T4in lactate
%model.S(55,70) = -0.065 ; %T5 formate
model.S(55,71) = -0 ; %T6 CO2
model.S(55,72) = -0 ; %T6in CO2
model.S(55,74) = -0.151 ; %T8 SO4
model.S(55,75) = -0 ; %T9 H2O
model.S(55,76) = -0 ; %T9in H2O
model.S(55,77) = -0.151 ; %T10 PO4
model.S(55,78) = 0 ; %T11 pyruvate, no cost  

% Solve model for optimization criterion
FBAsolution = optimizeCbModel(model,'min'); %min G1 glucose
FBAflux = FBAsolution.x;

end