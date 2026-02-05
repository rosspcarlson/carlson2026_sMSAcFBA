%OLMREpaper1, Over Lord Membrane Real Estate paper 1
%June/July 2024
%top management code: cycles through experimental data series, calls FBA function with
%constraints, and logs output
% Ross P. Carlson

%aerobic experimental data indexing
%rows: 1-23, cstr/batch, glucose limited, consensus data MG1655
%rows: 24-46, engineered PtsG/batch, glucose limited, consensus data NCM3722

i = 24 ; %experimental condition counter, above***
imax = 46 ; %max i for experimental conditions***
zmax = 20 ; %fMSA range considered 
rowtrack = 97 ; % number of model reactions plus added functions

outmF = zeros(rowtrack, 1);  
outsF = zeros(rowtrack, 1);


while i < imax+1
    
[mF,sF]=functionEc4800Paper1(i, zmax, rowtrack); %used for both cost and no cost
 
outmF = [outmF mF]; %mean of 100 fluxes
outsF = [outsF sF]; %standard deviation of 100 fluxes

i = i+ 1;

end

outmF(:,1) = [];
outsF(:,1) = [];


Descr = ["250602Ec4800NCM71sum.xlsx"]; 

file = Descr; 
                       
        tab = 'outmF'; 
        
        xlswrite(file,outmF,tab);
              
        tab = 'outsF'; 
        
        xlswrite(file,outsF,tab);