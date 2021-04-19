%FindFitness function
% Compare to the expected sfc curve.
% 1 point for having non-negative values
% 1 point for having a ratio greater than 1.1
% 1 point each for each values within 30% range of the expected value

function [output] = exponential_Sfc_M2(concentration)

sfc = [0.0221047858304046,0.0352278449085936,0.0554442306835930,0.0856265102987635,0.128650047912544,0.186128995225380,0.256703382176614,0.335146725042178,0.413861738245364,0.486504776127235,0.551283143017491,0.612319565692932,0.680641940101599,0.779779499174902,0.973241271495922,1.55605548046255,28.9050764890140,-1.05925239429263,-0.395170903303485,-0.195716342758060,-0.107459416457579,-0.0619694996030804,-0.0366855347004706,-0.0220424831085947,-0.0133598512571004,-0.00813946086272128,-0.00497449900841864,-0.00304599569637634,-0.00186729685926280,-0.00114552835080721,-0.000703051693030268];

score = 0;
tolerance = 0.4; % tolerance range from 0 to 1

if all(concentration>=0)
    score = 1;
   
    ratio = concentration(3)/concentration(1);
    
    if ratio>=1.1
        score = score + 1; % score = 2
    end
    
    y = [concentration(3),concentration(6),concentration(9),concentration(12)]; %% why are we only checking the first 12 cells? if that is the case, should not we adjust the range of 15 cells instead?
    
    if score==2
        for j = 1:4
            if y(j) >= sfc(j*3)*(1-tolerance) && y(j)<= sfc(j*3)*(1+tolerance) %eventually reduce it
                score = score + 1; % best score is 6
            end
        end
    end
    
    output = score;
    
end

