%FindFitness function
% Compare to the expected sfc curve.
% 1 point for having non-negative values
% 1 point for having a ratio greater than 1.1
% 1 point each for each values within 30% range of the expected value
function [output] = exponential_Sfc_M(concentration)

sfc = [0.00760273448820800,0.0101864566941910,0.0136139222327332,0.0181341842428832,0.0240499942911999,0.0317148437179244,0.0415178479478744,0.0538499604451974,0.0690461005234863,0.0873033787501381,0.108587347833183,0.132554059276025,0.158527016470830,0.185561770269678,0.212600363327769,0.238675783939125,0.263101534166405,0.285593520864504,0.306312555044299,0.325858585715879,0.345273694231898,0.366125561978866,0.390776157991274,0.423056762731367,0.469965559726712,0.546503189197296,0.693131717513635,1.07359413237585,4.00601130890706,-1.50471828461896,-0.527969279625462];

score = 0;
tolerance = 0.4; % tolerance range from 0 to 1

if all(concentration>=0)
    score = 1;
   
    ratio = concentration(5)/concentration(1);
    
    if ratio>=1.1
        score = score + 1; % score = 2
    end
    
    y = [concentration(5),concentration(10),concentration(15),concentration(20)]; % adding concentration 25 does not change the score because the last one will always match
    
    if score==2
        for j = 1:4
            if y(j) >= sfc(j*5)*(1-tolerance) && y(j)<= sfc(j*5)*(1+tolerance) %eventually reduce it
                score = score + 1; % best score is 6
            end
        end
    end
    
    output = score;
    
end

