function [output] = modelAll_v3_maptorange(digraphsFinal, parameters)
% Input: digraphs and parameters
% Output: A cell array of calculated fitness scores and idx indication the
% gene responsible for the highest score


%for k = 1 % for testing uncomment

%%%%% THIS IS THE OPTIMISED VERSION TO REDUCE NUMBER %%%%%%
%%%%%%%%% OF CASES MUST BE TRIED %%%%%%%%%%%%%%%%%


% parfor k=1:15 % number of digraphs
%     for i = 1:30 % number of parameters %% van modify this
%         score_M0 = findFitness_Sfc_M_maptorange(modelEuler_v3(zeros(1,96),parameters(i,:),digraphsFinal{k}));
%         if all(score_M0 > 4)
%             score_M1 = findFitness_Sfc_M_maptorange(modelEulerM1_v3(zeros(1,96),parameters(i,:),digraphsFinal{k}));
%             if all(score_M1 >3)
%                 score_M2 = findFitness_Sfc_M2_maptorange(modelEulerM2_v3(zeros(1,96),parameters(i,:),digraphsFinal{k}));
%                 final_score_G1= (score_M0(1)+score_M1(1)+score_M2(1));
%                 final_score_G2= (score_M0(2)+score_M1(2)+score_M2(2));
%                 final_score_G3= (score_M0(3)+score_M1(3)+score_M2(3));
%                 [ffs, index] = max([final_score_G1, final_score_G2, final_score_G3]);
%                 diffEqs{k,i} = [ffs, index];
%             end
%         end   
%     end
% end
% 
% output = diffEqs;
% end


%%%%%%%%%% THIS IS USED TO SEE THE GRAPH FIRST %%%%%%%%%%%%
params = 1;
digraphs = 1;
diffEqs = {};
% M0_G1 = zeros([params*digraphs,25]);
% M0_G2 = zeros([params*digraphs,25]);
% M0_G3 = zeros([params*digraphs,25]);
% M1_G1 = zeros([params*digraphs,25]);
% M1_G2 = zeros([params*digraphs,25]);
% M1_G3 = zeros([params*digraphs,25]);
% M2_G1 = zeros([params*digraphs,15]);
% M2_G2 = zeros([params*digraphs,15]);
% M2_G3 = zeros([params*digraphs,15]);

parfor k = 1:params %% change this back to parfor later!!!!
    for i=1:digraphs
        genes_M0 = modelEuler_v3(zeros(1,96),parameters(i,:),digraphsFinal{k}, 1, 5,5);    
        score_M0 = findFitness_Sfc_M_maptorange(genes_M0);
      
        
        genes_M1 = modelEulerM1_v3(zeros(1,96),parameters(i,:),digraphsFinal{k}, 1, 5,5);
        score_M1 = findFitness_Sfc_M_maptorange(genes_M1);
        
        
        genes_M2 = modelEulerM2_v3(zeros(1,96),parameters(i,:),digraphsFinal{k}, 1, 5,5);
        score_M2 = findFitness_Sfc_M2_maptorange(genes_M2);
        
        final_score_G1= (score_M0(1)+score_M1(1)+score_M2(1));
        final_score_G2= (score_M0(2)+score_M1(2)+score_M2(2));
        final_score_G3= (score_M0(3)+score_M1(3)+score_M2(3));


        [ffs, index] = max([final_score_G1, final_score_G2, final_score_G3]);
        diffEqs{k,i} = [ffs, index];


%         M0{k,i} = genes_M0;
%         M1{k,i} = genes_M1;
%         M2{k,i} = genes_M2;
%         
%         sfc1 = [0.00760273448820800,0.0101864566941910,0.0136139222327332,0.0181341842428832,0.0240499942911999,0.0317148437179244,0.0415178479478744,0.0538499604451974,0.0690461005234863,0.0873033787501381,0.108587347833183,0.132554059276025,0.158527016470830,0.185561770269678,0.212600363327769,0.238675783939125,0.263101534166405,0.285593520864504,0.306312555044299,0.325858585715879,0.345273694231898,0.366125561978866,0.390776157991274,0.423056762731367,0.469965559726712,0.546503189197296,0.693131717513635,1.07359413237585,4.00601130890706,-1.50471828461896,-0.527969279625462];
%         start_index = 1;
%         end_index = 25;
%         min_val = sfc1(start_index);
%         max_val = sfc1(end_index);
%         
%         datag1 = genes_M0(1:3:96);
%         datag2 = genes_M0(2:3:96);
%         datag3 = genes_M0(3:3:96);
%         if (datag1(start_index) ~= datag1(end_index)) % if the minimum and maximum values are the same, using maptorange would give an error
%             dataG1_M0 = maptorange(datag1(start_index:end_index), [datag1(start_index), datag1(end_index)], [min_val, max_val]);
%         else
%             dataG1_M0 = datag1;
%         end
%         if (datag2(start_index) ~= datag2(end_index))
%             dataG2_M0 = maptorange(datag2(start_index:end_index), [datag2(start_index), datag2(end_index)], [min_val, max_val]);
%         else
%             dataG2_M0 = datag2;
%         end
% 
%         if (datag3(start_index) ~= datag3(end_index))
%             dataG3_M0 = maptorange(datag3(start_index:end_index), [datag3(start_index), datag3(end_index)], [min_val, max_val]);
%         else
%             dataG3_M0 = datag3;
%         end 
%         M0_G1(k*digraphs+i,:) = dataG1_M0;
%         M0_G2(k*digraphs+i,:) = dataG2_M0;
%         M0_G3(k*digraphs+i,:) = dataG3_M0;
%         
%         
%         datag1 = genes_M1(1:3:96);
%         datag2 = genes_M1(2:3:96);
%         datag3 = genes_M1(3:3:96);
%         if (datag1(start_index) ~= datag1(end_index)) % if the minimum and maximum values are the same, using maptorange would give an error
%             dataG1_M1 = maptorange(datag1(start_index:end_index), [datag1(start_index), datag1(end_index)], [min_val, max_val]);
%         else
%             dataG1_M1 = datag1;
%         end
%         if (datag2(start_index) ~= datag2(end_index))
%             dataG2_M1 = maptorange(datag2(start_index:end_index), [datag2(start_index), datag2(end_index)], [min_val, max_val]);
%         else
%             dataG2_M1 = datag2;
%         end
% 
%         if (datag3(start_index) ~= datag3(end_index))
%             dataG3_M1 = maptorange(datag3(start_index:end_index), [datag3(start_index), datag3(end_index)], [min_val, max_val]);
%         else
%             dataG3_M1 = datag3;
%         end
%         M1_G1(k*digraphs+i,:) = dataG1_M1;
%         M1_G2(k*digraphs+i,:) = dataG2_M1;
%         M1_G3(k*digraphs+i,:) = dataG3_M1;
%         
%         sfc2 = [0.0221047858304046,0.0352278449085936,0.0554442306835930,0.0856265102987635,0.128650047912544,0.186128995225380,0.256703382176614,0.335146725042178,0.413861738245364,0.486504776127235,0.551283143017491,0.612319565692932,0.680641940101599,0.779779499174902,0.973241271495922,1.55605548046255,28.9050764890140,-1.05925239429263,-0.395170903303485,-0.195716342758060,-0.107459416457579,-0.0619694996030804,-0.0366855347004706,-0.0220424831085947,-0.0133598512571004,-0.00813946086272128,-0.00497449900841864,-0.00304599569637634,-0.00186729685926280,-0.00114552835080721,-0.000703051693030268];
%         start_index = 1;
%         end_index = 15;
%         min_val = sfc2(start_index);
%         max_val = sfc2(end_index);
%         datag1 = (1:3:96);
%         datag2 = genes_M2(2:3:96);
%         datag3 = genes_M2(3:3:96);
% 
%         if (datag1(start_index) ~= datag1(end_index)) % if the minimum and maximum values are the same, using maptorange would give an error
%             dataG1_M2 = maptorange(datag1(start_index:end_index), [datag1(start_index), datag1(end_index)], [min_val, max_val]);
%         else
%             dataG1_M2 = datag1;
%         end
%         if (datag2(start_index) ~= datag2(end_index))
%             dataG2_M2 = maptorange(datag2(start_index:end_index), [datag2(start_index), datag2(end_index)], [min_val, max_val]);
%         else
%             dataG2_M2 = datag2;
%         end
% 
%         if (datag3(start_index) ~= datag3(end_index))
%             dataG3_M2 = maptorange(datag3(start_index:end_index), [datag3(start_index), datag3(end_index)], [min_val, max_val]);
%         else
%             dataG3_M2 = datag3;
%         end
%         M2_G1(k*digraphs+i,:) = dataG1_M2;
%         M2_G2(k*digraphs+i,:) = dataG2_M2;
%         M2_G3(k*digraphs+i,:) = dataG3_M2;
        
        
        
    end
       
        
end



output = diffEqs;
end
    
        
        
        


