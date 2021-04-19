function [output] = findFitness_Sfc_M_maptorange(concentrations)
%Find the fitness scores
sfc1 = [0.00760273448820800,0.0101864566941910,0.0136139222327332,0.0181341842428832,0.0240499942911999,0.0317148437179244,0.0415178479478744,0.0538499604451974,0.0690461005234863,0.0873033787501381,0.108587347833183,0.132554059276025,0.158527016470830,0.185561770269678,0.212600363327769,0.238675783939125,0.263101534166405,0.285593520864504,0.306312555044299,0.325858585715879,0.345273694231898,0.366125561978866,0.390776157991274,0.423056762731367,0.469965559726712,0.546503189197296,0.693131717513635,1.07359413237585,4.00601130890706,-1.50471828461896,-0.527969279625462];
start_index = 1;
end_index = 25;
min_val = sfc1(start_index);
max_val = sfc1(end_index);
fitness = [];

datag1 = concentrations(1:3:96);
datag2 = concentrations(2:3:96);
datag3 = concentrations(3:3:96);

%%%%%%%%%%%%%%%%%%%%% van has modified it so that it is scaled according to
%%%%%%%%%%%%%%%%%%%%% sfc of the first 25 cells

if (datag1(start_index) ~= datag1(end_index)) % if the minimum and maximum values are the same, using maptorange would give an error
    dataG1 = maptorange(datag1(start_index:end_index), [datag1(start_index), datag1(end_index)], [min_val, max_val]);
else
    dataG1 = datag1;
end
if (datag2(start_index) ~= datag2(end_index))
    dataG2 = maptorange(datag2(start_index:end_index), [datag2(start_index), datag2(end_index)], [min_val, max_val]);
else
    dataG2 = datag2;
end

if (datag3(start_index) ~= datag3(end_index))
    dataG3 = maptorange(datag3(start_index:end_index), [datag3(start_index), datag3(end_index)], [min_val, max_val]);
else
    dataG3 = datag3;
end

fitness(1) = exponential_Sfc_M(dataG1);
fitness(2) = exponential_Sfc_M(dataG2);
fitness(3) = exponential_Sfc_M(dataG3);

output = fitness;
end

