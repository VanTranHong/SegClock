function [output] = map_M2_torange(concentrations)
%Find the fitness scores for each gene


sfc2 = [0.0221047858304046,0.0352278449085936,0.0554442306835930,0.0856265102987635,0.128650047912544,0.186128995225380,0.256703382176614,0.335146725042178,0.413861738245364,0.486504776127235,0.551283143017491,0.612319565692932,0.680641940101599,0.779779499174902,0.973241271495922,1.55605548046255,28.9050764890140,-1.05925239429263,-0.395170903303485,-0.195716342758060,-0.107459416457579,-0.0619694996030804,-0.0366855347004706,-0.0220424831085947,-0.0133598512571004,-0.00813946086272128,-0.00497449900841864,-0.00304599569637634,-0.00186729685926280,-0.00114552835080721,-0.000703051693030268];


%%%%%%%%%%%%%%%%%%%%% van has modified it so that it is scaled according to
%%%%%%%%%%%%%%%%%%%%% sfc of the first 25 cells
start_index = 1;
end_index = 15;
min_val = sfc2(start_index);
max_val = sfc2(end_index);
datag1 = concentrations(1:3:96);
datag2 = concentrations(2:3:96);
datag3 = concentrations(3:3:96);

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
output = [dataG1;dataG2;dataG3];
end