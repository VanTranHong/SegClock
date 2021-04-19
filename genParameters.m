function [outputArg1] = genParameters()

params_num = 1000;
parameters = zeros(params_num,12);
for k = 1:params_num % change this number according to how many parameters you want to generate
    V = zeros(1,12);
    for j=1:9
     i = randi(10000);
     R = 10*rand();% van has changed this
       V(j) = (0.9995.^(i))*R;    
    end
    for a = 1:3
       R = 0.05*rand();% van has changed this, instead of 10
       V(9+a) = (0.9996.^10000)*R;
    end
    parameters(k,:) = V;
end

outputArg1 = parameters;
end

