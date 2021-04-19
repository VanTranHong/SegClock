function [scores, maximum] = maxscore(job)

output = job{1,1};
[nr,nc] = size(output);
for i = 1:nr
    for j=1:nc
        scores(i,j) = output{i,j}(1);
    end
end

maximum = max(max(scores));

end
        