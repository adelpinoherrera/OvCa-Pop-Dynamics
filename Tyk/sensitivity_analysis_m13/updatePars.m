function p = updatePars(p, parsName, parsValue)
for ii = 1:length(parsName)
p.(parsName{ii}) = parsValue(ii);
end
end