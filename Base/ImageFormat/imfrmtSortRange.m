function range = imfrmtSortRange(range)

fn = fieldnames(range);
for i = 1:length(fn)
    f = fn{i};
    range.(f) = sort(range.(f));
end

end