function frmt = imfrmtFlip(refFrmt, frmt)
% flips al formats in frmt to have the case of the format in refFrmt if present


[id, pos] = ismember(lower(refFrmt), lower(frmt));
pos =pos(id);
frmt(pos) = refFrmt(id);

end