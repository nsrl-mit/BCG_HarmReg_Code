function cOut = findStrCell(sm,cIn)
cOut = cIn(cell2mat(cellfun(@(x)~isempty(strfind(x,sm)),cIn,'un',false)));