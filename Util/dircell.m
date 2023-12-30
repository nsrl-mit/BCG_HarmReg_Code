function dCell = dircell(directory)
lsDirectory = ls(directory);
dotPad = 2;
dCell = cell(size(lsDirectory,1)-dotPad,1);
kk = 1;
% keyboard
for k = (dotPad+1):size(lsDirectory,1)
    dirContent = deblank(lsDirectory(k,:));         %   grab directory listing
    if strncmp(dirContent(end-3:end),'.mat',4);     %   check if it's a .mat file
        dCell{kk} = dirContent;
        kk = kk+1;
    end
end