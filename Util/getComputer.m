function compName = getComputer()
%%  GETCOMPUTER returns the name of the current computer.
%   
%   compName = getComputer()
%       returns the string compname containing the name of the current
%       computer.
%

if isunix()
    unixGetCNameExp         = 'hostname -fs';
    [~,compName]            = unix(unixGetCNameExp);
else
    dosGetCNameExp          = 'ipconfig /all';
    [~,ipconfig]            = dos(dosGetCNameExp);
    startInd                = 67;    %   change this later
    compName                = strtok(ipconfig(startInd:end));
end
compName                    = deblank(compName);
