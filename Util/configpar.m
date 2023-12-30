function parexists = configpar(numpools)
%CONFIGPAR Sets up parallel processing by starting matlabpools. Calling configpar with no parameters launches the maximum available matlabpools. A max of 8 pools is allowed under normal licenses.
%   
%   Usage:
%   parexists = configpar
%   parexists = configpar(numpools)
%
%   Input:
%   numpools: Number of specified pools to start. Max of 12 is allowed.
%   Not specifying this automatically launches max available pools.
%
%   Output:
%   parexists: boolean to see if parallel is installed, same as haspar()
%   function
%
%   Example:
%
%       %Start up 3 matlabpools
%       parexists=configpar(3);
%
%   See also haspar
%
%   Copyright 12/16/2010 Michael J. Prerau, Ph.D. 
%   
%   Last modified 01/07/2011
%********************************************************************

parexists=haspar;

%If the toolbox is installed start pools
if parexists
    out=parcluster;
    max_pools=out.NumWorkers;
    %Start specified number of threads
    if nargin==1
        %Close if different number of threads open
        if matlabpool('size') > 0 && matlabpool('size') ~= numpools
            matlabpool close force;
        end
        
        %Do not let the number of pools exceed the core number
        matlabpool(min(numpools,max_pools));
    elseif matlabpool('size') ~= max_pools
        pool_size=matlabpool('size');
        
        %Check if fewer than max threads or 12 (max allowed by matlab) are open and close
        if pool_size > 0 && pool_size < max_pools
            matlabpool close force;
        end
        
        %Start max threads if not already running
        if pool_size ~= max_pools
            matlabpool('open',max_pools);
        end
    end
else
     Disp('Parallel toolbox not installed');
end
end