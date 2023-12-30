function userString = getUser
%%  getUser returns a string identifying the current user.
%   
%   userString = getUser
%   
%   
%   \////\ Michael Nolan 092911022015

%%  return 
if isunix()
    [~,userString] = unix('echo $USER');
else
    [~,userString] = dos('echo %USERNAME%');
end
%   remove whitespace from tail of string
userString = deblank(userString);
