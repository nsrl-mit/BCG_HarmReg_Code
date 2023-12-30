function askGmailConfig

yn = input('configure gmail settings? (necessary for ifttt progress updates)(y/n):','s');
if ~(strcmp(yn,'y')||strcmp(yn,'Y')||strcmp(yn,'n')||strcmp(yn,'N'))
    clear yn
    askGmailConfig
else
    if strcmp(yn,'y')||strcmp(yn,'Y')
        address = input('address?:','s');
        passwrd = input('password?:','s');
        gmailConfig(address,passwrd);
    end
end