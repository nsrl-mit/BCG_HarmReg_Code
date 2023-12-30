function gmailConfig(address, passWord)

%%  - - -- --- ----- -------- ------------- -------- ----- --- -- - -
%   
%   configure gmail smtp settings
%   mickey
%   
%   gmailConfig(userName, passWord)
%   
%   Give it your gmail address (i.e. "you.yeahyou@gmail.com") and password
%   and it configures matlab's sendmail function for you. Easy!
%   
%   created:    2014-08-08-1021
%   updated:    2014-08-08-1021
%   
%   - - -- --- ----- -------- ------------- -------- ----- --- -- - -

setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail',address);
setpref('Internet','SMTP_Username',address);
setpref('Internet','SMTP_Password',passWord);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class','javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

