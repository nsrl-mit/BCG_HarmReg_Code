function iftttUpdate(message)

tag = '#matlab';
iftttAddress = 'trigger@recipe.ifttt.com';
fullSubject = [ tag ' ' message ];
sendmail(iftttAddress,fullSubject,'');