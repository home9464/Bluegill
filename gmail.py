"""user gmail as SMTP server to send email notifications to users
"""
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email import Encoders
import os
import imp

class Gmail:
    def __init__(self,to,job_name):
        """to= ('user1@bar.com','user2@bar.com')
        """
        self.util = imp.load_source('util','util.py')
        self.to = to
        self.job_name = job_name
            
    def send(self,
             subject,
             body_msg=None,
             body_file=None,
             additional_subject=None,
             attach=None):
        #subject = "Job %s finished" % self.job_name
        #else:
        #    subject = "Job %s failed" % self.job_name

        if additional_subject:
            subject = "%s, %s" % (subject,additional_subject)
        
        body=[]
        if body_msg:
            body.append(body_msg)
            body.append('\n')
            body.append('='*50)
            body.append('\n')
            
        if body_file:
            try:
                fh = open(body_file)
                body.append(fh.read())
                fh.close()
            
            except:
                pass
        
        for recipient in self.to:
            try: 
                msg = MIMEMultipart()
                msg['From'] = "Bluegill"        
                msg['To'] = recipient        
                msg['Subject'] = subject
                #msg.attach(MIMEText(body))
                msg.attach(MIMEText('\n'.join(body)))
                
                if attach:
                    part = MIMEBase('application', 'octet-stream')
                    part.set_payload(open(attach,'rb').read())
                    Encoders.encode_base64(part)
                    part.add_header('Content-Disposition','attachment; filename="%s"' % os.path.basename(attach))
                    msg.attach(part)
                
                #mailServer = smtplib.SMTP("smtp.gmail.com", 587)
                mailServer = smtplib.SMTP(self.util.SMTP_SERVER)
                #mailServer.ehlo()
                #mailServer.starttls()
                #mailServer.ehlo()
                #mailServer.login(GMAIL_USER, GMAIL_PWD)
                #mailServer.sendmail(GMAIL_USER, recipient, msg.as_string())
                mailServer.sendmail("Bluegill", recipient, msg.as_string())
                # Should be mailServer.quit(), but that crashes...
                mailServer.close()
            except:
                pass