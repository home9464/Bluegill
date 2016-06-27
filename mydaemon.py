#!/usr/bin/python

"""
crontab -e

#every hour
* */1 * * * python mydaemon.py; [ $? -eq 1 ] && python bluegillframework >> log.txt &
"""

import subprocess,re,sys
cmd_check = 'ps ax | grep bluegillframework'
pattern_check = re.compile(".*?python\s.*bluegillframework")



def shell_exec(cmd,shell=True):
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=shell)
    std_out,std_err = p.communicate()
    if std_out:
        return std_out.strip().split('\n')
    else:
        return []

def checkout():
    ret = shell_exec(cmd_check)
    is_running = False
    for r in ret:
        if pattern_check.match(r):
            is_running = True
    if not is_running:
        return 1
    else:
        return 0

if __name__=='__main__':
    sys.exit(checkout())    
