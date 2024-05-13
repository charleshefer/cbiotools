#!/usr/bin/env python
import subprocess
import sys
import time

jobid = sys.argv[1]

retries = 20
i = 0

while i in range(20):
    time.sleep(1)
    #get the status
    output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())
    running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
    if "COMPLETED" in output:
        print("success")
        sys.exit()
    elif any(r in output for r in running_status):
        print("running")
        sys.exit()
    i = i + 1
#tried 20 times, no response, or it failed
print("failed")
