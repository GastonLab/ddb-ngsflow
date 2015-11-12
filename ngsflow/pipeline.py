__author__ = 'dgaston'

# -*- coding: utf-8 -*-

import sys
import subprocess as sub


def run_and_log_command(command, logfile):
    """Run a command and log StdErr to file"""

    with open(logfile, "wb") as err:
        sys.stdout.write("Executing %s and writing to logfile %s\n" % (command, logfile))
        err.write("Command: %s\n" % command)
        p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        output = p.communicate()
        code = p.returncode
        if code:
            raise RuntimeError("An error occurred when executing the commandline: {}. "
                               "Please check the logfile {} for details\n".format(command, logfile))
