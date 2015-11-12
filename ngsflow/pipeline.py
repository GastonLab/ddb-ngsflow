__author__ = 'dgaston'

# -*- coding: utf-8 -*-

import sys
import timeit
import subprocess as sub


def run_and_log_command(command, logfile):
    """Run a command and log StdErr to file"""

    start_time = timeit.default_timer()
    with open(logfile, "wb") as err:
        sys.stdout.write("Executing %s and writing to logfile %s\n" % (command, logfile))
        err.write("Command: %s\n" % command)
        p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        output = p.communicate()
        code = p.returncode
        if code:
            sys.stdout.write("An error occurred. Please check %s for details\n" % logfile)
            sys.stdout.write("%s\n" % output)
            sys.stderr.write("An error occurred. Please check %s for details\n" % logfile)

        elapsed = timeit.default_timer() - start_time
        sys.stdout.write("Elapsed Time: %s sec\n" % elapsed)

        return code
