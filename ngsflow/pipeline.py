"""
.. module:: pipeline
   :platform: Unix, OSX
   :synopsis: A wrapper module for executing commands and logging their output.

.. moduleauthor:: Daniel Gaston <daniel.gaston@dal.ca>


"""

import sys
import subprocess as sub


def run_and_log_command(command, logfile):
    """This function uses the python subprocess method to run the specified command and writes all error to the
    specified logfile

    :param command: The command-line command to execute.
    :type name: str.
    :param logfile: The logfile to output error messages to.
    :type logfile: str.
    :returns:  Nothing
    :raises: RuntimeError

    """

    with open(logfile, "wb") as err:
        sys.stdout.write("Executing %s and writing to logfile %s\n" % (command, logfile))
        err.write("Command: %s\n" % command)
        p = sub.Popen(command, stdout=sub.PIPE, stderr=err, shell=True)
        output = p.communicate()
        code = p.returncode
        if code:
            raise RuntimeError("An error occurred when executing the commandline: {}. "
                               "Please check the logfile {} for details\n".format(command, logfile))
