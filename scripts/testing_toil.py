import argparse
from toil.job import Job


class HelloWorld(Job):
    def __init__(self, message):
        Job.__init__(self,  memory="2G", cores=2, disk="3G")
        self.message = message

    def run(self, fileStore):
        fileStore.logToMaster("Hello, world!, I have a message: %s"
                              % self.message)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    Job.Runner.addToilOptions(parser)
    options = parser.parse_args()
    options.logLevel = "INFO"
    Job.Runner.startToil(HelloWorld("woot"), options)