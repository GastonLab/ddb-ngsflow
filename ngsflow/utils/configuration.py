"""
.. module:: configuration
   :platform: Unix, OSX
   :synopsis: A module of wrapper methods for configuring pipeline runs.

.. moduleauthor:: Daniel Gaston <daniel.gaston@gmail.com>


"""

__author__ = 'dgaston'

# -*- coding: utf-8 -*-

import sys
import requests

from collections import defaultdict

import ConfigParser


def configure_from_pipeline_service(url, port, pipe_name):
    """Query the DDBio-PipelineService for configuration parameters"""

    payload = {'name': pipe_name}
    r = requests.get("{url}:{port}/pipelines".format(url=url, port=port), params=payload)
    run_parameters = r.json()

    return run_parameters


def configure_runtime(infile):
    """Parse the configuration settings from a file"""

    configuration = defaultdict()
    config = ConfigParser.SafeConfigParser()
    config.read(infile)

    try:
        config.options('settings')
    except ConfigParser.NoSectionError:
        sys.stderr.write("No section: settings in file\n")
        sys.exit()

    try:
        config.options('resources')
    except ConfigParser.NoSectionError:
        sys.stderr.write("No section: resources in file\n")
        sys.exit()

    # Set all options specified in file
    for option in config.options('settings'):
            configuration[option] = config.get('settings', option)

    for resource in config.options('resources'):
            configuration[resource] = config.get('resources', resource)

    # Configure individual tools
    for section in config.sections():
        if section != 'settings' and section != 'resources':
            tool = section
            options = config.options(tool)
            tool_dict = dict()

            # Set all specified options
            for option in options:
                tool_dict[option] = config.get(tool, option)

            configuration[tool] = tool_dict

    return configuration


def configure_samples(infile):
    """Parse the Sample configuration file"""

    samples = dict()

    config = ConfigParser.SafeConfigParser()
    config.read(infile)

    for sample in config.sections():
        sample_dict = dict()
        for option in config.options(sample):
            sample_dict[option] = config.get(sample, option)
        samples[sample] = sample_dict

    return samples
