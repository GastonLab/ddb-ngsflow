__author__ = 'dgaston'

# -*- coding: utf-8 -*-

import sys

from collections import defaultdict

import ConfigParser


def configure_runtime(infile):
    """Parse the configuration settings from a file"""

    configuration = defaultdict
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
    for option in config.options('tools'):
            configuration[option] = config.get('tools', option)

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

            # Check for specific options and if not specified set to defaults
            if 'num_cores' not in tool_dict.keys():
                tool_dict['num_cores'] = config['num_cores']

            if 'max_mem' not in tool_dict.keys():
                tool_dict['max_mem'] = config['max_mem']

            configuration[tool] = tool_dict

    return configuration
