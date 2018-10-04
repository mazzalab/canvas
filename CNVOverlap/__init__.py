#!/usr/bin/env python

"""
TITLE
"""
import os
import sys
import getopt
import re

__author__ = "Mauro Truglio"
__copyright__ = "Copyright 201X"
__credits__ = [""]
__version__ = "0.0.1"
__maintainer__ = "Mauro Truglio"
__email__ = "bioinformatics@css-mendel.it"
__status__ = "Development"
__date__ = ""
__license__ = u"""
  Copyright (C) 2016-2017  Mauro Truglio <m.truglio@css-mendel.it>
  Viale Regina Margherita 261, 00198 Rome, Italy

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301 USA
  """


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """
    return [atoi(c) for c in re.split('(\d+)', text)]


def print_help():
    print("Help")


if __name__ == '__main__':

    if len(sys.argv) == 1:
        print_help()

    opts, args = getopt.getopt(sys.argv[1:], "hc:")

    for opt, arg in opts:
        if opt == '-h':
            print_help()
        if opt in "-c":
            foo = arg

