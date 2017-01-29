"""
  Phasego -- Automatic calculation and plot of phase diagram

  Copyright (C) 2012-2016 by Zhong-Li Liu

  This program is free software; you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software Foundation
  version 3 of the License.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
  PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  E-mail: zl.liu@163.com
"""


# Read in control arguments
def readinput():
    maininfile = open('phasego.in', 'r')
    line = maininfile.readline()
    global indict
    indict = {}
    while line:
        line = maininfile.readline()
        llist = line.split('=')
        if llist != ['\n'] and llist != ['']:
            if llist[0][0] != '#':
                inputlist = [i.strip().split() for i in llist]
                if not inputlist[1]:
                    with open('log.Phasego', 'a') as logfile:
                        print >> logfile, "Please give the value(s) for: %s" % inputlist[0][0]

                    exit(1)
                else:
                    indict[inputlist[0][0]] = inputlist[1]

    return indict

indict = readinput()