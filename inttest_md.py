#!/usr/bin/env python3

import md
import os, sys


md.run_md()

if(os.path.isfile("cu.traj")):
    sys.exit(0)
else:
    sys.exit(1)
