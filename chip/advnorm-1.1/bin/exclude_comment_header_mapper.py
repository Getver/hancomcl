#!/usr/bin/env python

############################
## Copyright (c) 2017 Thermo Fisher Scientific.  All Rights Reserved. This source
## code and the methods embodied within it represent the confidential
## intellectual property of Thermo Fisher Scientific, and this intellectual
## property is communicated on a confidential basis to the recipient.  Neither
## this code nor the methods may not be disclosed to any other party in any
## form without the express written permission of Thermo Fisher Scientific.  This
## copyright and notice must be preserved in any derivative work.
############################

import sys

for line_s in sys.stdin:

   if (line_s[0] == '#'): continue
   if (line_s[0:11] == "probeset_id"): continue
   sys.stdout.write(line_s)

exit(0)
