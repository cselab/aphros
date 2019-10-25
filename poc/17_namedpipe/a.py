#!/usr/bin/env python

import os
import errno


def Send(msg):
    with open('in', 'w') as f:
        f.write("{:}\n".format(msg))

Send(1)
Send("asdf")
Send("exit")
