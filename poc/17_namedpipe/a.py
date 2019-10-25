#!/usr/bin/env python

def Send(msg):
    with open('in', 'w') as f:
        f.write("{:}\n".format(msg))

Send(1)
Send("asdf")
Send("exit")
