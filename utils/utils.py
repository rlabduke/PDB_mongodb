import sys

def broadcast(m, log = sys.stdout) :
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)
