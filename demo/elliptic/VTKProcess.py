#! /usr/bin/env python
import sys
import string
from array import array
VecScal = string.atoi(sys.argv[1]) # 0 save no vectors/scalars, 1 save scalars, 2 save vectors
fpIn=open(sys.argv[2], 'r')
fpOut=open(sys.argv[3], 'w')
NumOfVertices=-1
while 1: # through whole file
    line = fpIn.readline()
    if not line: # already at end of file
        break
    if line[0] != '#' or line[1] == ' ': # this is not a line to edit
        fpOut.write(line)
        continue
    if line[1] == '1':
        parts = line.split(' ')
        for i in range (1,6,2): # gives 1,3,5
            if parts[i] == 'DimOfProb':
                DimOfProb = string.atoi(parts[i+1])
                approx = array('f')
                exact = array('f')
                for j in range (DimOfProb):
                    approx.append(0)
                    exact.append(0)
            elif parts[i] == 'DimOfWorld':
                DimOfWorld = string.atoi(parts[i+1])            
            elif parts[i] == 'ExactSol':
                ExactSol = string.atoi(parts[i+1])
    elif line[1] == '2':
        if VecScal == 1: 
            if line[3] == 'L':
                fpOut.write(line[3:])
            else:
                fpOut.write('SCALARS diff float 1\n')
    elif line[1] == '3':
        if VecScal == 2: 
            fpOut.write('VECTORS diff float\n')
    elif line[1] == '4':
        if VecScal == 0: # I want to save nothing
            continue
        parts = line[3:].split()
        for i in range(DimOfProb):
            approx[i] = string.atof(parts[i])
            if ExactSol:
                exact[i] = string.atof(parts[i+DimOfProb])
        if VecScal == 1: # writing scalar values
            fpOut.write('%(r1)f\n' % {'r1': approx[0]-exact[0]})
        elif VecScal == 2: # writing vector values
            fpOut.write('%(r1)f %(r2)f %(r3)f\n' % {'r1': approx[0]-exact[0], 'r2': approx[1]-exact[1], 'r3': approx[2]-exact[2]})
fpOut.write("\n")
fpIn.seek(0)
while 1: # through whole file
    line = fpIn.readline()
    if not line: # already at end of file
        break
    if line[0] != '#' or line[1] == ' ': # this is not a line to edit
        continue
    if line[1] == '1':
        continue
    elif line[1] == '2':
        if VecScal == 1: 
            if line[3] == 'L':
                fpOut.write(line[3:])
            else:
                fpOut.write('SCALARS approx float 1\n')
    elif line[1] == '3':
        if VecScal == 2: 
            fpOut.write('VECTORS approx float\n')
    elif line[1] == '4':
        if VecScal == 0: # I want to save nothing
            continue
        parts = line[3:].split()
        for i in range(DimOfProb):
            approx[i] = string.atof(parts[i])
            if ExactSol:
                exact[i] = string.atof(parts[i+DimOfProb])
        if VecScal == 1: # writing scalar values
            fpOut.write('%(r1)f\n' % {'r1': approx[0]})
        elif VecScal == 2: # writing vector values
            fpOut.write('%(r1)f %(r2)f %(r3)f\n' % {'r1': approx[0], 'r2': approx[1], 'r3': approx[2]})
fpOut.write("\n")
fpIn.seek(0)
while 1: # through whole file
    line = fpIn.readline()
    if not line: # already at end of file
        break
    if line[0] != '#' or line[1] == ' ': # this is not a line to edit
        continue
    if line[1] == '1':
        continue
    elif line[1] == '2':
        if VecScal == 1: 
            if line[3] == 'L':
                fpOut.write(line[3:])
            else:
                fpOut.write('SCALARS exact float 1\n')
    elif line[1] == '3':
        if VecScal == 2: 
            fpOut.write('VECTORS exact float\n')
    elif line[1] == '4':
        if VecScal == 0: # I want to save nothing
            continue
        parts = line[3:].split()
        for i in range(DimOfProb):
            approx[i] = string.atof(parts[i])
            if ExactSol:
                exact[i] = string.atof(parts[i+DimOfProb])
        if VecScal == 1: # writing scalar values
            fpOut.write('%(r1)f\n' % {'r1': exact[0]})
        elif VecScal == 2: # writing vector values
            fpOut.write('%(r1)f %(r2)f %(r3)f\n' % {'r1': exact[0], 'r2': exact[1], 'r3': exact[2]})
