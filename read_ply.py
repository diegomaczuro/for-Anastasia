import numpy as np

f = open('surface.ply', 'r')
a = []
vertex_ = False
face_ = False
index = 0

for line in f:
    if line.find('element vertex') != -1:
        print 'evrika'
        vertex_ = True
        line1 = line.split(' ')
        num = int(line1[2])
        print num
        vertex = np.zeros((num, 3))
        index = 0


    if vertex_ and index <= num and index != 0:
        line1 = line.split(' ')
        #n = float(line1[0])
        print line1
        #vertex[index][0] = float(line1[0])
        #vertex[index][1] = float(line1[1])
       # vertex[index][2] = float(line1[2])
    index += 1



f.close()

