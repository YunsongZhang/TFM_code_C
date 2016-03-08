#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


NetBonds=np.loadtxt("./data/RelaxedNet.txt")
numBonds,n=NetBonds.shape
tmp=np.ones((numBonds,2))*float('nan')

PlotData=np.concatenate((NetBonds,tmp),axis=1)
#print PlotData

x1=PlotData[:,0].reshape(numBonds,1)
x2=PlotData[:,2].reshape(numBonds,1)
x3=PlotData[:,4].reshape(numBonds,1)
y1=PlotData[:,1].reshape(numBonds,1)
y2=PlotData[:,3].reshape(numBonds,1)
y3=PlotData[:,5].reshape(numBonds,1)
"""
x1.reshape(numBonds,1)
x2.reshape(numBonds,1)
x3.reshape(numBonds,1)
y1.reshape(numBonds,1)
y2.reshape(numBonds,1)
y3.reshape(numBonds,1)
"""
#print x1
X=np.concatenate((x1,x2,x3),axis=1)
Y=np.concatenate((y1,y2,y3),axis=1)
#print X

Xdata=X.reshape(3*numBonds,1)
Ydata=Y.reshape(3*numBonds,1)
#print Xdata

#fig=plt.figure()
#ax=fig.add_subplot(111)

line2=Line2D(Xdata.T,Ydata.T,color='red',linestyle='solid',label='relaxed')

fig=plt.figure()
ax=fig.add_subplot(111)

ax.add_line(line2)


ax.set_xlim(min(Xdata),max(Xdata))
ax.set_ylim(min(Ydata),max(Ydata))



plt.title("initial net VS relaxed net")
#plt.show()

"""
relaxed net finished, now for predicted net
"""

NetBonds=np.loadtxt("./data/predictedNet.txt")
numBonds,n=NetBonds.shape
tmp=np.ones((numBonds,2))*float('nan')

PlotData=np.concatenate((NetBonds,tmp),axis=1)
#print PlotData

x1=PlotData[:,0].reshape(numBonds,1)
x2=PlotData[:,2].reshape(numBonds,1)
x3=PlotData[:,4].reshape(numBonds,1)
y1=PlotData[:,1].reshape(numBonds,1)
y2=PlotData[:,3].reshape(numBonds,1)
y3=PlotData[:,5].reshape(numBonds,1)
"""
x1.reshape(numBonds,1)
x2.reshape(numBonds,1)
x3.reshape(numBonds,1)
y1.reshape(numBonds,1)
y2.reshape(numBonds,1)
y3.reshape(numBonds,1)
"""
#print x1
X=np.concatenate((x1,x2,x3),axis=1)
Y=np.concatenate((y1,y2,y3),axis=1)
#print X

Xdata=X.reshape(3*numBonds,1)
Ydata=Y.reshape(3*numBonds,1)
#print Xdata

#fig=plt.figure()
#ax=fig.add_subplot(111)

line3=Line2D(Xdata.T,Ydata.T,color='green',linestyle='dotted',label='predicted')

ax.add_line(line3)


ax.set_xlim(min(Xdata),max(Xdata))
ax.set_ylim(min(Ydata),max(Ydata))




plt.legend(handles=[line3],loc=3)


pts1=np.loadtxt("./data/linkedPts.txt")
p1=plt.scatter(pts1[:,0],pts1[:,1],s=20,c='red',label='relaxed')

pts2=np.loadtxt("./data/linkedPtsPredicted.txt")
p2=plt.scatter(pts2[:,0],pts2[:,1],s=10,c='green',label='predicted')


plt.legend(handles=[line2,line3,p1,p2],loc=2,ncol=1,borderaxespad=-0.,prop={'size':6})





#plt.show()
plt.savefig('net_pred')
