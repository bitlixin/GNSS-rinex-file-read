# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:29:37 2019

@author: Xin
"""
import math
import numpy as np
import datetime
import matplotlib.pylab as pl
from PositionSim import *
import copy
Omegae = 7.2921150e-5
Mname = r'C:\Users\Administrator\Desktop\毕业设计\数据\201900203\abmf002d.19n'
#input('请输入卫星电文文件路径和文件名：')
Oname = r'C:\Users\Administrator\Desktop\毕业设计\数据\201900203\abmf002d.19o'
#input('请输入观测文件文件路径和文件名：')
sats = messageRead(Mname)
T,dictOData,Xr = odataRead(Oname)
Xr.append(0)
# 开始解算
XrMean = np.array([0.,0.,0.])
ErDis = []
ErDis2 = []
ErDis3 = []
ErDis4 = []
ErDis5 = []
ErDis6 = []
zheng = 0
fu = 0
for t in T[0:]:
    satNew = copy.deepcopy(sats)
    # 留下离t时刻最近的卫星
    i = 0
    while i < len(satNew) - 1:
        j = i + 1
        while j < len(satNew):
            if satNew[i].prn == satNew[j].prn:
                t1 = abs(t - satNew[i].toe)
                t2 = abs(t - satNew[j].toe)
                if(t1 < t2):
                    del(satNew[j])
                    continue
                else:
                    del(satNew[i])
                    continue
            j += 1
        i += 1

    for i in range(len(satNew)):
        satNew[i].satPosCal(t-0.06667)
        satNew[i].dtdisCal(t)        

    PseDis = []
    k = 0
    Xs = np.zeros([1,3])
    for i in dictOData[t][0]:
        if('G' in i):   #观测文件中的GPS卫星
            #print(i)
            prn = int(i[1:])
            
            for j in range(len(satNew)):
                if(satNew[j].prn == prn):
                    dTdis = satNew[j].dtdis
                    break

            if('P1' in dictOData[t][1][i]):
                PseDis.append(dictOData[t][1][i]['P1'] + dTdis)

#            elif('P2' in dictOData[t][1][i]):
#                PseDis.append(dictOData[t][1][i]['P2'] + dTdis)
                
            elif('C1' in dictOData[t][1][i]):
                PseDis.append(dictOData[t][1][i]['C1'] + dTdis)
                
            else: 
                print("ERROR！未找到测量值！\n")
                
            for j in range(len(satNew)):
                if(satNew[j].prn == prn):
                    #print('ok')
                    satNew[j].satPosCal(t - PseDis[-1]/3e8)
                    X = satNew[j].posArray()
                    [r,lon,lat] = xyz2lonlat(X)
                    lon -= PseDis[-1]/3e8*Omegae
                    X = np.array(lonlat2xyz([r,lon,lat]))
                    h = highAngle(Xr,X)
                    if(h > 5):
                        if(k == 0):
                            Xs[0] += X
                            k += 1
                        else:
                            Xs = np.vstack((Xs,X))
                            k += 1
                    else:
                        del(PseDis[k])

            
    W = np.eye(k)
    for i in range(len(Xs)):
        d = Xs[i,:] - Xr[0:3]
        dis = np.linalg.norm(d)/100000
        #h = highAngle(Xr,Xs[i])
        W[i,i] = 1/dis**2

    PseDis = np.array(PseDis)
    h = np.zeros(PseDis.shape)
    for i in range(len(PseDis)):
        h[i] = highAngle(Xr,Xs[i])
        
        
    
    if(k > 4):
        M = 4
        N = k
        Xr_real = Xr
        X_int = np.array(Xr_real)
        # print(Xs)

        Xr_cal = POS_cal(N,M,Xs,PseDis,X_int)
        Xr_cal2 = POS_cal2(N,M,Xs,PseDis,X_int,W)
        ErDis.append(np.linalg.norm(Xr_cal[0:3] - np.array(Xr[0:3])))

    XsNew,PseNew,dis = satCheck(PseDis,Xr,Xs)
    print('实际坐标',Xr)
    print('计算坐标',Xr_cal)
    print('第一次误差',(Xr - Xr_cal))
    print('误差',np.linalg.norm(Xr - Xr_cal))
    PDOP = PDOPcal(Xs,Xr)
    print('PDOP:',PDOP)
    XsNew,PseNew,dis1 = satCheck(PseDis,Xr_cal,Xs)
    
    W = np.eye(k)
    m = np.median(dis[1,:])
    for i in range(len(Xs)):
        a = dis[i,0] - m
        if(a == 0):
            a = 0.05
        W[i,i] = 1/a**2
    
    #在Xr未知的情况下
    #加权系数与计算出的伪距误差的倒数
    Xr_cal2 = POS_cal2(N,M,Xs,PseDis,X_int,W)
    print('第二次误差',(Xr - Xr_cal2))
    print('误差',np.linalg.norm(Xr - Xr_cal2),'\n')
    ErDis2.append(np.linalg.norm(Xr_cal2[0:3] - np.array(Xr[0:3])))
    
    
    #将误差最大的卫星和伪距删去
    Xr_cal3 = POS_cal(N-2,M,XsNew,PseNew,X_int)
    print('第三次误差',(Xr - Xr_cal3))
    print('误差',np.linalg.norm(Xr - Xr_cal3),'\n')
    ErDis3.append(np.linalg.norm(Xr_cal3[0:3] - np.array(Xr[0:3])))
    
    #以高度角为加权
    height = []
    W = np.eye(k)
    for i in range(len(Xs)):
        h = highAngle(Xr,Xs[i])
        height.append(h)
        W[i,i] = h**2
    
    Xr_cal4 = POS_cal2(N,M,Xs,PseDis,X_int,W)
    print('第四次误差',(Xr - Xr_cal4))
    print('误差',np.linalg.norm(Xr - Xr_cal4),'\n')
    ErDis4.append(np.linalg.norm(Xr_cal4[0:3] - np.array(Xr[0:3])))
    
    #减少误差最大的2项
    W = np.eye(len(PseNew))
    for i in range(len(XsNew)):
        h = highAngle(Xr,XsNew[i])
        W[i,i] = h**2
    
    Xr_cal5 = POS_cal2(N-2,M,XsNew,PseNew,X_int,W)
    print('第5次误差',(Xr - Xr_cal5))
    print('误差',np.linalg.norm(Xr - Xr_cal5),'\n')
    ErDis5.append(np.linalg.norm(Xr_cal5[0:3] - np.array(Xr[0:3])))
    
    #减少高度角最低的一项
    XsNew2,PseNew2,height2 = delHeightMin(PseDis,Xr_cal,Xs)
    W = np.eye(len(PseNew2))
    for i in range(len(XsNew)):
        W[i,i] = height2[i]**2
        
    Xr_cal6 = POS_cal2(N-2,M,XsNew2,PseNew2,X_int,W)
    print('第6次误差',(Xr - Xr_cal6))
    print('误差',np.linalg.norm(Xr - Xr_cal6),'\n')
    ErDis6.append(np.linalg.norm(Xr_cal6[0:3] - np.array(Xr[0:3])))
    
ErDis = np.array(ErDis)
ErDis2 = np.array(ErDis2)
ErDis3 = np.array(ErDis3)
ErDis4 = np.array(ErDis4)
ErDis5 = np.array(ErDis5)
ErDis6 = np.array(ErDis6)
pl.figure(1)
pl.plot(ErDis,color="blue")
pl.plot(ErDis2,color="green")
pl.plot(ErDis3,color="red")
pl.plot(ErDis4,color="yellow")
pl.plot(ErDis5,color="orange")
pl.plot(ErDis6,color="black")
pl.show()