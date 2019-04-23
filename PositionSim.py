import math
import numpy as np
import datetime

c = 299792458#299708539.6#
u = 3.986004418e14
Omegae = 7.2921150e-5

class Sat:
    """卫星类"""
    def __init__(self, x = 0,y = 0,z = 0):
        self.x = x
        self.y = y
        self.z = z

    def posArray(self):
        a = np.array([self.x,self.y,self.z])
        return a

    def GetPara(self,para):
        self.toe    =   para[0]
        self.sqrA   =   para[1]
        self.e      =   para[2]
        self.w      =   para[3]
        self.M0     =   para[4]
        self.Omega0 =   para[5]
        self.i0     =   para[6]
        self.n1     =   para[7]
        self.Omega1 =   para[8]
        self.IDOT   =   para[9]
        self.Cuc    =   para[10]
        self.Cus    =   para[11]
        self.Crc    =   para[12]
        self.Crs    =   para[13]
        self.Cic    =   para[14]
        self.Cis    =   para[15]
        self.a0     =   para[16]
        self.a1     =   para[17]
        self.a2     =   para[18]
        self.tgd    =   para[19]
        self.prn    =   para[20]


    def satPosCal(self,t):
        """
        计算卫星在t时刻的位置
        t为周积秒
        """
        delt = t - self.toe
        A = self.sqrA**2
        n = (u/A**3)**0.5 + self.n1
        Mk = self.M0 + n * delt
        while(Mk < 0):
            Mk += 2*math.pi
        while(Mk >2*math.pi):
            Mk -= 2*math.pi

        E0 = Mk
        E1 = Mk + self.e*math.sin(E0)
        while(abs(E1 - E0) >0.0000001):
            E0 = E1
            E1 = Mk + self.e*math.sin(E0)

        self.Ek = E1
        vk = math.atan((1-self.e**2)**0.5*math.sin(self.Ek)/(math.cos(self.Ek)-self.e))
        if(math.sin(self.Ek)*vk<0):
            if(vk>0):
                vk -= math.pi
            else:
                vk += math.pi
        phy = vk + self.w
        #计算校正项
        uk1 = self.Cus*math.sin(2*phy) + self.Cuc*math.cos(2*phy)
        rk1 = self.Crs*math.sin(2*phy) + self.Crc*math.cos(2*phy)
        ik1 = self.Cis*math.sin(2*phy) + self.Cic*math.cos(2*phy)
        uk = phy + uk1
        rk = A*(1 - self.e*math.cos(self.Ek)) + rk1
        ik = self.i0 + self.IDOT*delt + ik1
        #计算轨道平面的位置
        xk = rk*math.cos(uk)
        yk = rk*math.sin(uk)
        #计算升交点赤经
        Omegak = self.Omega0 + (self.Omega1 - Omegae)*delt - Omegae*self.toe
        self.x = xk*math.cos(Omegak) - yk*math.cos(ik)*math.sin(Omegak)
        self.y = xk*math.sin(Omegak) + yk*math.cos(ik)*math.cos(Omegak)
        self.z = yk*math.sin(ik)


    def dtdisCal(self,t):
        # 计算卫星钟差
        F = -4.442807633e-10
        deltat = F*self.e*self.sqrA*np.sin(self.Ek)
        self.dt = self.a0 + self.a1*(t - self.toe) + self.a2*(t - self.toe) - self.tgd + deltat
        self.dtdis = c*self.dt
        self.toe2 = self.toe - self.dt

def messageRead(fname):
    #读取导航电文文件
    f = open(fname,'r')
    #读到头文件结尾
    while(True):
        txt = f.readline()
        if 'END OF HEADER' in txt:
            break
    k = 0
    para = []

    #读取数据
    while(True):
        for i in range(8):
            txt = f.readline()
            if(txt == ''):
                break
            if(i == 0):
                l = [txt[0:2].split()[0],txt[2:5].split()[0],txt[5:8].split()[0],txt[8:11].split()[0],txt[11:14].split()[0]
                     ,txt[14:17].split()[0],txt[17:22].split()[0],txt[22:41].split()[0],txt[41:60].split()[0],txt[60:79].split()[0]]
                para.append(l)
            elif(i == 7):
                para [k] += [txt[3:22].split()[0]]
            else:
                l = [txt[3:22].split()[0],txt[22:41].split()[0],txt[41:60].split()[0],txt[60:79].split()[0]]
                para[k] += l
        k += 1
        if(txt == ''):
            break

    len1 = len(para)


    p = []
    sats = []
    for i in range(len1):
        prn     =   int(para[i][0])
        toe     =   float(para[i][18].replace('D','e'))
        sqrA    =   float(para[i][17].replace('D','e'))
        e       =   float(para[i][15].replace('D','e'))
        w       =   float(para[i][24].replace('D','e'))
        M0      =   float(para[i][13].replace('D','e'))
        Omega0  =   float(para[i][20].replace('D','e'))
        i0      =   float(para[i][22].replace('D','e'))
        n1      =   float(para[i][12].replace('D','e'))
        Omega1  =   float(para[i][25].replace('D','e'))
        IDOT    =   float(para[i][26].replace('D','e'))
        Cuc     =   float(para[i][14].replace('D','e'))
        Cus     =   float(para[i][16].replace('D','e'))
        Crc     =   float(para[i][23].replace('D','e'))
        Crs     =   float(para[i][11].replace('D','e'))
        Cic     =   float(para[i][19].replace('D','e'))
        Cis     =   float(para[i][21].replace('D','e'))
        a0      =   float(para[i][7].replace('D','e'))
        a1      =   float(para[i][8].replace('D','e'))
        a2      =   float(para[i][9].replace('D','e'))
        tgd     =   float(para[i][32].replace('D','e'))
        
        p.append([toe,sqrA,e,w,M0,Omega0,i0,n1,Omega1,IDOT,Cuc,Cus,Crc,Crs,Cic,Cis,a0,a1,a2,tgd,prn])
        sats.append(Sat())
        sats[i].GetPara(p[i])
    f.close()
    return sats

def odataRead(fname):
    #读取观测数据文件
    fid = open(fname,'r');
    i = 0
    while(True):
        l = fid.readline()
        if(l == ''):
            break;

        if(l.find('# / TYPES OF OBSERV') != -1):  
            i +=1
            if i ==1:
                list1 = l[:60].split()
                dcnt = int(list1.pop(0))    #观测值类型数量
            else:
                list2 = l[:60].split()
                list1.extend(list2)
        #list1:观测值类型
        if(l.find('APPROX POSITION XYZ') != -1):
            temp = l[:60].split()
            Xr = [float(temp[0]),float(temp[1]),float(temp[2])]

        if(l.find('END OF HEADER') != -1):
            break
    '''
    数据结构：
    字典
    key：时间
    value： [satsprn , data]
    '''
    #观测数据
    dictOData = {}
    T = []
    while(True):
        #时间
        t = fid.readline(29)
        if t == '':
            break
        time = t.split()
        if(time.pop(-1) != '0'):
            print('该历元数据出错！')

        time[0] = '20'+ time[0]
        dtime = []
        dtime = datetime.datetime(int(time[0]),int(time[1]),
                                  int(time[2]),int(time[3]),int(time[4]),
                                  int(float(time[5])))
        dtime1 = datetime.datetime(int(time[0]),int(time[1]),
                                  int(time[2]))
        # print('观测日期：',dtime)
        deltaT = (dtime.weekday()+1)*86400+(dtime-dtime1).seconds   #与历元开始
        if(deltaT > 7*86400):   #周日的时间大于7*86400
            deltaT -= 7*86400
        T.append(deltaT)
            #卫星数量
        t = fid.readline(3)
        cnt = int(t)
            #卫星名
        sat = []
        t = fid.readline(cnt*3)
        t = t[:-1]
        while(len(t)<cnt*3):
            t2 = fid.readline()
            l1 = t2.split()
            t2 = l1[0]
            t = t + t2
            
        for i in range(cnt):
            sat.append(t[i*3:i*3+3])

            #读测量值

        d = {}  #观测值字典
        for i in range(cnt):    #不同卫星
            k = 0
            d1 = {}
            for j in range(dcnt):   #读卫星的不同TYPES OF OBSERV
                t = fid.readline(16)    #读16个字符
                l = t[:16].split()
                if(l != []):
                    d1[list1[k]] = float(l[0])
                    

                if(t[-1] == '\n'):  #若读出的字符串含有换行符
                    if(k == dcnt-1):#最后一个数据
                        break
                    
                    k = k + 5 - k%5 #这时一行为满5个数据
                    if(dcnt<k):     #最后一行
                        break
                    continue
                else:
                    k += 1
                    if(k%5==0): #若读至一行结尾
                        t = fid.readline()
                    elif(k == dcnt):
                        t = fid.readline()
                        break
            d[sat[i]] = d1
        dictOData[deltaT] = [sat,d]

    fid.close()

    return T,dictOData,Xr

def POS_cal(N,M,Xs,PseDis,X_int):
    #计算接收机坐标
    N = len(Xs)
    Xr = np.array(X_int)
    Xs = np.array(Xs)
    G = np.zeros([N,M])
    Bs = np.zeros([N,1])
    for j in range(100):
        for i in range(N):
            R = np.sqrt((Xr[0] - Xs[i,0])**2 + (Xr[1] - Xs[i,1])**2 + (Xr[2] - Xs[i,2])**2)
            G[i,:3] = (Xr[:3]-Xs[i,0:3])/R
            Bs[i] = PseDis[i] - R - Xr[3]
            G[i,3] = 0.5

        GI = np.transpose(G)
        H = np.matmul(GI,G)
        I = np.linalg.inv(H)
        BI = np.matmul(GI,Bs)
        DeltaX = np.matmul(I,BI)
        result = Xr + np.transpose(DeltaX[:,0])
        MagX = np.sqrt(np.vdot(DeltaX[:,0],DeltaX[:,0]))

        if(MagX <1e-6):
            break
        else:
            Xr = result
    return Xr

def POS_cal2(N,M,Xs,PseDis,X_int,W):
    #计算接收机坐标
    Xr = np.array(X_int)
    Xs = np.array(Xs)
    G = np.zeros([N,M])
    Bs = np.zeros([N,1])
    for j in range(100):
        for i in range(N):
            R = np.sqrt((Xr[0] - Xs[i,0])**2 + (Xr[1] - Xs[i,1])**2 + (Xr[2] - Xs[i,2])**2)
            G[i,:3] = (Xr[:3]-Xs[i,0:3])/R
            Bs[i] = PseDis[i] - R - Xr[3]
            G[i,3] = 1

        GI = np.transpose(G)
        H1 = np.matmul(GI,W)
        H = np.matmul(H1,G)
        I = np.linalg.inv(H)
        B1 = np.matmul(W,Bs)
        BI = np.matmul(GI,B1)
        DeltaX = np.matmul(I,BI)
        
        result = Xr + np.transpose(DeltaX[:,0])
        MagX = np.sqrt(np.vdot(DeltaX[:,0],DeltaX[:,0]))

        if(MagX <1e-6):
            break
        else:
            Xr = result
    return Xr


def highAngle(Xr,Xs):
    # 计算卫星高度角
    X1 = np.array(Xr[0:3])
    X2 = np.array(Xs[0:3])
    dx = X2 - X1

    cosangle = np.dot(X1,dx)/(np.linalg.norm(X1)*np.linalg.norm(dx))
    #print(cosangle)
    a = np.pi/2 - math.acos(cosangle)
    a = a/np.pi*180
    return a

def PDOPcal(Xs,Xr):
    N = len(Xs)
    G = np.mat(np.zeros([N,4]))
    for i in range(N):
        R = np.sqrt((Xr[0] - Xs[i,0])**2 + (Xr[1] - Xs[i,1])**2 + (Xr[2] - Xs[i,2])**2)
        G[i,:3] = (Xr[:3]-Xs[i,0:3])/R
        G[i,3] = 1
    H = (G.T*G).I
    PDOP = (H[1,1] + H[2,2] + H[3,3])**0.5
    return PDOP
    

def satCheck(PseDis,Xr_cal,Xs):
    Xr = Xr_cal
    N = len(Xs)
    d = np.mat(np.zeros([N,1]))
    PseDis = np.array(PseDis)
    Xr_new = np.broadcast_to(Xr[0:3],(N,3))
    dd = Xr_new - Xs
    for i in range(N):
        dis = np.linalg.norm(dd[i,:])
        d[i,0] = PseDis[i] - dis
    dis = d
    for i in range(2):
        m = max(abs(d))
        for j in range(len(d)):
            if(abs(d[j,0]) == m):
                d = np.vstack((d[:j,0],d[j+1:,0]))
                Xs = np.vstack((Xs[:j,:],Xs[j+1:,:]))
                PseDis = np.hstack((PseDis[0:j],PseDis[j+1:]))
                break
        
    return Xs,PseDis,dis

def delHeightMin(PseDis,Xr,Xs):
    N = len(Xs)
    height = []
    for i in range(N):
        h = highAngle(Xr,Xs[i])
        height.append(h)
    
    for i in range(2):
        m = min(height)
        for j in range(len(height)):
            if(height[j] == m):
                del(height[j])
                Xs = np.vstack((Xs[:j,:],Xs[j+1:,:]))
                PseDis = np.hstack((PseDis[0:j],PseDis[j+1:]))
                break
    return Xs,PseDis,height


def xyz2lonlat(X):
    # 直角坐标系换算经纬度
    lon = math.atan(X[1]/X[0])
    if(X[0] < 0):
        lon += math.pi
    if(lon < 0):
        lon += 2*math.pi
    
    t = X[2]/math.sqrt(X[0]**2+X[1]**2)
    lat = math.atan(t)
    
    r = np.linalg.norm(X)
    return [r,lon,lat]

def lonlat2xyz(X):
    r,lon,lat = X
    x = r * math.cos(lat)*math.cos(lon)
    y = r * math.cos(lat)*math.sin(lon)
    z = r * math.sin(lat)
    return [x,y,z]

