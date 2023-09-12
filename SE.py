import math
import OutTxt
import numpy as np
import scipy.linalg
from numpy.linalg import inv
from OutTxt import Real
from math import pi,sqrt
from numpy import arctan,sin,cos,tan,arccos
def ADSE(U,Angle,Vd,Id,Kt,W,fi,Udc,Idc,derta,M,Ydc,Yac,NodeData,LineData,LCC_NodeData,VSC_NodeData,DC_LineData,AC_LCC,AC_VSC,Tol,**Option):
#------------------------------------------LCC不平衡量----------------------------------------#
    Nc = LCC_NodeData.shape[0]
    kr = 0.995
    Vt = np.zeros([Nc,1])                      # 交流电压
    for i in range(Nc):
        Vt[i] = U[AC_LCC[i]]
    Uc=np.append(Vd,Udc)                       # LCC 与 MMC 直流电压
    N_LCC = LCC_NodeData[:,12]                 # 换流器组数
    N_VSC = VSC_NodeData[:,12]                 
    N = np.append(N_LCC,N_VSC)
    X = LCC_NodeData[:,5]                  #换相电抗
    Xc= LCC_NodeData[:,6]                  #无功补偿电纳 注意是电纳
    Pss = np.zeros([Nc,1])
    Qss = np.zeros([Nc,1])
    Pd = np.zeros([Nc,1])
    # 控制参数——伪量测
    Kontrol1 = LCC_NodeData[:,3] 
    Kontrol2 = LCC_NodeData[:,4]               
    #量测值
    Psm = LCC_NodeData[:,13]
    Qsm = LCC_NodeData[:,14]
    Pdm = LCC_NodeData[:,15]
    VdM = LCC_NodeData[:,16]
    IdM = LCC_NodeData[:,17] 
    # 10个量测方程
    Delta_D1 = np.zeros([Nc,1])   # 交流侧有功量测             
    Delta_D2 = np.zeros([Nc,1])   # 交流侧无功量测 
    Delta_D3 = np.zeros([Nc,1])   # 直流有功量测
    Delta_D4 = np.zeros([Nc,1])   # 直流电压量测
    Delta_D5 = np.zeros([Nc,1])   # 直流电流量测
    Delta_D6 = np.zeros([Nc,1])   # 直流电压伪量测            
    Delta_D7 = np.zeros([Nc,1])   # 直流电压伪量测
    Delta_D8 = np.zeros([Nc,1])   # 直流网络方程伪量测
    Delta_D9 = np.zeros([Nc,1])   # 控制方程伪量测
    Delta_D10 = np.zeros([Nc,1])  # 控制方程伪量测

    for i in range(Nc):                        # 计算Delta_
        Pss[i] = N_LCC[i]*Vd[i]*Id[i]
        Qss[i] = N_LCC[i]*Vd[i]*Id[i]*np.tan(fi[i])
        # s-Vt[i]*Vt[i]*Xc #注意
        Pd[i] = N_LCC[i]*Vd[i]*Id[i]
        Delta_D1[i] = Psm[i]-Pss[i]
        Delta_D2[i] = Qsm[i]-Qss[i]
        Delta_D3[i] = Pdm[i]-Pd[i]
        Delta_D4[i] = VdM[i]-Vd[i]
        Delta_D5[i] = IdM[i]-Id[i]
        Delta_D6[i] = Vd[i]-2.7*Kt[i]*Vt[i]*W[i]+1.9*X[i]*Id[i] # 伪量测
        Delta_D7[i] = Vd[i]-2.7*kr*Kt[i]*Vt[i]*cos(fi[i])
        if LCC_NodeData[i,2]==1:  # 定电流 定控制角
            Delta_D8[i] = Id[i]-np.sum(Ydc[i,:]*Uc*N)
            Delta_D9[i] =Id[i]-Kontrol1[i]
            Delta_D10[i] = W[i]-cos(Kontrol2[i])       
        else:                     # 定电压 定控制角
            Delta_D8[i] = -Id[i]-np.sum(Ydc[i,:]*Uc*N)
            Delta_D9[i] =Vd[i]-Kontrol1[i]
            Delta_D10[i] = W[i]-cos(Kontrol2[i])            
#------------------------------------------VSC不平衡量----------------------------------------#
    VSC_Num = VSC_NodeData.shape[0]
    R = VSC_NodeData[:,6]
    Xl = VSC_NodeData[:,7]
    a = arctan(R/Xl)
    Y = 1/np.sqrt(R*R+Xl*Xl)
    Usi = np.zeros([VSC_Num,1])
    for i in range(VSC_Num):
        Usi[i] = U[AC_VSC[i]]                 # 交流母线电压 
    
    Pv = np.zeros([VSC_Num,1])
    Qv = np.zeros([VSC_Num,1])
    Pdc = np.zeros([VSC_Num,1])
    PM = VSC_NodeData[:,13]#量测值
    QM = VSC_NodeData[:,14]
    PdcM = VSC_NodeData[:,15]
    UM = VSC_NodeData[:,16]
    IM = VSC_NodeData[:,17]
    #伪量测 
    Pcontronl = VSC_NodeData[:,3]                     # 换流器控制值(迭代过程中不变)
    Ucontronl = VSC_NodeData[:,8] 
    Qcontronl = VSC_NodeData[:,4]
    #-----8量测-----#
    Deltad1 = np.zeros([VSC_Num,1])   # 交流侧有功量测
    Deltad2 = np.zeros([VSC_Num,1])   # 交流侧无功量测
    Deltad3 = np.zeros([VSC_Num,1])   # 直流功率量测
    Deltad4 = np.zeros([VSC_Num,1])   # 直流网络方程伪量测
    Deltad5 = np.zeros([VSC_Num,1])   # 直流电压量测
    Deltad6 = np.zeros([VSC_Num,1])   # 直流电流量测
    Deltad7 = np.zeros([VSC_Num,1])   # 控制方程伪量测
    Deltad8 = np.zeros([VSC_Num,1])   # 控制方程伪量测
    iter = 0 
    for i in range(VSC_Num):                   # 求解功率不平衡量
        Pv[i] =  (sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i]) + N_VSC[i]*Usi[i]*Usi[i]*Y[i]*sin(a[i])
        Qv[i] = -(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i]) + N_VSC[i]*Usi[i]*Usi[i]*Y[i]*cos(a[i])                     
        Pdc[i] = N_VSC[i]*Udc[i]*Idc[i]
        Deltad1[iter] = PM[i]-Pv[i] 
        Deltad2[iter] = QM[i]-Qv[i] 
        Deltad3[iter] = PdcM[i]-Pdc[i] 
        Deltad4[iter] = Idc[i]-np.sum((Ydc[i+Nc,:]*Uc*N))
        Deltad5[iter] = UM[i]- Udc[i] 
        Deltad6[iter] = IM[i]- Idc[i]
         #--------1-PQ,2-UQ-----------#
        if  VSC_NodeData[i,2]==1:
            Deltad7[iter] = Pv[i]-Pcontronl[i] 
            Deltad8[iter] = Qv[i]-Qcontronl[i]
        else:
            Deltad7[iter] = Udc[i]-Ucontronl[i] 
            Deltad8[iter] = Qv[i]-Qcontronl[i]
        iter = iter+1
    #  直流线路功率量测
    n_DC=Ydc.shape[0]
    n_Line=DC_LineData.shape[0]
    xdc= [int(i) for i in DC_LineData[:,0]]
    ydc= [int(i) for i in DC_LineData[:,1]]
    Pijdc_Real= DC_LineData[:,3]
    Pijdc = np.zeros([n_Line,1])
    DeltaPijdc=np.zeros([n_Line,1])
    for i in range(n_Line):
        q = xdc[i]-1
        r = ydc[i]-1
        Pijdc[i]=-N[q]*Uc[q]*Ydc[q,r]*(N[q]*Uc[q]-N[r]*Uc[r])
        DeltaPijdc[i] = Pijdc_Real[i]-Pijdc[i] 
    # print(Pijdc)    
#-----------------------------------------------交流不平衡量--------------------------------------#                                                               # 节点数目                                                   
    NumNode = Yac.shape[0]
    nlines =  LineData.shape[0]
    G = Yac.real
    B = Yac.imag
    # 节点电压、注入功率及线路功率量测
    U_Real = NodeData[:,7]
    P_Real = NodeData[:,5]                                                 
    Q_Real = NodeData[:,6]
    Pij_Real= LineData[:,6]
    Qij_Real= LineData[:,7]                                               # 节点输入无功功率量测
    # 节点注入功率计算值
    P = np.zeros([NumNode,1])
    Q = np.zeros([NumNode,1])
    Pij= np.zeros([nlines,1])
    Qij= np.zeros([nlines,1])
    DeltaU = np.zeros([NumNode,1])
    DeltaP = np.zeros([NumNode,1])
    DeltaQ = np.zeros([NumNode,1])
    DeltaPij = np.zeros([nlines,1])
    DeltaQij = np.zeros([nlines,1])
    
    # 求解不平衡量
    for i in range(NumNode):                 
        if i in AC_LCC:
            jz = np.where(LCC_NodeData[:,0]==(i+1))
            if LCC_NodeData[jz,2]==1:
                P[i] = U[i]*np.sum(U*(G[i,:]*np.cos(Angle[i]-Angle) +  B[i,:]*np.sin(Angle[i]-Angle))) + Pss[jz]
                Q[i] = U[i]*np.sum(U*(G[i,:]*np.sin(Angle[i]-Angle) -  B[i,:]*np.cos(Angle[i]-Angle))) + Qss[jz]
            else:
                P[i] = U[i]*np.sum(U*(G[i,:]*np.cos(Angle[i]-Angle) +  B[i,:]*np.sin(Angle[i]-Angle))) - Pss[jz]
                Q[i] = U[i]*np.sum(U*(G[i,:]*np.sin(Angle[i]-Angle) -  B[i,:]*np.cos(Angle[i]-Angle))) + Qss[jz]   
        elif i in AC_VSC:
            jzz = np.where(VSC_NodeData[:,0]==(i+1))
            P[i] = U[i]*np.sum(U*(G[i,:]*np.cos(Angle[i]-Angle) +  B[i,:]*np.sin(Angle[i]-Angle))) + Pv[jzz]
            Q[i] = U[i]*np.sum(U*(G[i,:]*np.sin(Angle[i]-Angle) -  B[i,:]*np.cos(Angle[i]-Angle))) + Qv[jzz] 
        else:
            P[i] = U[i]*np.sum(U*(G[i,:]*np.cos(Angle[i]-Angle) +  B[i,:]*np.sin(Angle[i]-Angle)))
            Q[i] = U[i]*np.sum(U*(G[i,:]*np.sin(Angle[i]-Angle) -  B[i,:]*np.cos(Angle[i]-Angle)))
        
        DeltaU[i] = U_Real[i]-U[i]
        DeltaP[i] = P_Real[i]-P[i]     
        DeltaQ[i] = Q_Real[i]-Q[i] 
    
    #  交流线路功率量测 
    x = [int(i) for i in LineData[:,0]]
    start_bus = x
    y = [int(i) for i in LineData[:,1]] 
    end_bus = y
    for i in range(nlines):
        q = start_bus[i]-1
        r = end_bus[i]-1
        Pij[i]=-((U[q])**2)*G[q,r]+U[q]*U[r]*(np.cos(Angle[q]-Angle[r])*G[q,r]+np.sin(Angle[q]-Angle[r])*B[q,r])
        Qij[i]=((U[q])**2)*B[q,r]+U[q]*U[r]*(np.sin(Angle[q]-Angle[r])*G[q,r]-np.cos(Angle[q]-Angle[r])*B[q,r])
        DeltaPij[i] = Pij_Real[i]-Pij[i]     
        DeltaQij[i] = Qij_Real[i]-Qij[i] 
    
    DeltaD = np.vstack([DeltaU,DeltaP,DeltaQ,DeltaPij,DeltaQij,Delta_D1,Delta_D2,Delta_D3,Delta_D4,Delta_D5,Delta_D6,Delta_D7,Delta_D8,Delta_D9,Delta_D10,Deltad1,Deltad2,Deltad3,Deltad4,Deltad5,Deltad6,Deltad7,Deltad8,DeltaPijdc])  # 不平衡量
    # Option['string'] = '不平衡量为：\n'
    # Real(DeltaD,**Option)
    # MaxError = np.max(np.abs(DeltaD))
#----------------------------------------------雅克比矩阵-------------------------------------------#
#-----------------------------------------------交流-----------------------------------------------#
    # 状态变量：功角Angle 电压U
    # 电压量测方程偏导
    Hv = np.zeros([NumNode,NumNode]) 
    Hd = -np.diag(U)  #偏导数乘以U 修正量中分母有U
    # 注入功率量测方程偏导
    H = np.zeros([NumNode,NumNode])
    Ns = np.zeros([NumNode,NumNode]) 
    J = np.zeros([NumNode,NumNode])
    L = np.zeros([NumNode,NumNode])
    # 线路传输功率率
    Hpij1 = np.zeros((nlines,NumNode))
    Hpij2 = np.zeros((nlines,NumNode))
    Hqij1 = np.zeros((nlines,NumNode))
    Hqij2 = np.zeros((nlines,NumNode))
    H_iter = -1                         # H代表行
    for i in range(NumNode):
        N_iter = -1                     # N代表列               
        H_iter = H_iter+1               
        for j in range(NumNode):
                N_iter = N_iter+1   
                if i != j:             
                    Angleij = Angle[i]-Angle[j]
                    H[H_iter,N_iter] = -U[i]*U[j]*(G[i,j]*np.sin(Angleij)-B[i,j]*np.cos(Angleij))
                    J[H_iter,N_iter] = U[i]*U[j]*(G[i,j]*np.cos(Angleij)+B[i,j]*np.sin(Angleij))
                    Ns[H_iter,N_iter] = -U[i]*U[j]*(G[i,j]*np.cos(Angleij)+B[i,j]*np.sin(Angleij))
                    L[H_iter,N_iter] = -U[i]*U[j]*(G[i,j]*np.sin(Angleij)-B[i,j]*np.cos(Angleij))
                else:
                    H[H_iter,N_iter] = U[i]*np.sum(U*(G[i,:]*np.sin(Angle[i]-Angle) -  B[i,:]*np.cos(Angle[i]-Angle)))+U[i]**2*B[i,i]
                    J[H_iter,N_iter] = -U[i]*np.sum(U*(G[i,:]*np.cos(Angle[i]-Angle) +  B[i,:]*np.sin(Angle[i]-Angle)))+G[i,i]*U[i]**2
                    Ns[H_iter,N_iter] = -U[i]*np.sum(U*(G[i,:]*np.cos(Angle[i]-Angle) +  B[i,:]*np.sin(Angle[i]-Angle)))-G[i,i]*U[i]**2
                    L[H_iter,N_iter] = -U[i]*np.sum(U*(G[i,:]*np.sin(Angle[i]-Angle) -  B[i,:]*np.cos(Angle[i]-Angle)))+B[i,i]*U[i]**2

    line_tracker = np.array([x,y])    
    for i in range(nlines):
        temp = line_tracker[:,i]
        temp[0] -= 1
        temp[1] -= 1 #此处对line_tracker的所有行列坐标进行了减一操作，后面无需在进行此操作
        for j in range(NumNode):
            if j == temp[0]:
                Hpij1[i,j] = -U[temp[0]]*U[temp[1]]*(-G[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]])+B[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]]))
                Hqij1[i,j] = -U[temp[0]]*U[temp[1]]*(G[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]])+B[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]]))
            if j == temp[1]:
                Hpij1[i,j] = -U[temp[1]]*U[temp[1]]*(G[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]])-B[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]]))
                Hqij1[i,j] = -U[temp[1]]*U[temp[1]]*(-G[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]])-B[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]]))        
            if j == temp[0]:
                Hpij2[i,j] = -U[temp[0]]*(-2*U[temp[0]]*G[temp[0],temp[1]]+U[temp[1]]*(G[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]])+B[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]])))
                Hqij2[i,j] = -U[temp[0]]*(2*U[temp[0]]*B[temp[0],temp[1]]+U[temp[1]]*(G[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]])-B[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]])))
            if j == temp[1]:
                Hpij2[i,j] = -U[temp[1]]*(U[temp[0]]*(G[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]])+B[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]])))
                Hqij2[i,j] = -U[temp[1]]*(U[temp[0]]*(G[temp[0],temp[1]]*np.sin(Angle[temp[0]]-Angle[temp[1]])-B[temp[0],temp[1]]*np.cos(Angle[temp[0]]-Angle[temp[1]])))
    Jaccobi = np.vstack([np.hstack([Hv,Hd]),np.hstack([H,Ns]),np.hstack([J,L]),np.hstack([Hpij1,Hpij2]),np.hstack([Hqij1,Hqij2])])
#  修正接MMC的交流雅克比
    for i in range(VSC_Num):                
            Jaccobi[NumNode+AC_VSC[i],NumNode+AC_VSC[i]] = Jaccobi[NumNode+AC_VSC[i],NumNode+AC_VSC[i]]-(sqrt(6)/4)*N_VSC[i]*M[i]*Udc[i]*Usi[i]*Y[i]*sin(derta[i]-a[i])-2*N_VSC[i]*Usi[i]*Usi[i]*Y[i]*sin(a[i])                     # P对U
            Jaccobi[2*NumNode+AC_VSC[i],NumNode+AC_VSC[i]] = Jaccobi[2*NumNode+AC_VSC[i],NumNode+AC_VSC[i]]+(sqrt(6)/4)*N_VSC[i]*M[i]*Udc[i]*Usi[i]*Y[i]*cos(derta[i]-a[i])-2*N_VSC[i]*Usi[i]*Usi[i]*Y[i]*cos(a[i]) # Q对U
#------------------------------------------------LCC------------------------------------------------#   
    # 状态变量：Vd Id Kt W fi
    #Psm
    F11 = -np.diag(N_LCC*Id)
    F12 = -np.diag(N_LCC*Vd)
    F13 = np.zeros([Nc,Nc])
    F14 = np.zeros([Nc,Nc])
    F15 = np.zeros([Nc,Nc])
    #Qsm
    F21 = -np.diag(N_LCC*Id*tan(fi.reshape(Nc)))
    F22 = -np.diag(N_LCC*Vd*tan(fi.reshape(Nc)))
    F23 = np.zeros([Nc,Nc])
    F24 = np.zeros([Nc,Nc])
    F25 = -np.diag(N_LCC*Vd*Id*(1/cos(fi.reshape(Nc)))*(1/cos(fi.reshape(Nc))))
    #Pdm
    F31 = -np.diag(N_LCC*Id)
    F32 = -np.diag(N_LCC*Vd)
    F33 = np.zeros([Nc,Nc])
    F34 = np.zeros([Nc,Nc])
    F35 = np.zeros([Nc,Nc])
    #Vdm
    F41 = -np.eye(Nc)
    F42 = np.zeros([Nc,Nc])
    F43 = np.zeros([Nc,Nc])
    F44 = np.zeros([Nc,Nc])
    F45 = np.zeros([Nc,Nc])
    #Idm
    F51 = np.zeros([Nc,Nc])
    F52 = -np.eye(Nc)
    F53 = np.zeros([Nc,Nc])
    F54 = np.zeros([Nc,Nc])
    F55 = np.zeros([Nc,Nc])
    #电压wei量测
    F61 = np.eye(Nc)
    F62 = np.diag(1.9*X)
    F63 = -np.diag(Vt.reshape(Nc)*W*2.7)
    F64 = -np.diag(Kt*Vt.reshape(Nc)*2.7)
    F65 = np.zeros([Nc,Nc])
    F71 = np.eye(Nc)
    F72 = np.zeros([Nc,Nc])
    F73 = -np.diag(kr*Vt.reshape(Nc)*cos(fi.reshape(Nc))*2.7)
    F74 = np.zeros([Nc,Nc])
    F75 = np.diag(kr*Kt*Vt.reshape(Nc)*sin(fi.reshape(Nc))*2.7)
    #网络方程wei量测
    F81 = np.zeros([Nc,Nc])
    F82 = np.eye(Nc)
    for i in range(Nc):
        if LCC_NodeData[i,2]==2:
            F82[i,i]=F82[i,i]*(-1)
    F83 = np.zeros([Nc,Nc])
    F84 = np.zeros([Nc,Nc])
    F85 = np.zeros([Nc,Nc])
    #控制伪量测
    F91 = np.zeros([Nc,Nc])
    F92 = np.zeros([Nc,Nc])
    for i in range(Nc):
        if LCC_NodeData[i,2]==2:
            F91[i,i]=1  #定电压
        else:             
            F92[i,i]=1  #定电流
    F93 = np.zeros([Nc,Nc])
    F94 = np.zeros([Nc,Nc])
    F95 = np.zeros([Nc,Nc])
    F101 = np.zeros([Nc,Nc])
    F102 = np.zeros([Nc,Nc])
    F103 = np.zeros([Nc,Nc])
    F104 = np.eye(Nc)  #控制角
    F105 = np.zeros([Nc,Nc])
    F = np.vstack([np.hstack([F11,F12,F13,F14,F15]),np.hstack([F21,F22,F23,F24,F25]),np.hstack([F31,F32,F33,F34,F35]),np.hstack([F41,F42,F43,F44,F45]),np.hstack([F51,F52,F53,F54,F55]),np.hstack([F61,F62,F63,F64,F65]),np.hstack([F71,F72,F73,F74,F75]),np.hstack([F81,F82,F83,F84,F85]),np.hstack([F91,F92,F93,F94,F95]),np.hstack([F101,F102,F103,F104,F105])]) 
#----------------------------------------------VSC-------------------------------------------------#
    # 状态变量：Udc Idc derta M
    # Psm
    D11=np.zeros([VSC_Num,VSC_Num])
    D12=np.zeros([VSC_Num,VSC_Num]) 
    D13=np.zeros([VSC_Num,VSC_Num])
    D14=np.zeros([VSC_Num,VSC_Num])
    # Qsm
    D21=np.zeros([VSC_Num,VSC_Num])
    D22=np.zeros([VSC_Num,VSC_Num]) 
    D23=np.zeros([VSC_Num,VSC_Num])
    D24=np.zeros([VSC_Num,VSC_Num])
    # Pdm
    D31=np.zeros([VSC_Num,VSC_Num])
    D32=np.zeros([VSC_Num,VSC_Num]) 
    D33=np.zeros([VSC_Num,VSC_Num])
    D34=np.zeros([VSC_Num,VSC_Num])
    # 直流网络方程伪量测
    D41=np.zeros([VSC_Num,VSC_Num])
    D42=np.eye(VSC_Num)             
    D43=np.zeros([VSC_Num,VSC_Num]) 
    D44=np.zeros([VSC_Num,VSC_Num])
    # Vdm
    D51=-np.eye(VSC_Num)   # 单位阵
    D52=np.zeros([VSC_Num,VSC_Num]) # 0
    D53=np.zeros([VSC_Num,VSC_Num]) # 0阵 
    D54=np.zeros([VSC_Num,VSC_Num]) # 0
    # Idm
    D61=np.zeros([VSC_Num,VSC_Num])
    D62=-np.eye(VSC_Num)   # 单位阵 
    D63=np.zeros([VSC_Num,VSC_Num]) # 0阵 
    D64=np.zeros([VSC_Num,VSC_Num]) # 0
    #------P or U------#
    D71=np.zeros([VSC_Num,VSC_Num])
    D72=np.zeros([VSC_Num,VSC_Num])
    D73=np.zeros([VSC_Num,VSC_Num]) # 0阵 
    D74=np.zeros([VSC_Num,VSC_Num]) # 0
    #------Q------#
    D81=np.zeros([VSC_Num,VSC_Num])
    D82=np.zeros([VSC_Num,VSC_Num]) # 0
    D83=np.zeros([VSC_Num,VSC_Num]) # 0阵 
    D84=np.zeros([VSC_Num,VSC_Num]) # 0 
    for i in range(VSC_Num):
        D11[i,i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Y[i]*sin(derta[i]-a[i])
        D13[i,i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i])
        D14[i,i]=-(sqrt(6)/4)*N_VSC[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i])
        D21[i,i]=(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Y[i]*cos(derta[i]-a[i])
        D23[i,i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i])
        D24[i,i]=(sqrt(6)/4)*N_VSC[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i])
        D31[i,i]=-N_VSC[i]*Idc[i]
        D32[i,i]=-N_VSC[i]*Udc[i]
        if VSC_NodeData[i,2]==2:
            D71[i,i]=1
    D = np.vstack([np.hstack([D11,D12,D13,D14]),np.hstack([D21,D22,D23,D24]),np.hstack([D31,D32,D33,D34]),np.hstack([D41,D42,D43,D44]),np.hstack([D51,D52,D53,D54]),np.hstack([D61,D62,D63,D64]),np.hstack([D71,D72,D73,D74]),np.hstack([D81,D82,D83,D84])])      
#-----------------------------------------------整合雅克比矩阵---------------------------------------#
    J_D = scipy.linalg.block_diag(F,D)
    #   补充直流线路功率偏导
    DCline_tracker = np.array([xdc,ydc])    
    Pdcij=np.zeros([n_Line,5*Nc+4*VSC_Num])
    for i in range(n_Line):
        temp = DCline_tracker[:,i]
        temp[0] -= 1
        temp[1] -= 1 
        for j in range(Nc):
            if j == temp[0]:    
                Pdcij[i,j] = -(-2*N[temp[0]]*N[temp[0]]*Uc[temp[0]]+N[temp[0]]*N[temp[1]]*Uc[temp[1]])*Ydc[temp[0],temp[1]]  #   VSC
            if j == temp[1]:
                Pdcij[i,j] = -N[temp[0]]*N[temp[1]]*Uc[temp[0]]*Ydc[temp[0],temp[1]]
        for j in range(VSC_Num):
            if j == (temp[0]-Nc):    
                Pdcij[i,5*Nc+j] = -(-2*N[temp[0]]*N[temp[0]]*Uc[temp[0]]+N[temp[0]]*N[temp[1]]*Uc[temp[1]])*Ydc[temp[0],temp[1]]  #   VSC
            if j == (temp[1]-Nc):
                Pdcij[i,5*Nc+j] = -N[temp[0]]*N[temp[1]]*Uc[temp[0]]*Ydc[temp[0],temp[1]]
    J_DD = np.vstack([J_D,Pdcij])
    J_J = scipy.linalg.block_diag(Jaccobi,J_DD)
#-----------------------------------------------交流对LCC非对角--------------------------------------#
    for i in range(Nc):
        if LCC_NodeData[i,2]==1:
            J_J[NumNode+AC_LCC[i],2*NumNode+i]=-N_LCC[i]*Id[i]                                              # P对Udc偏导,-1是减掉平衡节点
            J_J[NumNode+AC_LCC[i],2*NumNode+Nc+i]=-N_LCC[i]*Vd[i]                                           # P对Idc偏导
            J_J[2*NumNode+AC_LCC[i],2*NumNode+i]=-N_LCC[i]*Id[i]*tan(fi[i])                                 # Q对Udc偏导
            J_J[2*NumNode+AC_LCC[i],2*NumNode+Nc+i]=-N_LCC[i]*Vd[i]*tan(fi[i])                              # Q对Idc偏导
            J_J[2*NumNode+AC_LCC[i],2*NumNode+4*Nc+i]=-N_LCC[i]*Vd[i]*Id[i]*(1/cos(fi[i]))*(1/cos(fi[i]))   
        else:                                            
            J_J[NumNode+AC_LCC[i],2*NumNode+i]=N_LCC[i]*Id[i]                                                      
            J_J[NumNode+AC_LCC[i],2*NumNode+Nc+i]=N_LCC[i]*Vd[i]   
            J_J[2*NumNode+AC_LCC[i],2*NumNode+i]=-N_LCC[i]*Id[i]*tan(fi[i])                               
            J_J[2*NumNode+AC_LCC[i],2*NumNode+Nc+i]=-N_LCC[i]*Vd[i]*tan(fi[i])                              
            J_J[2*NumNode+AC_LCC[i],2*NumNode+4*Nc+i]=-N_LCC[i]*Vd[i]*Id[i]*(1/cos(fi[i]))*(1/cos(fi[i]))           
        #-------------------------------------------------------------------------------------------#
        J_J[3*NumNode+2*nlines+5*Nc+i,NumNode+AC_LCC[i]]=-2.7*Kt[i]*Vt[i]*W[i]                                    
        J_J[3*NumNode+2*nlines+6*Nc+i,NumNode+AC_LCC[i]]=-2.7*kr*Kt[i]*Vt[i]*cos(fi[i])                        
#-----------------------------------------------交流对MMC非对角--------------------------------------#
    for i in range(VSC_Num):
        J_J[NumNode+AC_VSC[i],5*Nc+2*NumNode+i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Y[i]*sin(derta[i]-a[i])                                    # P对Udc偏导
        J_J[NumNode+AC_VSC[i],5*Nc+2*NumNode+2*VSC_Num+i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i])                   # P对a偏导
        J_J[NumNode+AC_VSC[i],5*Nc+2*NumNode+3*VSC_Num+i]=-(sqrt(6)/4)*N_VSC[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i])                        # P对M偏导
        J_J[2*NumNode+AC_VSC[i],5*Nc+2*NumNode+i]=(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Y[i]*cos(derta[i]-a[i])                                   # Q对Udc偏导
        J_J[2*NumNode+AC_VSC[i],5*Nc+2*NumNode+2*VSC_Num+i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i])                 # Q对a偏导
        J_J[2*NumNode+AC_VSC[i],5*Nc+2*NumNode+3*VSC_Num+i]=(sqrt(6)/4)*N_VSC[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i])                       # Q对M偏导
        #-------------------------------------------------------------------------------------------#
        J_J[3*NumNode+2*nlines+10*Nc+i,NumNode+AC_VSC[i]]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Udc[i]*Usi[i]*Y[i]*sin(derta[i]-a[i])-2*N_VSC[i]*Usi[i]*Usi[i]*Y[i]*sin(a[i])          # P对Us偏导
        J_J[3*NumNode+2*nlines+10*Nc+VSC_Num+i,NumNode+AC_VSC[i]]=(sqrt(6)/4)*N_VSC[i]*M[i]*Udc[i]*Usi[i]*Y[i]*cos(derta[i]-a[i])-2*N_VSC[i]*Usi[i]*Usi[i]*Y[i]*cos(a[i])   # Q对Us偏导
#---------------------------------------------LCC与MMC直流网络方程偏导-------------------------------# 
    for i in range(Nc):
        for j in range(Nc):
            J_J[3*NumNode+2*nlines+7*Nc+i,2*NumNode+j]=-N_LCC[j]*Ydc[i,j]                                    #注意：LCC标号靠前
        for j in range(VSC_Num):
            J_J[3*NumNode+2*nlines+7*Nc+i,2*NumNode+5*Nc+j]=-N_VSC[j]*Ydc[i,Nc+j]
    for i in range(VSC_Num):
        for j in range(Nc):
            J_J[3*NumNode+2*nlines+10*Nc+3*VSC_Num+i,2*NumNode+j]=-N_LCC[j]*Ydc[Nc+i,j]
        for j in range(VSC_Num):
            J_J[3*NumNode+2*nlines+10*Nc+3*VSC_Num+i,2*NumNode+5*Nc+j]=-N_VSC[j]*Ydc[Nc+i,Nc+j]
#-------------------------------------------------求解----------------------------------------------#
    # Option['string'] = 'jacobi矩阵为：\n'
    # Real(J_J,**Option)
    Gain = J_J.T.dot(J_J)
    # Gain
    Gain[0,0] = 10000000
    niag = inv(Gain)
    Delta = inv(Gain).dot(J_J.T).dot(DeltaD)
    # print(DeltaD)
    # Option['string'] = '方程组求解结果：\n'
    # Real(Delta,**Option)
#-------------------------------------------------修正------------------------------------------#
    Angle = Angle- Delta[0:NumNode].reshape(NumNode)                                                           
    for i in range(NumNode):
        U[i] = U[i]-U[i]*Delta[NumNode+i]
    Vd = Vd-Delta[2*NumNode:2*NumNode+Nc].reshape(Nc)                                   
    Id = Id-Delta[2*NumNode+Nc:2*NumNode+2*Nc].reshape(Nc)
    Kt = Kt-Delta[2*NumNode+2*Nc:2*NumNode+3*Nc].reshape(Nc)
    W= W-Delta[2*NumNode+3*Nc:2*NumNode+4*Nc].reshape(Nc)
    fi = fi-Delta[2*NumNode+4*Nc:2*NumNode+5*Nc]
    Udc = Udc-Delta[2*NumNode+5*Nc:2*NumNode+5*Nc+VSC_Num].reshape(VSC_Num)         
    Idc = Idc-Delta[2*NumNode+5*Nc+VSC_Num:2*NumNode+5*Nc+2*VSC_Num].reshape(VSC_Num)
    derta = derta-Delta[2*NumNode+5*Nc+2*VSC_Num:2*NumNode+5*Nc+3*VSC_Num].reshape(VSC_Num)
    M = M-Delta[2*NumNode+5*Nc+3*VSC_Num:2*NumNode+5*Nc+4*VSC_Num].reshape(VSC_Num)

    Option['string'] = '更新之后的电压幅值为：\n'
    Real(U,**Option)
    Option['string'] = '相角为：\n'
    Real(Angle,**Option)
    Option['string'] = '\nLCC更新之后的直流电压为：\n'
    Real(Vd,**Option)
    Option['string'] = '直流电流为：\n'
    Real(Id,**Option)
    Option['string'] = '换流变变比：\n'
    Real(Kt,**Option)
    Option['string'] = '控制角\n：'
    Real(W,**Option)
    # Real(57.3*arccos(W),**Option)
    Option['string'] = '功率因数：\n'
    Real(cos(fi),**Option)
    Option['string'] = '\nVSC更新之后的直流电压为：\n'
    Real(Udc,**Option)
    Option['string'] = '直流电流为：\n'
    Real(Idc,**Option)
    Option['string'] = '功角：\n'
    Real(57.3*derta,**Option)
    Option['string'] = '调制比：\n'
    Real(M,**Option)
    MaxError = np.max(np.abs(Delta))
    return(U,Angle,P,Q,Pij,Qij,Vd,Id,Kt,W,fi,Pss,Qss,Pd,Udc,Idc,derta,M,Pv,Qv,Pdc,Pijdc,MaxError)