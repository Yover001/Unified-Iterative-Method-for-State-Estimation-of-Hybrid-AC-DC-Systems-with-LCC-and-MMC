import datetime
import numpy as np
from math import pi,sqrt
from numpy import arctan,sin,cos,arccos
from SE import ADSE
from Polar import AC_Polar,LCC_Polar,VSC_Polar
from Data import GetNodeData,GetLineData,LCC_GetNodeData,VSC_GetNodeData,DC_GetLineData,GetY,GetYdc
from OutTxt import SingleTxt,StringTxt,Real
Out_Path = 'Cal-Process.txt'
starttime = datetime.datetime.now()
#------------------------------------------网络信息读取------------------------------------------------#
# 获取直流系统信息
RootPath = 'C:\\Users\\lenovo\\Desktop\\python代码\\状态估计\\SE-统一\\Data\\'
VSC_Node = RootPath+'VSC_NodeData.txt'
LCC_Node = RootPath+'LCC_NodeData.txt'    
DC_Line = RootPath+'DC_LineData.txt'
VSC_NodeData = VSC_GetNodeData(VSC_Node,show=1)
LCC_NodeData = LCC_GetNodeData(LCC_Node,show=1)
DC_LineData = DC_GetLineData(DC_Line,show=1)
Ydc=GetYdc(LCC_NodeData,VSC_NodeData,DC_LineData,path=Out_Path,width=6)
#获取交流系统信息
FilePath_Node = RootPath+'NodeData.txt'    
FilePath_Line = RootPath+'LineData.txt'
NodeData = GetNodeData(FilePath_Node,show=1)
LineData = GetLineData(FilePath_Line,show=1)
Y = GetY(NodeData,LineData,path=Out_Path,width=6) 

AC_LCC=LCC_NodeData[:,0]                            #相连的交流节点(LCC)
AC_LCC=AC_LCC.astype(int)-1
AC_VSC=VSC_NodeData[:,0]                            #相连的交流节点(VSC)
AC_VSC=AC_VSC.astype(int)-1
#----------------------------------------------交直流状态估计----------------------------------------#
StringTxt(path=Out_Path,string='交直流状态估计',fmt='w')
S0 = '-'
for i in range(130):
    S0 = S0+'-' 
StringTxt(path=Out_Path,string=S0)
#-初始化   
U,Angle = AC_Polar(NodeData[:,3],NodeData[:,4],path=Out_Path,width=9)
Vd,Id,Kt,W,fi = LCC_Polar(LCC_NodeData,path=Out_Path,width=9) 
Udc,Idc,derta,M = VSC_Polar(VSC_NodeData,path=Out_Path,width=9) 
# 迭代
Iter = 0
Tol = 1e-5
MaxIter = 7
while True:
	Iter = Iter + 1
	U,Angle,P,Q,Pij,Qij,Vd,Id,Kt,W,fi,Pss,Qss,Pd,Udc,Idc,derta,M,Pv,Qv,Pdc,Pdcij,MaxError= ADSE(U,Angle,Vd,Id,Kt,W,fi,Udc,Idc,derta,M,Ydc,Y,NodeData,LineData,LCC_NodeData,VSC_NodeData,DC_LineData,AC_LCC,AC_VSC,Tol,path=Out_Path,width=9)
	print(MaxError)
	if Iter>MaxIter or MaxError<Tol:
		break
	# 结束交流循环
if MaxError<Tol:
	SingleTxt(Iter-1,path=Out_Path,string='交直流迭代完成，更新次数为：')
	SingleTxt(MaxError,path=Out_Path,string='最大误差为：')
#-------------------------------------------AC状态估计结果---------------------------------------------#
	Real(U,path=Out_Path,string=S0+'\nAC电压：\n')
	Real(P,path=Out_Path,string='注入有功：\n')
	Real(Q,path=Out_Path,string='注入有功：\n')
	Real(Pij,path=Out_Path,string='线路有功：\n')
	Real(Qij,path=Out_Path,string='线路有功：\n')
#-------------------------------------------LCC状态估计结果---------------------------------------------#
	Real(Vd,path=Out_Path,string=S0+'\nLCC直流电压：\n')
	Real(Id,path=Out_Path,string='LCC直流电流：\n')
	Real(Pss,path=Out_Path,string='交流侧有功：\n')
	Real(Qss,path=Out_Path,string='交流侧无功：\n')
	Real(Pd,path=Out_Path,string='直流功率：\n')
#-------------------------------------------VSC状态估计结果---------------------------------------------#
	Real(Udc,path=Out_Path,string=S0+'\nVSC直流电压：\n')
	Real(Idc,path=Out_Path,string='VSC直流电流：\n')
	Real(Pv,path=Out_Path,string='交流侧有功：\n')
	Real(Qv,path=Out_Path,string='交流侧无功：\n')
	Real(Pdc,path=Out_Path,string='直流功率：\n')
	Real(Pdcij,path=Out_Path,string='线路功率：\n')
else:
	SingleTxt(MaxError,path=Out_Path,string='结果不收敛!')
	print ('不收敛')
# TIME
endtime = datetime.datetime.now()
print (endtime - starttime)
