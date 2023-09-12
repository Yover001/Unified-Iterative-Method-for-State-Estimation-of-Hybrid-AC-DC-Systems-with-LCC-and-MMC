import numpy as np
import OutTxt   
#----------------------------NodeData-------------------------#
def GetNumNode(FilePath_Node):
    with open(FilePath_Node,'r') as file:  
        NumNode = 0
        for line in file: 
            NumNode = NumNode+1  # 得到节点数目
        print('交流节点数目：',NumNode)
    return(NumNode)
def GetNodeData(FilePath_Node,*NumNode,**PlotOption):
    if not NumNode:  # 空的
        NumNode = GetNumNode(FilePath_Node)
    NodeData = np.zeros([NumNode,8]) # 初始化
    i = 0
    with open(FilePath_Node,'r') as file:  
        for line in file:
            NodeData[i,:] = np.array([float(i) for i in line.split()])
            i = i+1
    # if 'show' in PlotOption and PlotOption['show'] == 1:
    #     print('NodeData=',NodeData)
    return(NodeData) 

#----------------------------LineData------------------------#
def GetNumLine(FilePath_Line):
    with open(FilePath_Line,'r') as file:  
        NumLine = 0
        for line in file: 
            NumLine = NumLine+1  # 得到节点数目
        print('交流支路数目：',NumLine)
    return(NumLine)
def GetLineData(FilePath_Line,*NumLine,**PlotOption):
    if not NumLine:  # 空的
        NumLine = GetNumLine(FilePath_Line)
    LineData = np.zeros([NumLine,8]) # 初始化
    i = 0
    with open(FilePath_Line,'r') as file:  
        for line in file:
            LineData[i,:] = np.array([float(i) for i in line.split()])
            i = i+1
    # if 'show' in PlotOption and PlotOption['show'] == 1:
    #     print('LineData = ',LineData)
    return(LineData)

#-------------------------形成节点导纳矩阵------------------------#
def GetY(NodeData,LineData,**Option):
    NumNode = NodeData.shape[0]
    NumLine = LineData.shape[0]
    Y = np.zeros([NumNode,NumNode])+np.zeros([NumNode,NumNode])*1j
    G0 = NodeData[:,1]  # 节点对地电导
    B0 = NodeData[:,2]  # 节点对地电纳
    for i in range(NumLine):
        Node1 = int(LineData[i,0]-1)
        Node2 = int(LineData[i,1]-1)
        # print(Node1,Node2)
        R = LineData[i,2] 
        X = LineData[i,3] 
        if LineData[i,5]==0:   # 普通线路，无变压器
            B_2 = LineData[i,4] 
            Y[Node1,Node1] = Y[Node1,Node1]+B_2*1j+1/(R+1j*X) 
            Y[Node2,Node2] = Y[Node2,Node2]+B_2*1j+1/(R+1j*X) 
            Y[Node1,Node2] = Y[Node1,Node2]-1/(R+1j*X) 
            Y[Node2,Node1] = Y[Node2,Node1]-1/(R+1j*X) 
        else:  # 有变压器支路
            K  = LineData[i,5] 
            YT = 1/(R+1j*X) 
            Y[Node1,Node1] = Y[Node1,Node1]+(K-1)/K*YT+YT/K 
            Y[Node2,Node2] = Y[Node2,Node2]+(1-K)/K**2*YT+YT/K 
            Y[Node1,Node2] = Y[Node1,Node2]-1/K*YT 
            Y[Node2,Node1] = Y[Node2,Node1]-1/K*YT 
    for i in range(NumNode):
        Node = int(NodeData[i,0]-1)   # 第一列为节点编号
        Y[Node,Node] = Y[Node,Node]+G0[i]+1j*B0[i] 
    return(Y)

#----------------------------LCC_NodeData-------------------------#
def LCC_GetNumNode(DC_FilePath_Node):
    with open(DC_FilePath_Node,'r') as file: 
        DC_NumNode = 0
        for line in file: 
            DC_NumNode = DC_NumNode+1 
        print('LCC节点数目：',DC_NumNode)
    return(DC_NumNode)
def LCC_GetNodeData(DC_FilePath_Node,*DC_NumNode,**PlotOption):
    if not DC_NumNode:  # 空的
        NumNode = LCC_GetNumNode(DC_FilePath_Node)
    DC_NodeData = np.zeros([NumNode,18])
    i = 0
    with open(DC_FilePath_Node,'r') as file:  
        for line in file:
            DC_NodeData[i,:] = np.array([float(i) for i in line.split()])
            i = i+1
    # if 'show' in PlotOption and PlotOption['show'] == 1:
    #     print('LCC_NodeData=',DC_NodeData)
    return(DC_NodeData) 

#----------------------------VSC_NodeData-------------------------#
def VSC_GetNumNode(FilePath_Node):
    with open(FilePath_Node,'r') as file:  
        DC_NumNode = 0 
        for line in file: 
            DC_NumNode = DC_NumNode+1  # 得到节点数目
        print('VSC节点数目：',DC_NumNode)
    return(DC_NumNode)
def VSC_GetNodeData(FilePath_Node,*DC_NumNode,**PlotOption):
    if not DC_NumNode:  # 空的
        DC_NumNode = VSC_GetNumNode(FilePath_Node)
    DC_NodeData = np.zeros([DC_NumNode,18]) # 初始化
    i = 0
    with open(FilePath_Node,'r') as file:  
        for line in file:
            DC_NodeData[i,:] = np.array([float(i) for i in line.split()])
            i = i+1
    # if 'show' in PlotOption and PlotOption['show'] == 1:
    #     print('VSC_NodeData=',DC_NodeData)
    return(DC_NodeData) 

#----------------------------DC_LineData------------------------#
def DC_GetNumLine(FilePath_Line):
    with open(FilePath_Line,'r') as file:  
        DC_NumLine = 0
        for line in file: 
            DC_NumLine = DC_NumLine+1  # 得到节点数目
        print('直流支路数目：',DC_NumLine)
    return(DC_NumLine)
def DC_GetLineData(FilePath_Line,*DC_NumLine,**PlotOption):
    if not DC_NumLine:  # 空的
        DC_NumLine = DC_GetNumLine(FilePath_Line)
    DC_LineData = np.zeros([DC_NumLine,4]) # 初始化
    i = 0
    with open(FilePath_Line,'r') as file:  
        for line in file:
            DC_LineData[i,:] = np.array([float(i) for i in line.split()])
            i = i+1
    # if 'show' in PlotOption and PlotOption['show'] == 1:
    #     print('DC_LineData = ',DC_LineData)
    return(DC_LineData)

#-------------------------直流网络导纳矩阵------------------------#
def GetYdc(LCC_NodeData,VSC_NodeData,DC_LineData,**Option):
    DC_NumNode = VSC_NodeData.shape[0]+LCC_NodeData.shape[0]
    DC_NumLine = DC_LineData.shape[0]
    Ydc = np.zeros([DC_NumNode,DC_NumNode])
    for i in range(DC_NumLine):
        Node1 = int(DC_LineData[i,0]-1)
        Node2 = int(DC_LineData[i,1]-1)
        R = DC_LineData[i,2] 
        Ydc[Node1,Node1] = Ydc[Node1,Node1]+1/R 
        Ydc[Node2,Node2] = Ydc[Node2,Node2]+1/R 
        Ydc[Node1,Node2] = Ydc[Node1,Node2]-1/R 
        Ydc[Node2,Node1] = Ydc[Node2,Node1]-1/R 
    return(Ydc)