# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:19:17 2018

@author: yy
"""

import os

import logging
#backup 无效 logging.basicConfig(level=logging.DEBUG)#temp log
logger = logging.getLogger("developing")#debug log

#执行cmd
def cmd(cmdString):
    logger.debug(cmdString.split(" "))
    import subprocess
    try:
        what=subprocess.check_output(cmdString.split(" "))
        logger.debug(what)
    except subprocess.CalledProcessError as e:
        logger.error(e.output)
        
#sgrna
def sgrnacas9(taskSpacePath,geneFilePath,referenceFilePath,sgrnaAppFilePath,nxmAppFilePath,pam_value):
    #backup 
    #import uuid
    #taskSpacePath="/media/huangteng/data_cas9/"+str(uuid.uuid5(uuid.uuid1(), 'sgrnaFolder'))+"/"
    
    #处理基因文件，加长度
    gene5lenghFilePath=os.path.join(taskSpacePath,"gene5lengh.txt")#基因文件，带基因长度
    from .computingProcess import addSeqLength
    addSeqLength(geneFilePath,gene5lenghFilePath)
    
    #目标基因中sgRNA预测
    cmdString_sgrna="perl "+sgrnaAppFilePath+" -i "+gene5lenghFilePath+" -g "+referenceFilePath+" -n 0 -o b -t s -v l -p "+taskSpacePath
    cmd(cmdString_sgrna)
	
    #sgRNA可编辑位点预测
    total_sgRNA=taskSpacePath+'sgRNAcas9.report_20.b/report_protospacer_single.txt'#输入文件
    sgRNA_editablePredict=os.path.join(taskSpacePath,"sgRNA_editable_predict.txt")#处理结果文件
    cmdString_sgrna2="perl "+nxmAppFilePath+" -i "+total_sgRNA+" -o "+sgRNA_editablePredict
    cmd(cmdString_sgrna2)
	
    #筛选包含可编辑位点的sgRNA
    sgRNA_editable=os.path.join(taskSpacePath,"sgRNA_editable.txt")#结果文件
    from .computingProcess import selectSgRNA
    selectSgRNA(sgRNA_editablePredict,sgRNA_editable)
    
    #生成off-target预测输入文件
    off_targetInput=os.path.join(taskSpacePath,"off_targetInput.txt")#结果文件，off-target输入文件
    from .computingProcess import CasOFFinderInput
    CasOFFinderInput(sgRNA_editable,off_targetInput,referenceFilePath,pam_value)
	
    #off-target预测
    Cas_OFFinderOutputPath=taskSpacePath+"CasOffinderOutput.txt"#结果文件
    cmdString_sgrna3="./componenttib/func/cas-offinder "+off_targetInput+" C "+Cas_OFFinderOutputPath
    cmd(cmdString_sgrna3)
	
    #生成结果文件.txt
    result_FilePath=os.path.join(taskSpacePath,"sgRNA_result_editable.txt")
    from .computingProcess import riskEvaluation
    from .computingProcess import hamming_distance
    riskEvaluation(sgRNA_editable,Cas_OFFinderOutputPath,result_FilePath)
    
    #生成结果文件.csv
    import csv
    sgRNA_result_editableFilePath=os.path.join(taskSpacePath,"sgRNA_result_editable.csv")#处理结果文件
    with open(sgRNA_result_editableFilePath, 'wb') as oFile:
        spamwriter = csv.writer(oFile, delimiter=str(','))
        with open(result_FilePath, 'rU') as iFile:
            for line in iFile:
                line = line.strip().split('\t')
                spamwriter.writerow(line)
    sgrna_reportFilePath=sgRNA_result_editableFilePath
    return sgrna_reportFilePath
