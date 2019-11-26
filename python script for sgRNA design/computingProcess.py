# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 15:37:33 2019

@author: yy
"""

#处理目标基因文件，加长度
def addSeqLength(inputFilePath,outputFilePath):
    content=""
    with open(outputFilePath,"w") as oFile:
        with open(inputFilePath,"r") as iFile:
            for line in iFile.readlines():
                content+=line
        contentList=content.split(">")[1:]#分割为基因
        #计数
        #backup 正则方式，也对
        #import re
        #pattern=re.compile(r'[ATGCatgc]')
        for gene in contentList:
            lengthCount=0
            name=gene.split("\n")[0].strip()
            sequences=gene.split("\n")[1:]
            for sequence in sequences:
                #print(len(pattern.findall(sequence)))#backup 正则方式，也对
                lengthCount+=len(sequence.strip())
            name=">"+str(lengthCount)+"_"+name
            oFile.write(name+"\n")
            oFile.writelines(sequences)
            oFile.write("\n")
        
#Select sgRNA with editable site
def selectSgRNA(inputFilePath,outputFilePath):
    with open(outputFilePath, 'w') as oFile:
        with open(inputFilePath, 'r') as iFile:
            for line in iFile:
                linelist = line.strip().split('\t')
                if len(linelist) ==10:
                    if int(linelist[7]) == 1:
                        lst1=linelist[0].split('_')
                        geneName=linelist[0].strip(lst1[0]).strip(lst1[-1]).strip('_').strip(lst1[-2]).strip('_')
                        sgrnaID=geneName+'_'+lst1[0]+'_'+lst1[-2]+'_'+lst1[-1]
                        pamLen=int(linelist[4])-20
                        linewithoutIDlist=linelist[1:4]+[linelist[3][:20],linelist[3][-pamLen:]]+linelist[4:]
                        linewithoutID='\t'.join(linewithoutIDlist)
                        oFile.write(sgrnaID+'\t'+linewithoutID+'\n')
                elif len(linelist) ==14:
                    if int(linelist[7]) == 1 or int(linelist[11]) == 1:
                        lst1=linelist[0].split('_')
                        geneName=linelist[0].strip(lst1[0]).strip(lst1[-1]).strip('_').strip(lst1[-2]).strip('_')
                        sgrnaID=geneName+'_'+lst1[0]+'_'+lst1[-2]+'_'+lst1[-1]
                        pamLen=int(linelist[4])-20
                        linewithoutIDlist=linelist[1:4]+[linelist[3][:20],linelist[3][-pamLen:]]+linelist[4:]
                        linewithoutID='\t'.join(linewithoutIDlist)
                        oFile.write(sgrnaID+'\t'+linewithoutID+'\n')

#Generate input file for Cas-OFFinder
def CasOFFinderInput(inputFilePath,outputFilePath,referenceFilePath,pam_value):
    with open(outputFilePath,"w") as oFile:
        with open(inputFilePath,"r") as iFile:
            oFile.write(referenceFilePath+'\n')
            PAM=pam_value.split('-')[-1]
            if PAM=='NGG':
                oFile.write('NNNNNNNNNNNNNNNNNNNNNRG'+'\n')
            else:
                oFile.write('NNNNNNNNNNNNNNNNNNNN'+PAM+'\n')
            for line in iFile.readlines():
                sgrnaSeq=line.split('\t')[3]+'\t'+'5'
                oFile.write(sgrnaSeq+'\n')

#等长度碱基序列的错配数即汉明距离（Compute Hamming distance）
def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))

#Generate result with off-target risk evaluation
def riskEvaluation(sgRNAEditableFilePath,offtargetOutputFilePath,resultFilePath):
    with open(resultFilePath,"w") as oFile:
        with open(offtargetOutputFilePath,"r") as iFile1:
            with open(sgRNAEditableFilePath,"r") as iFile2:
                title="sgRNA_ID\tStart\tEnd\tCRISPR_target_sequence(5'-3')\tsgRNA_Sequence\tPAM_Sequence\tCRISPR_target_sequence_Length(nt)\tGC%\t1st_target_codon\tEditable_for_forming_an_early_stop_codon\tPosition_of_C_inside_gRNA\tPosition_of_C_inside_the_gene (%)\t2nd_target_codon\tEditable_for_forming_an_early_stop_codon\tPosition_of_C_inside_gRNA\tPosition_of_C_inside_the_gene (%)\tM0_OT\tM1_OT\tM2_OT\tM3_OT\tM4_OT\tM5_OT\tTotal_No.of_OT\tM0_POT\tM1_POT\tM2_POT\tM3_POT\tM4_POT\tM5_POT\tTotal_No.of_POT\tRisk_evaluation\n"
                oFile.write(title)
                offTargetEvaluation={}
                for line1 in iFile1:
                    sgRNA=line1.split('\t')[0]
                    offTargetSeq=line1.split('\t')[3]
                    OT=sgRNA[:20]
                    offOT=offTargetSeq[:20]
                    mismatch_OT=hamming_distance(OT,offOT)
                    M_OT='M'+str(mismatch_OT)+'_POT'
                    POT=sgRNA[:8]
                    offPOT=offTargetSeq[:8]
                    mismatch_POT=hamming_distance(POT,offPOT)
                    M_POT='M'+str(mismatch_POT)+'_POT2'
                    seed=sgRNA[8:20]
                    offseed=offTargetSeq[8:20]
                    mismatch_seed=hamming_distance(seed,offseed)
                    if sgRNA in offTargetEvaluation:
                        offTargetEvaluation[sgRNA][M_OT]+=1
                        if mismatch_seed==0:
                            offTargetEvaluation[sgRNA][M_POT]+=1
                    else:
                        sgrnaOffTarget={'M0_POT':0,'M1_POT':0,'M2_POT':0,'M3_POT':0,'M4_POT':0,'M5_POT':0,'Total_No.of_OT':0,'M0_POT2':0,'M1_POT2':0,'M2_POT2':0,'M3_POT2':0,'M4_POT2':0,'M5_POT2':0,'Total_No.of_SOT':0,'Risk_evaluation':''}
                        sgrnaOffTarget[M_OT]+=1
                        if mismatch_seed==0:
                            sgrnaOffTarget[M_POT]+=1
                        offTargetEvaluation[sgRNA]=sgrnaOffTarget
                for el in offTargetEvaluation:
                    TotalNoofOT=int(offTargetEvaluation[el]['M0_POT'])+int(offTargetEvaluation[el]['M1_POT'])+int(offTargetEvaluation[el]['M2_POT'])+int(offTargetEvaluation[el]['M3_POT'])+int(offTargetEvaluation[el]['M4_POT'])+int(offTargetEvaluation[el]['M5_POT'])-1
                    offTargetEvaluation[el]['Total_No.of_OT']+=TotalNoofOT
                    TotalNoofSOT=int(offTargetEvaluation[el]['M0_POT2'])+int(offTargetEvaluation[el]['M1_POT2'])+int(offTargetEvaluation[el]['M2_POT2'])+int(offTargetEvaluation[el]['M3_POT2'])+int(offTargetEvaluation[el]['M4_POT2'])+int(offTargetEvaluation[el]['M5_POT2'])-1
                    offTargetEvaluation[el]['Total_No.of_SOT']+=TotalNoofSOT
                    if offTargetEvaluation[el]['M0_POT']==0:
                        offTargetEvaluation[el]['Risk_evaluation']+='Discard'
                    elif offTargetEvaluation[el]['M0_POT']>=2:
                        offTargetEvaluation[el]['Risk_evaluation']+='Repeat_sites_or_bad'
                    elif offTargetEvaluation[el]['M1_POT']==0 and offTargetEvaluation[el]['M1_POT2']==0 and offTargetEvaluation[el]['M2_POT']==0 and offTargetEvaluation[el]['M2_POT2']==0 and offTargetEvaluation[el]['M3_POT2']==0 and offTargetEvaluation[el]['M4_POT2']==0 and offTargetEvaluation[el]['M5_POT2']==0:
                        offTargetEvaluation[el]['Risk_evaluation']+='Best'
                    elif offTargetEvaluation[el]['M1_POT']==0 and offTargetEvaluation[el]['M1_POT2']==0 and offTargetEvaluation[el]['M2_POT']==0 and offTargetEvaluation[el]['M2_POT2']==0 and offTargetEvaluation[el]['M3_POT2']>=0 and offTargetEvaluation[el]['M4_POT2']>=0 and offTargetEvaluation[el]['M5_POT2']>=0:
                        offTargetEvaluation[el]['Risk_evaluation']+='Low_risk'
                    elif offTargetEvaluation[el]['M1_POT']>=1 and offTargetEvaluation[el]['M1_POT2']>=1:
                        offTargetEvaluation[el]['Risk_evaluation']+='High_risk'
                    elif offTargetEvaluation[el]['M2_POT']>=1 and offTargetEvaluation[el]['M2_POT2']>=1:
                        offTargetEvaluation[el]['Risk_evaluation']+='Moderate_risk'
                    else:
                        offTargetEvaluation[el]['Risk_evaluation']+='Moderate_risk'
                for line2 in iFile2:
                    selectedsgrna=line2[:-1]
                    line3=selectedsgrna.split('\t')
                    sgrnaSeq2=line3[3]
                    if sgrnaSeq2 in offTargetEvaluation:
                        if len(line3)==12:
                            offtarget1=''
                            el2=offTargetEvaluation[sgrnaSeq2]
                            offtarget1+=('\t'+str(el2['M0_POT'])+'\t'+str(el2['M1_POT'])+'\t'+str(el2['M2_POT'])+'\t'+str(el2['M3_POT'])+'\t'+str(el2['M4_POT'])+'\t'+str(el2['M5_POT'])+'\t'+str(el2['Total_No.of_OT'])+'\t'+str(el2['M0_POT2'])+'\t'+str(el2['M1_POT2'])+'\t'+str(el2['M2_POT2'])+'\t'+str(el2['M3_POT2'])+'\t'+str(el2['M4_POT2'])+'\t'+str(el2['M5_POT2'])+'\t'+str(el2['Total_No.of_SOT'])+'\t'+str(el2['Risk_evaluation']))
                            oFile.write(selectedsgrna+'\t\t\t\t'+offtarget1+'\n')
                        else:
                            offtarget2=''
                            el3=offTargetEvaluation[sgrnaSeq2]
                            offtarget2+=('\t'+str(el3['M0_POT'])+'\t'+str(el3['M1_POT'])+'\t'+str(el3['M2_POT'])+'\t'+str(el3['M3_POT'])+'\t'+str(el3['M4_POT'])+'\t'+str(el3['M5_POT'])+'\t'+str(el3['Total_No.of_OT'])+'\t'+str(el3['M0_POT2'])+'\t'+str(el3['M1_POT2'])+'\t'+str(el3['M2_POT2'])+'\t'+str(el3['M3_POT2'])+'\t'+str(el3['M4_POT2'])+'\t'+str(el3['M5_POT2'])+'\t'+str(el3['Total_No.of_SOT'])+'\t'+str(el3['Risk_evaluation']))
                            oFile.write(selectedsgrna+offtarget2+'\n')
                    else:
                        oFile.write(selectedsgrna+'\n')
