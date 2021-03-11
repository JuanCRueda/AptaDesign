#!/usr/bin/env python3

import pandas as pd
from subprocess import Popen, PIPE, run, STDOUT
from random import randint, choice
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import sys
from os import path, mkdir
import logging

def AptaDesign(fasta_file='',n_conserved_seqs=8,Hybridation_target='',n_pool=100,n_gen=1000,n_candidates=10,min_length=10,hyperdiverse_generations=3,max_consecutive_hyperdiverse=10,max_consecutive_score=10,visualize='True',break_score=0.99,output_path='./',output_name='Final_candidates'):
    print('Welcome to AptaDesign v1.0.2-l')
    print('------------------------------------------')
    print('If you use this software tool for a publication, please cite:')
    print('Rueda-Silva, Juan Carlos. (2021). AptaDesign')
    print('------------------------------------------')
    print('Copyright: 2021 Juan Carlos Rueda Silva')
    print('------------------------------------------')
    if not(path.isfile(fasta_file)):
        print('Please provide the path to a valid fasta/fastq file')
        sys.exit()
    if not(path.isdir(output_path)):
        print('Please provide a valid directory for saving the output')
    if  not(('.fasta' in fasta_file) or ('.fastq' in fasta_file)):
        print('No fasta or fastq file provided, please try again')
        sys.exit()
    save_path=output_path+'/'+output_name
    if path.isdir(save_path):
        print('Please choose another name, this name is already used')
        sys.exit()
    n_conserved_seqs=check_num_input(n_conserved_seqs)
    if not(n_conserved_seqs):
        print('Please enter a valid number of mottifs to use as evaluation')
        sys.exit()
    n_pool=check_num_input(n_pool)
    if not(n_pool):
        print('Please enter a valid pool size')
        sys.exit()
    n_gen=check_num_input(n_gen)
    if not(n_gen):
        print('Please enter a valid number of generations')
        sys.exit()
    n_candidates=check_num_input(n_candidates)
    if not(n_candidates):
        print('Please enter a valid number of candidates to select')
        sys.exit()
    min_length=check_num_input(min_length)
    if not(min_length):
        print('Please enter a valid minimum length')
        sys.exit()
    hyperdiverse_generations=check_num_input(hyperdiverse_generations)
    if not(hyperdiverse_generations):
        print('Please enter a valid number of subgenerations per hyperdiverse period')
        sys.exit()
    max_consecutive_hyperdiverse=check_num_input(max_consecutive_hyperdiverse)
    if not(max_consecutive_hyperdiverse):
        print('Please enter a valid maximum number of consecutive hyperdiverse periods')
        sys.exit()
    max_consecutive_score=check_num_input(max_consecutive_score)
    if not(max_consecutive_score):
        print('Please enter a valid maximum number of consecutive equal scores')
        sys.exit()
    break_score=check_float_input(break_score)
    if not(break_score) or break_score>1 or break_score<0:
        print('Please enter a valid break score')
        sys.exit()

    mkdir(save_path)
    logging.basicConfig(filename=save_path+'/'+output_name+'.log',filemode='w',format='%(asctime)s -%(message)s',level=logging.INFO)
    logging.info('Starting program')
    df_conserved_seqs,max_strc_MFE=Get_conserved_seqs(fasta_file,n_conserved_seqs)
    print('------------------------------------------')
    print('Begining with Aptamer Design by Directed Evolution')
    logging.info('Begining with Aptamer Design by Directed Evolution')
    print('------------------------------------------')
    print('Defining evaluation equation')
    logging.info('Defining evaluation equation')
    total_motif_score=df_conserved_seqs['Score'].sum()
    df_conserved_seqs['Weights']=df_conserved_seqs['Score']/total_motif_score
    print('Evaluation equation defined')
    logging.info('Evaluation equation defined')
    print('------------------------------------------')
    print('Generating starting pool')
    logging.info('Generating starting pool')
    print('------------------------------------------')
    if bool(Hybridation_target):
        target_seq=DNA_check(Hybridation_target.lower())
        if bool(target_seq):
            target_dna=Seq(target_seq)
            starting_sequence=str(target_dna.reverse_complement())
            max_MFE=abs(MFE_Hybridization(starting_sequence,target_seq))
            pool=initial_pool_with_seq(n_pool,starting_sequence,min_length)
        else:
            print('Target sequence provided is not valid, please try again')
            logging.info('Target sequence provided is not valid, please try again')
            sys.exit()
    else:
        pool=initial_pool_gen(n_pool,min_length)
    candidates=pd.DataFrame({'sequences':['','','','','','','','','',''],'score':[0,0,0,0,0,0,0,0,0,0]})
    print('Starting pool generated')
    logging.info('Starting pool generated')
    print('------------------------------------------')
    if visualize in ['True','true','T','t','Yes','yes','Y','y','TRUE','YES',True]:
        visualize=True
    else:
        visualize=False
    if visualize:
        generations=[]
        scores=[]
        line=[]
        line2=[]
    print('Starting Evolution')
    logging.info('Starting Evolution')
    print('------------------------------------------')
    c_hyperdiverse=0
    score=0
    c_score=0
    for generation in range(n_gen):
        logging.info('Score: '+str(score))
        print('Current generation: '+str(generation+1)+' of '+str(n_gen))
        logging.info('Current generation: '+str(generation+1)+' of '+str(n_gen))
        candidates_past=candidates
        past_score=score
        if bool(Hybridation_target):
            candidates=pool_evaluation_with_target(pool,n_candidates,df_conserved_seqs,max_strc_MFE,target_seq,max_MFE)
        else:
            candidates=pool_evaluation(pool,n_candidates,df_conserved_seqs,max_strc_MFE)
        score=np.max(candidates.score)
        if visualize:
            generations.append(generation+1)
            scores.append(score)
            line=plotter(generations,scores,line,pause_time=0.1)
        if score==past_score:
            c_score+=1
            if c_score>max_consecutive_score:
                print('------------------------------------')
                print('Detected maximum with low possibilities of overcoming')
                print('Stopping Evolution')
                logging.info('Detected maximum with low possibilities of overcoming')
                logging.info('Stopping Evolution')
                print('------------------------------------')
                break
        else:
            c_score=0
        if score>break_score:
            print('------------------------------------')
            print('Minimum break score reached')
            print('Stopping evolution')
            logging.info('Minimum break score reached')
            logging.info('Stopping evolution')
            print('------------------------------------')
            break
        if candidate_similarities(candidates,candidates_past)>0.7:
            c_hyperdiverse+=1
            if c_hyperdiverse>max_consecutive_hyperdiverse:
                print('------------------------------------')
                print('Reached plateau unlikely to overcome to the current settings')
                print('Stopping Evolution')
                logging.info('Reached plateau unlikely to overcome to the current settings')
                logging.info('Stopping Evolution')
                print('------------------------------------')
                break
            logging.info('Starting hyperdiverse period')
            if bool(Hybridation_target):
                pool=new_pool_explosion_with_target(candidates,n_pool,n_candidates,df_conserved_seqs,hyperdiverse_generations,min_length,max_strc_MFE,target_seq,max_MFE)
                candidates=pool_evaluation_with_target(pool,n_candidates,df_conserved_seqs,max_strc_MFE,target_seq,max_MFE)
                pool=new_pool_gen(candidates,n_pool,n_candidates,min_length)
            else:
                pool=new_pool_explosion(candidates,n_pool,n_candidates,df_conserved_seqs,hyperdiverse_generations,min_length,max_strc_MFE)
                candidates=pool_evaluation(pool,n_candidates,df_conserved_seqs,max_strc_MFE)
                pool=new_pool_gen(candidates,n_pool,n_candidates,min_length)
            logging.info('End of hyperdiverse period')
        else:
            c_hyperdiverse=0
            pool=new_pool_gen(candidates,n_pool,n_candidates,min_length)
    if bool(Hybridation_target):
        final_candidates=pool_evaluation_with_target(pool,n_candidates,df_conserved_seqs,max_strc_MFE,target_seq,max_MFE)
    else:
        final_candidates=pool_evaluation(pool,n_candidates,df_conserved_seqs,max_strc_MFE)
    logging.info('Creating output files')
    final_candidates.to_excel(save_path+'/'+output_name+'_results.xlsx')
    if visualize:
        plt.savefig(save_path+'/'+output_name+'_fig.png',bbox_inches='tight')
    print('------------------------------------')
    print('Succesfully executed!')
    print('------------------------------------')
    print('Output located in: '+save_path)
    print('------------------------------------')
    print('Thank you for using AptaDesign!')


def Get_conserved_seqs(fasta_file,n_conserved_seqs=8):
    print('Selecting conserved motifs')
    logging.info('Selecting conserved motifs')
    print('------------------------------------------')
    if '.fasta' in fasta_file:
        df_fasta=fasta_to_df(fasta_file)
    else:
        df_fasta=fastq_to_df(fasta_file)
    print('Generating pool of possible motifs')
    logging.info('Generating pool of possible motifs')
    print('------------------------------------------')
    max_strc_MFE=abs(np.min(df_fasta['MFE']))
    df_seqs=Get_Candidate_seqs(df_fasta)
    print('Pool of possible motifs generated')
    logging.info('Pool of possible motifs generated')
    print('------------------------------------------')
    df_conserved_seqs=Conserved_Seq_Evaluate(df_fasta,df_seqs,n_conserved_seqs)
    print('Selected conserved motifs:')
    print(df_conserved_seqs.loc[:,['Sequence','Structure']])
    logging.info('Selected conserved motifs:')
    logging.info('Selected conserved motifs:')
    print('------------------------------------------')
    return df_conserved_seqs, max_strc_MFE

def fasta_to_df(fasta_file):
    ids=[]
    seqs=[]
    strs=[]
    MFEs=[]
    with open(fasta_file) as fasta:
        while True:
            line=fasta.readline().rstrip()
            if '>' in line:
                seq=fasta.readline().rstrip()
                seq=DNA_check(seq)
                if bool(seq):
                    ids.append(line[1:])
                    seqs.append(seq.lower())
                    strc,MFE=Structure_Aptamer(seq)
                    strs.append(strc)
                    MFEs.append(MFE)
            if len(line)==0:
                break
    if len(seqs)==0:
        print('No valid sequences detected in the provided fasta file, please try again')
        sys.exit()
    dict_fasta={'Ids':ids,'Sequence':seqs,'Structure':strs,'MFE':MFEs}
    df_fasta=pd.DataFrame(dict_fasta)
    return df_fasta

def fastq_to_df(fasta_file):
    ids=[]
    seqs=[]
    strs=[]
    MFEs=[]
    with open(fasta_file) as fasta:
        while True:
            line=fasta.readline().rstrip()
            if '@' in line:
                seq=fasta.readline().rstrip()
                seq=DNA_check(seq)
                if bool(seq):
                    ids.append(line[1:])
                    seqs.append(seq.lower())
                    strc,MFE=Structure_Aptamer(seq)
                    strs.append(strc)
                    MFEs.append(MFE)
            if len(line)==0:
                break
    if len(seqs)==0:
        print('No valid sequences detected in the provided fasta file, please try again')
        sys.exit()
    dict_fasta={'Ids':ids,'Sequence':seqs,'Structure':strs,'MFE':MFEs}
    df_fasta=pd.DataFrame(dict_fasta)
    return df_fasta

def DNA_check(seq):
    seq=seq.lower()
    result_seq=''
    for n in seq:
        if n in ['a','t','c','g','n']:
            result_seq+=n
        elif n=='u':
            result_seq+='t'
        elif not(n==' '):
            return False
    return result_seq
        

def Build_df_seqs(df_fasta,min_conserved_length,max_conserved_length):
    seq_df=pd.DataFrame(columns=['Sequence','Present_in','Structure'])
    for ind in list(df_fasta.index.values):
        seq=df_fasta.loc[ind,'Sequence']
        strc=df_fasta.loc[ind,'Structure']
        for n in range(min_conserved_length,max_conserved_length+1):
            superior_lim=len(seq)-n
            if superior_lim<1:
                break
            else:
                for m in range(superior_lim):
                    if m+n>len(seq):
                        break
                    else:
                        sub_seq=seq[m:m+n]
                        pos=seq.find(sub_seq)
                        sub_strc=strc[pos:pos+n]
                        if ('(((' in sub_strc or ')))' in sub_strc) and '..' in sub_strc:
                            if sub_seq in list(seq_df['Sequence']):
                                index=seq_df.index[seq_df['Sequence']==sub_seq]
                                if sub_strc==str(seq_df.loc[index,'Structure']):
                                    seq_df.loc[index,'Present_in']=seq_df.loc[index,'Present_in']+1
                                else:
                                    seq_df.loc[str(ind)+'_'+str(n)+'_'+str(m)]=[sub_seq,1,strc[pos:pos+n]] 
                            else:
                                seq_df.loc[str(ind)+'_'+str(n)+'_'+str(m)]=[sub_seq,1,sub_strc] 
    return seq_df

def Get_Candidate_seqs(df_fasta):
    seq_df=pd.DataFrame(columns=['Sequence','Present_in','Structure'])
    seq_num=0
    for ind in list(df_fasta.index.values):
        seq_num+=1
        print('Obtaining motifs from sequence '+str(seq_num)+' of '+str(len(list(df_fasta.index.values))))
        logging.info('Obtaining motifs from sequence '+str(seq_num)+' of '+str(len(list(df_fasta.index.values))))
        seq=df_fasta.loc[ind,'Sequence']
        strc=df_fasta.loc[ind,'Structure']
        start_n=strc.find('(')
        c=0
        while True:
            if start_n==len(seq):
                print('Number of motifs generated for sequence '+str(seq_num)+': '+str(c))
                logging.info('Number of motifs generated for sequence '+str(seq_num)+': '+str(c))
                print('------------------------------------------')
                break
            else:
                new_s=strc[start_n:].find('...')
                if new_s==-1:
                    print('Number of motifs generated for sequence '+str(seq_num)+': '+str(c))
                    logging.info('Number of motifs generated for sequence '+str(seq_num)+': '+str(c))
                    print('------------------------------------------')
                    break
                elif strc[new_s+start_n:].find('(')==-1 and strc[new_s+start_n:].find('.)')==-1:
                    print('Number of motifs generated for sequence '+str(seq_num)+': '+str(c))
                    logging.info('Number of motifs generated for sequence '+str(seq_num)+': '+str(c))
                    print('------------------------------------------')
                    break
            for n in range(start_n,len(seq)):
                if '...'==strc[n:n+3]:
                    c+=1
                    v2_end=n
                    n_start_1=strc[n:].find('(')
                    n_start_2=strc[n:].find(')')
                    if n_start_1==-1 or n_start_1>n_start_2:
                        next_start=n_start_2+n
                    else:
                        next_start=n_start_1+n
                    break
            sub_seq=seq[start_n:next_start]
            sub_strc=strc[start_n:next_start]
            if sub_seq in list(seq_df['Sequence']):
                index=seq_df.index[seq_df['Sequence']==sub_seq]
                ind2=index[0]
                if sub_strc==str(seq_df.loc[ind2,'Structure']):
                    seq_df.loc[ind2,'Present_in']=seq_df.loc[ind2,'Present_in']+1
                else:
                    seq_df.loc[str(ind)+'_'+str(c)]=[sub_seq,1,sub_strc]
            else:
                seq_df.loc[str(ind)+'_'+str(c)]=[sub_seq,1,sub_strc]
            v2_start_1=strc[:start_n-1].rfind('(')
            v2_start_2=strc[:start_n-1].rfind(')')
            if  v2_start_2>v2_start_1:
                v2_start=v2_start_2
            elif v2_start_1>v2_start_2:
                v2_start=v2_start_1
            else:
                v2_start=False
            if v2_start!=False:
                sub_seq_2=seq[v2_start+1:v2_end]
                sub_strc_2=strc[v2_start+1:v2_end]
                if len(sub_seq_2)>0:
                    c+=1
                    if sub_seq_2 in list(seq_df['Sequence']):
                        index=seq_df.index[seq_df['Sequence']==sub_seq_2]
                        ind2=index[0]
                        if sub_strc_2==str(seq_df.loc[ind2,'Structure']):
                            seq_df.loc[index,'Present_in']=seq_df.loc[ind2,'Present_in']+1
                        else:
                            seq_df.loc[str(ind)+'_'+str(c)+'_2']=[sub_seq_2,1,sub_strc_2]
                    else:
                        seq_df.loc[str(ind)+'_'+str(c)+'_2']=[sub_seq_2,1,sub_strc_2]
            
            start_n=next_start
    print('------------------------------------------')
    return seq_df
            

def Conserved_Seq_Evaluate(df_fasta,df_seqs,n_conserved_seqs):
    Edits_seq=[]
    Edits_str=[]
    c=0
    print('Total candidate motifs: '+str(len(list(df_seqs.index.values))))
    logging.info('Total candidate motifs: '+str(len(list(df_seqs.index.values))))
    for ind in list(df_seqs.index.values):
        c+=1
        print('Evaluating motif: '+str(c)+' of '+str(len(list(df_seqs.index.values))))
        logging.info('Evaluating motif: '+str(c)+' of '+str(len(list(df_seqs.index.values))))
        edit_seqs=0
        edit_strcs=0
        sub_seq=df_seqs.loc[ind,'Sequence']
        sub_strc=df_seqs.loc[ind,'Structure']
        for ind2 in list(df_fasta.index.values):
            seq=df_fasta.loc[ind2,'Sequence']
            strc=df_fasta.loc[ind2,'Structure']
            editDseq,pos=EditDistance(sub_seq,seq)
            editDstrc,_=EditDistance(sub_strc,strc[pos:])
            edit_seqs+=editDseq
            edit_strcs+=editDstrc
        Edits_seq.append(edit_seqs)
        Edits_str.append(edit_strcs)
    df_seqs['Edit_seq']=Edits_seq
    df_seqs['Edit_structure']=Edits_str
    df_seqs['Edit_seq_norm']=1/(1+df_seqs['Edit_seq'])
    df_seqs['Edit_structure_norm']=1/(1+df_seqs['Edit_structure'])
    df_seqs['Score']=0.5*df_seqs['Edit_seq_norm']+0.5*df_seqs['Edit_structure_norm']
    return df_seqs.nlargest(n_conserved_seqs,['Score'])
            

def EditDistance(x,y):
    D=pd.DataFrame(columns=range(len(y)+1),index=range(len(x)+1))
    D[0]=range(len(x)+1)
    D.loc[0,:]=0
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            HDist=D.loc[i,j-1]+1
            VDist=D.loc[i-1,j]+1
            if x[i-1]==y[j-1] or x[i-1]=='n' or y[j-1]=='n':
                DDist=D.loc[i-1,j-1]
            else:
                DDist=D.loc[i-1,j-1]+1
            D.loc[i,j]=min(HDist,VDist,DDist)
    EditD=min(D.loc[len(x),:])
    posV=len(x)
    for n in range(len(y),-1,-1):
        posH=n
        if D.loc[posV,posH]==EditD:
            break
    while posV!=0:
        if posH!=0:
            left=D.loc[posV,posH-1]
            diag=D.loc[posV-1,posH-1]
            up=D.loc[posV-1,posH]
            if diag<=left and diag<=up:
                posV-=1
                posH-=1
            elif left<diag and left<up:
                posH-=1
            else:
                posV-=1
        else:
            posV-=1
    pos=posH
    return EditD,pos

def Structure_Aptamer(seq):
    '''Return the MFE (kcal/mol) of the RNA secondary strucure of a given sequence'''
    MFE_calc=Popen(['/usr/bin/RNAfold'], stdin=PIPE, stdout=PIPE,stderr=STDOUT,shell=True)
    Result=MFE_calc.communicate(seq.encode())
    pos=Result[0].find(b' (')
    Structure=Result[0][len(seq)+2:pos]
    pos2=Result[0].find(b")\n',")
    MFE=float(Result[0][pos+3:pos2-1])
    Str=Structure.decode('UTF-8')
    return Str.strip(),MFE

def MFE_Hybridization(seq_aptamer,seq_target):
    '''Return the MFE (kcal/mol) of the hybridization of two given RNA sequences'''
    Result=run(['/usr/local/bin/RNAhybrid'+' -d xi '+seq_aptamer+' '+seq_target],shell=True,stdin=PIPE, stdout=PIPE,stderr=STDOUT)
    Result=Result.stdout.decode('utf-8')
    pos1=Result.find('mfe: ')
    pos2=Result.find(' kcal/mol')
    while pos1==-1 or pos2==-1:
        Result=run(['/usr/local/bin/RNAhybrid'+' -d xi '+seq_aptamer+' '+seq_target],shell=True,stdin=PIPE, stdout=PIPE,stderr=STDOUT)
        Result=Result.stdout.decode('utf-8')
        pos1=Result.find('mfe: ')
        pos2=Result.find(' kcal/mol')
    return float(Result[pos1+4:pos2])

def initial_pool_gen(n_pool,min_length):
    seqs=[]
    for n in range(n_pool):
        seq=''
        if min_length!=False:
            for m in range(min_length):
                seq=seq+choice('atgc')
        else:
            length=randint(12,100)
            for m in range(length):
                seq=seq+choice('atgc')
        seqs.append(seq)
    return pd.DataFrame(seqs,columns=['sequences'])

def initial_pool_with_seq(n_pool,seq,min_length):
    seq=seq.lower()
    seqs=[seq]
    for n in range(n_pool-1):
        new_seq=mutation(seq)
        while new_seq in seqs:
            new_seq=mutation(new_seq)
        seqs.append(new_seq)
    return pd.DataFrame(seqs,columns=['sequences'])

def pool_evaluation(pool,n_candidates,df_conserved_seqs,max_strc_MFE):
    pool['score']=0
    for ind in list(pool.index.values):
        candidate_str,candidate_str_MFE=Structure_Aptamer(pool.loc[ind,'sequences'])
        for ind2 in list(df_conserved_seqs.index.values):
            EditDist_seq,pos=EditDistance(df_conserved_seqs.loc[ind2,'Sequence'],pool.loc[ind,'sequences'])
            EditDist_strc,_=EditDistance(df_conserved_seqs.loc[ind2,'Structure'],candidate_str[pos:pos+len(df_conserved_seqs.loc[ind2,'Structure'])+1])
            pool.loc[ind,'EditDist_seq_'+str(ind2)]=EditDist_seq
            pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_seq_'+str(ind2)])
            pool.loc[ind,'EditDist_strc_'+str(ind2)]=EditDist_strc
            pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_strc_'+str(ind2)])
            pool.loc[ind,'score']=pool.loc[ind,'score']+(0.4*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm'])+(0.6*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm'])
        pool.loc[ind,'strc_MFE']=candidate_str_MFE
        if candidate_str_MFE<-max_strc_MFE:
            adj_MFE=candidate_str_MFE-(candidate_str_MFE+max_strc_MFE)
            pool.loc[ind,'strc_MFE_norm']=abs(adj_MFE)/max_strc_MFE
        elif candidate_str_MFE>0:
            pool.loc[ind,'strc_MFE_norm']=0
        else:
            pool.loc[ind,'strc_MFE_norm']=abs(candidate_str_MFE)/max_strc_MFE
    pool['score']=0.85*pool['score']+0.15*pool['strc_MFE_norm']
    return pool.nlargest(n_candidates,['score'])

def pool_evaluation_with_target(pool,n_candidates,df_conserved_seqs,max_strc_MFE,target_seq,max_MFE):
    pool['dMFE']=0
    pool['score']=0
    for ind in list(pool.index.values):
        candidate_str,candidate_str_MFE=Structure_Aptamer(pool.loc[ind,'sequences'])
        Hybrid_MFE=MFE_Hybridization(pool.loc[ind,'sequences'],target_seq)
        if Hybrid_MFE<0:
            if Hybrid_MFE<candidate_str_MFE:
                pool.loc[ind,'dMFE']=abs(abs(Hybrid_MFE)-abs(candidate_str_MFE))
            else:
                pool.loc[ind,'dMFE']=0.0
        else:
            pool.loc[ind,'dMFE']=0.0
        for ind2 in list(df_conserved_seqs.index.values):
            EditDist_seq,pos=EditDistance(df_conserved_seqs.loc[ind2,'Sequence'],pool.loc[ind,'sequences'])
            EditDist_strc,_=EditDistance(df_conserved_seqs.loc[ind2,'Structure'],candidate_str[pos:pos+len(df_conserved_seqs.loc[ind2,'Structure'])+1])
            pool.loc[ind,'EditDist_seq_'+str(ind2)]=EditDist_seq
            pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_seq_'+str(ind2)])
            pool.loc[ind,'EditDist_strc_'+str(ind2)]=EditDist_strc
            pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_strc_'+str(ind2)])
            pool.loc[ind,'score']=pool.loc[ind,'score']+(0.4*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm'])+(0.6*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm'])
        pool.loc[ind,'strc_MFE']=candidate_str_MFE
        if candidate_str_MFE<-max_strc_MFE:
            adj_MFE=candidate_str_MFE-(candidate_str_MFE+max_strc_MFE)
            pool.loc[ind,'strc_MFE_norm']=abs(adj_MFE)/max_strc_MFE
        elif candidate_str_MFE>0:
            pool.loc[ind,'strc_MFE_norm']=0
        else:
            pool.loc[ind,'strc_MFE_norm']=abs(candidate_str_MFE)/max_strc_MFE
    pool['dMFE_norm']=pool['dMFE']/max_MFE
    pool['score']=0.85*pool['score']+0.15*pool['strc_MFE_norm']
    pool['score']=0.85*pool['score']+0.15*pool['dMFE_norm']
    return pool.nlargest(n_candidates,['score'])

def plotter(generations,scores,line,pause_time=30):
    if line==[]:
        plt.ion()
        fig=plt.figure(figsize=[6.4,6.0])
        ax1=fig.add_subplot(1,1,1)
        line,=ax1.plot(generations,scores)
        ax1.set_ylabel('Max. score')
        ax1.set_xlabel('Generations')
        plt.show()
    line.set_data(generations,scores)
    fig_current=plt.gcf()
    ax_list=fig_current.get_axes()
    if np.min(scores)<=line.axes.get_ylim()[0] or np.max(scores)>=line.axes.get_ylim()[1]:
        ax_list[0].set_ylim([np.min(scores)-0.05,np.max(scores)+0.05])
    if np.max(generations)>=line.axes.get_xlim()[1]:
        ax_list[0].set_xlim([0,np.max(generations)+0.05*np.max(generations)])
    plt.pause(pause_time)
    return line

def mutation(seq):
    seq_list=list(seq)
    mutation_type=randint(0,2)
    n_to_mutate=randint(0,len(seq)-1)
    if mutation_type==0:
        del seq_list[n_to_mutate]
    elif mutation_type==1:
        seq_list.insert(n_to_mutate,choice('atgc'))
    else:
        if seq_list[n_to_mutate]=='a':
            seq_list[n_to_mutate]=choice('tgc')
        elif seq_list[n_to_mutate]=='t':
            seq_list[n_to_mutate]=choice('agc')
        elif seq_list[n_to_mutate]=='g':
            seq_list[n_to_mutate]=choice('atc')
        else:
            seq_list[n_to_mutate]=choice('atg')
    return ''.join(seq_list)

def n_addition_mut(seq):
    seq_list=list(seq)
    n_to_mutate=randint(0,len(seq)-1)
    seq_list.insert(n_to_mutate,choice('atgc'))
    return ''.join(seq_list)

def candidate_similarities(candidates,candidates_past):
    counter=0
    past_list=list(candidates_past.sequences)
    for seq in list(candidates.sequences):
        if seq in past_list:
            counter+=1
    return float(counter)/float(len(past_list))

def new_pool_gen(candidates,n_pool,n_candidates,min_length):
    div=int((n_pool/n_candidates)-1)
    seqs=[]
    for seq in list(candidates.sequences):
        if seq not in seqs:
            seqs.append(seq)
            div_s=div
        else:
            div_s=div+1
        for n in range(div_s):
            new_seq=mutation(seq)
            while len(new_seq)<min_length:
                new_seq=n_addition_mut(new_seq)
            while new_seq in seqs:
                new_seq=mutation(new_seq)
                while len(new_seq)<min_length:
                    new_seq=n_addition_mut(new_seq)
            seqs.append(new_seq)
    return pd.DataFrame(seqs,columns=['sequences'])

def new_pool_explosion(candidates,n_pool,n_candidates,df_conserved_seqs,hyperdiverse_generations,min_length,max_strc_MFE):
    print('------------------------------------------')
    print('Starting hyperdiverse subgenerations to liberate plateau')
    div=int((n_pool/n_candidates)-1)
    seqs=list(candidates.sequences)
    seqs_len=[len(x) for x in seqs]
    print('------------------------------------------')
    for n in range(hyperdiverse_generations):
        seqs_2=[]
        print('Current subgeneration: '+str(n+1))
        for seq in seqs:
            seqs_2.append(seq)
            for n in range(div):
                new_seq=mutation(seq)
                while len(new_seq)<min_length:
                    new_seq=n_addition_mut(new_seq)
                while (new_seq in seqs) or (new_seq in seqs_2):
                    new_seq=mutation(new_seq)
                    while len(new_seq)<min_length:
                        new_seq=n_addition_mut(new_seq)
                seqs_2.append(new_seq)
        seqs=seqs_2[:]
    print('------------------------------------------')
    print('End of hyperdiverse period')
    print('------------------------------------------')
    pool=pd.DataFrame(seqs,columns=['sequences'])
    pool['score']=0
    print('------------------------------------------')
    print('Selecting subset from hyperdiverse period')
    print('------------------------------------------')
    for ind in list(pool.index.values):
        candidate_str,candidate_str_MFE=Structure_Aptamer(pool.loc[ind,'sequences'])
        for ind2 in list(df_conserved_seqs.index.values):
            EditDist_seq,pos=EditDistance(df_conserved_seqs.loc[ind2,'Sequence'],pool.loc[ind,'sequences'])
            EditDist_strc,_=EditDistance(df_conserved_seqs.loc[ind2,'Structure'],candidate_str[pos:pos+len(df_conserved_seqs.loc[ind2,'Structure'])+1])
            pool.loc[ind,'EditDist_seq_'+str(ind2)]=EditDist_seq
            pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_seq_'+str(ind2)])
            pool.loc[ind,'EditDist_strc_'+str(ind2)]=EditDist_strc
            pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_strc_'+str(ind2)])
            pool.loc[ind,'score']=pool.loc[ind,'score']+(0.4*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm'])+(0.6*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm'])
        pool.loc[ind,'strc_MFE']=candidate_str_MFE
        if candidate_str_MFE<-max_strc_MFE:
            adj_MFE=candidate_str_MFE-(candidate_str_MFE+max_strc_MFE)
            pool.loc[ind,'strc_MFE_norm']=abs(adj_MFE)/max_strc_MFE
        elif candidate_str_MFE>0:
            pool.loc[ind,'strc_MFE_norm']=0
        else:
            pool.loc[ind,'strc_MFE_norm']=abs(candidate_str_MFE)/max_strc_MFE
    pool['score']=0.85*pool['score']+0.15*pool['strc_MFE_norm']
    print('------------------------------------------')
    print('Subset created, return to main algorithm')
    print('------------------------------------------')
    return pool.nlargest(n_pool,['score'])

def new_pool_explosion_with_target(candidates,n_pool,n_candidates,df_conserved_seqs,hyperdiverse_generations,min_length,max_strc_MFE,target_seq,max_MFE):
    print('------------------------------------------')
    print('Starting hyperdiverse subgenerations to liberate plateau')
    div=int((n_pool/n_candidates)-1)
    seqs=list(candidates.sequences)
    seqs_len=[len(x) for x in seqs]
    print('------------------------------------------')
    for n in range(hyperdiverse_generations):
        seqs_2=[]
        print('Current subgeneration: '+str(n+1))
        for seq in seqs:
            seqs_2.append(seq)
            for _ in range(div):
                new_seq=mutation(seq)
                while len(new_seq)<min_length:
                    new_seq=n_addition_mut(new_seq)
                while (new_seq in seqs) or (new_seq in seqs_2):
                    new_seq=mutation(new_seq)
                    while len(new_seq)<min_length:
                        new_seq=n_addition_mut(new_seq)
                seqs_2.append(new_seq)
        seqs=seqs_2[:]
    print('------------------------------------------')
    print('End of hyperdiverse period')
    print('------------------------------------------')
    pool=pd.DataFrame(seqs,columns=['sequences'])
    pool['score']=0
    print('------------------------------------------')
    print('Selecting subset from hyperdiverse period')
    print('------------------------------------------')
    for ind in list(pool.index.values):
        candidate_str,candidate_str_MFE=Structure_Aptamer(pool.loc[ind,'sequences'])
        Hybrid_MFE=MFE_Hybridization(pool.loc[ind,'sequences'],target_seq)
        if Hybrid_MFE<0:
            if Hybrid_MFE<candidate_str_MFE:
                pool.loc[ind,'dMFE']=abs(abs(Hybrid_MFE)-abs(candidate_str_MFE))
            else:
                pool.loc[ind,'dMFE']=0.0
        else:
            pool.loc[ind,'dMFE']=0.0
        for ind2 in list(df_conserved_seqs.index.values):
            EditDist_seq,pos=EditDistance(df_conserved_seqs.loc[ind2,'Sequence'],pool.loc[ind,'sequences'])
            EditDist_strc,_=EditDistance(df_conserved_seqs.loc[ind2,'Structure'],candidate_str[pos:pos+len(df_conserved_seqs.loc[ind2,'Structure'])+1])
            pool.loc[ind,'EditDist_seq_'+str(ind2)]=EditDist_seq
            pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_seq_'+str(ind2)])
            pool.loc[ind,'EditDist_strc_'+str(ind2)]=EditDist_strc
            pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm']=1/(1+pool.loc[ind,'EditDist_strc_'+str(ind2)])
            pool.loc[ind,'score']=pool.loc[ind,'score']+(0.35*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_seq_'+str(ind2)+'_norm'])+(0.5*df_conserved_seqs.loc[ind2,'Weights']* pool.loc[ind,'EditDist_strc_'+str(ind2)+'_norm'])
        pool.loc[ind,'strc_MFE']=candidate_str_MFE
        if candidate_str_MFE<-max_strc_MFE:
            adj_MFE=candidate_str_MFE-(candidate_str_MFE+max_strc_MFE)
            pool.loc[ind,'strc_MFE_norm']=abs(adj_MFE)/max_strc_MFE
        elif candidate_str_MFE>0:
            pool.loc[ind,'strc_MFE_norm']=0
        else:
            pool.loc[ind,'strc_MFE_norm']=abs(candidate_str_MFE)/max_strc_MFE
    pool['dMFE_norm']=pool['dMFE']/max_MFE
    pool['score']=0.85*pool['score']+0.15*pool['strc_MFE_norm']
    pool['score']=0.85*pool['score']+0.15*pool['dMFE_norm']
    print('------------------------------------------')
    print('Subset created, return to main algorithm')
    print('------------------------------------------')
    return pool.nlargest(n_pool,['score'])

def check_num_input(inpt):
    try:
        inpt=int(inpt)
        return inpt
    except ValueError:
        return False
    
def check_float_input(inpt):
    try:
        inpt=float(inpt)
        return inpt
    except ValueError:
        return False

def main():
    print('Welcome to Aptadesigner v1.0.2-l')
    print('------------------------------------------')
    print('Press enter to use the defaults')
    print('Please provide the path to your desired fasta or fastq file')
    while True:
        fasta_file=input()
        if not(path.isfile(fasta_file)):
            print('Please provide the path to a valid fasta/fastq file')
        elif not(('.fasta' in fasta_file) or ('.fastq' in fasta_file)):
            print('No fasta or fastq file provided, please try again')
        else:
            break
    print('Please insert the number of motifs to use for the evaluaion (default: 8)')
    while True:
        n_conserved_seqs=input()
        if n_conserved_seqs=='':
            n_conserved_seqs=8
            break
        else:
            n_conserved_seqs=check_num_input(n_conserved_seqs)
            if not(n_conserved_seqs):
                print('Please enter a valid pool size')
            else:
                break
    
    print('Please insert the pool size to use (default: 100)')
    while True:
        n_pool=input()
        if n_pool=='':
            n_pool=100
            break
        else:
            n_pool=check_num_input(n_pool)
            if not(n_pool):
                print('Please enter a valid pool size')
            else:
                break
    print('Please insert the number of generations to simulate (default: 1000)')
    while True:
        n_gen=input()
        if n_gen=='':
            n_gen=1000
            break
        else:
            n_gen=check_num_input(n_gen)
            if not(n_gen):
                print('Please enter a valid number of generations')
            else:
                break

    print('Please insert the number of candidates to generate (default. 10)')
    while True:
        n_candidates=input()
        if n_candidates=='':
            n_candidates=10
            break
        else:
            n_candidates=check_num_input(n_candidates)
            if not(n_candidates):
                print('Please enter a valid number of candidates')
            else:
                break
    print('Please insert the maximum number of generations without score inprovement (default. 10)')
    while True:
        max_consecutive_score=input()
        if max_consecutive_score=='':
            max_consecutive_score=10
            break
        else:
            max_consecutive_score=check_num_input(max_consecutive_score)
            if not(max_consecutive_score):
                print('Please enter a valid number of generations')
            else:
                break
    print('If you desire to enter an hybridation target please enter it now, otherwise press enter')
    Hybridation_target=input()
    if Hybridation_target!='':
        Hybridation_target=DNA_check(Hybridation_target.lower())
        if not(Hybridation_target):
            while True:
                print('Please enter a valid sequence')
                Hybridation_target=input()
                if Hybridation_target!='':
                    Hybridation_target=DNA_check(Hybridation_target.lower())
                    if Hybridation_target:
                        break
                else:
                    print('No hybridation target provided')
                    break
    print('Please insert the output path to save the results')
    while True:
        output_path=input()
        if not(path.isdir(output_path)):
            print('Please insert a valid path')
        else:
            break
    print('Please give a name to this project')
    while True:
        output_name=input()
        if len(output_name)>0:
            save_path=output_path+'/'+output_name
            if path.isdir(save_path):
                print('Please choose another name, this name is already used')
            else:
                break
    print('Do you want to generate graphs for this project?')
    visualize=input('T/F: ')
    print('Open advanzed options?')
    adv_opt=input('Y/N: ')
    if adv_opt in ['True','true','T','t','Yes','yes','Y','y','TRUE','YES']:
        print('Please enter the minimum legnth of the generated apatmers (default. 10)')
        while True:
            min_length=input()
            if min_length=='':
                min_length=10
                break
            else:
                min_length=check_num_input(min_length)
                if not(min_length):
                    print('Please enter a valid length')
                else:
                    break
        print('Please enter the number of subgenerations per hyperdiverse period (default. 3)')
        while True:
            hyperdiverse_generations=input()
            if hyperdiverse_generations=='':
                hyperdiverse_generations=3
                break
            else:
                hyperdiverse_generations=check_num_input(hyperdiverse_generations)
                if not(hyperdiverse_generations):
                    print('Please enter a valid number of generations')
                else:
                    break
        print('Print enter the maximum number tolerable of consecutive hyperdiverse periods (default. 10)')
        while True:
            max_consecutive_hyperdiverse=input()
            if max_consecutive_hyperdiverse=='':
                max_consecutive_hyperdiverse=10
                break
            else:
                max_consecutive_hyperdiverse=check_num_input(max_consecutive_hyperdiverse)
                if not(max_consecutive_hyperdiverse):
                    print('Please enter a valid number of generations')
                else:
                    break
        print('Please enter the break score (between 0 and 1) (default. 0.99)')
        while True:
            break_score=input()
            if break_score=='':
                break_score=0.99
                break
            else:
                break_score=check_float_input(break_score)
                if not(break_score):
                    print('Please enter a valid score (between 0 and 1)')
                elif break_score>1 or break_score<0:
                    print('Please enter a valid score (between 0 and 1)')
                else:
                    break
        AptaDesign(fasta_file=fasta_file,n_conserved_seqs=n_conserved_seqs,break_score=break_score,min_length=min_length,hyperdiverse_generations=hyperdiverse_generations,max_consecutive_hyperdiverse=max_consecutive_hyperdiverse,Hybridation_target=Hybridation_target,n_pool=n_pool,n_gen=n_gen,n_candidates=n_candidates,max_consecutive_score=max_consecutive_score,visualize=visualize,output_path=output_path,output_name=output_name)
    else:
        AptaDesign(fasta_file=fasta_file,n_conserved_seqs=n_conserved_seqs,Hybridation_target=Hybridation_target,n_pool=n_pool,n_gen=n_gen,n_candidates=n_candidates,max_consecutive_score=max_consecutive_score,visualize=visualize,output_path=output_path,output_name=output_name)


if __name__=='__main__':
    main()
        


    
