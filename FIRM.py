import argparse
from pssm import pssm
import pickle, gzip, os, sys, re, os, math, shutil
from copy import deepcopy
from subprocess import *
from random import sample
from multiprocessing import Pool, cpu_count, Manager
from collections import defaultdict
from miRvestigator import miRvestigator
import json


def miRNAInDict(miRNA, dict1):
    retMe = []
    for i in dict1.keys():
        if compareMiRNANames(miRNA, i):
            retMe.append(miRNAIDs[i])
    return retMe


def compareMiRNANames(a, b):
    if a==b:
        return 1
    if len(a) < len(b):
        re1 = re.compile(a+'[a-z]$')
        if re1.match(b):
            return 1
    else:
        re1 = re.compile(b+'[a-z]$')
        if re1.match(a):
            return 1
    return 0


# Function to run the meme function
def runWeeder(i):
    weeder(i)


# Run weeder and parse its output
# First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
def weeder(i=None, percTargets=50, revComp=False):
    global fastaFiles, weederPSSMs1

    seqFile = fastaFiles[i]
    print(seqFile)
    if not os.path.exists('tmp/weeder'):
        os.makedirs('tmp/weeder')

    # First run weederTFBS for 6bp motifs
    weederArgs = ' '+str(seqFile)+' HS3P small T50'
    if revComp==True:
        weederArgs += ' -S'
    errOut = open('tmp/weeder/stderr.out','w')
    weederProc = Popen("weederlauncher " + weederArgs, shell=True,stdout=PIPE,stderr=errOut)
    output = weederProc.communicate()

    # Now parse output from weeder
    PSSMs = []
    output = open(str(seqFile)+'.wee','r')
    outLines = [line for line in output.readlines() if line.strip()]
    hitBp = {}
    # Get top hit of 6bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[6] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Searching for motifs of length 8') == -1:
            break

    # Get top hit of 8bp look for "1)"
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('1) ') == -1:
            break
    hitBp[8] = outLine.strip().split(' ')[1:]

    # Scroll to where the 8bp reads wll be
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Your sequences:') == -1:
            break

    # Get into the highest ranking motifs
    seqDict = {}
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('**** MY ADVICE ****') == -1:
            break
        splitUp = outLine.strip().split(' ')
        seqDict[splitUp[1]] = splitUp[3].lstrip('>')

    # Get into the highest ranking motifs
    while 1:
        outLine = outLines.pop(0)
        if not outLine.find('Interesting motifs (highest-ranking)') == -1:
            break
    while 1:
        name = seqFile.split('/')[-1].split('.')[0] +'_'+ outLines.pop(0).strip() # Get match
        if not name.find('(not highest-ranking)') == -1:
            break
        # Get redundant motifs
        outLines.pop(0)
        redMotifs = [i for i in outLines.pop(0).strip().split(' ') if not i=='-']
        outLines.pop(0)
        outLines.pop(0)
        line = outLines.pop(0)
        instances = []
        while line.find('Frequency Matrix') == -1:
            splitUp = [i for i in line.strip().split(' ') if i]
            instances.append({'gene':seqDict[splitUp[0]],
                              'strand':splitUp[1],
                              'site':splitUp[2],
                              'start':splitUp[3],
                              'match':splitUp[4].lstrip('(').rstrip(')') })
            line = outLines.pop(0)
        # Read in Frequency Matrix
        outLines.pop(0)
        outLines.pop(0)
        matrix = []
        col = outLines.pop(0)
        while col.find('======') == -1:
            nums = [i for i in col.strip().split('\t')[1].split(' ') if i]
            colSum = 0
            for i in nums:
                colSum += int(i.strip())
            matrix += [[ float(nums[0])/float(colSum),
                         float(nums[1])/float(colSum),
                         float(nums[2])/float(colSum),
                         float(nums[3])/float(colSum)]]
            col = outLines.pop(0)
        weederPSSMs1.append(pssm(biclusterName=name,
                                 nsites=instances,eValue=hitBp[len(matrix)][1],
                                 pssm=matrix,
                                 genes=redMotifs))


def phyper(q, m, n, k):
    # Get an array of values to run
    rProc = Popen('R --no-save --slave', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    runMe = []
    for i in range(len(q)):
        runMe.append('phyper('+str(q[i])+','+str(m[i])+','+str(n[i])+','+str(k[i])+',lower.tail=F)')
    runMe = ('\n'.join(runMe)+'\n').encode()
    out = rProc.communicate(runMe)
    return [line.strip().split(' ')[1]
            for line in out[0].decode("utf-8").strip().split('\n') if line]


def clusterHypergeo(cluster):
    global db, dataset, datasetGenes, totalTargets, clusters, miRNATargetDict

    print('Cluster %d' % cluster)
    outFile = open('miRNA_'+db+'/'+str(dataset[0])+'_'+str(cluster)+'.csv','w')
    outFile.write('miRNA,Cluster.Targets,miRNA.Targets,Cluster.Genes,Total,P.Value\n')
    # k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
    # Take gene list and compute overlap with each miRNA
    allGenes = set(datasetGenes).intersection(set(totalTargets))
    genes = set(clusters[cluster]).intersection(set(allGenes))
    writeMe = []
    keys1 = miRNATargetDict.keys()
    m1s = []
    q = []
    m = []
    n = []
    k = []
    for m1 in keys1:
        m1s.append(m1)
        miRNAGenes = set(miRNATargetDict[m1]).intersection(allGenes)
        q.append(len(set(miRNAGenes).intersection(genes)))
        m.append(len(miRNAGenes))
        n.append(len(allGenes)-len(miRNAGenes))
        k.append(len(genes))
    results = phyper(q,m,n,k)
    for i in range(len(m1s)):
        writeMe.append(str(m1s[i]) + ',' + str(q[i]) + ',' + str(m[i]) + ',' + str(n[i]) + ',' + str(k[i]) + ',' + str(results[i]))
    outFile.write('\n'.join(writeMe))


# Sort two lists based on one of the lists
def qsortBasedOn(sortMe, basedOn):
    if not len(sortMe) == len(basedOn):
        return 'ERROR!'
    if len(basedOn) <= 1:
            return [sortMe, basedOn]
    pivot = basedOn.pop(0)
    pivotSM = sortMe.pop(0)
    greater = []
    lesser = []
    greaterSM = []
    lesserSM = []
    while len(basedOn) > 0:
        cur = basedOn.pop(0)
        curSM = sortMe.pop(0)
        if cur >= pivot:
            greater.append(cur)
            greaterSM.append(curSM)
        else:
            lesser.append(cur)
            lesserSM.append(curSM)
    greaterOut = qsortBasedOn(greaterSM, greater)
    lesserOut = qsortBasedOn(lesserSM, lesser)
    return [lesserOut[0] + [pivotSM] + greaterOut[0], lesserOut[1] + [pivot] + greaterOut[1]]


# Benjamini-Hochberg - takes a dictionary of { name: pValue, ... }
def benjaminiHochberg(dict1, tests, alpha=0.001):
    # First sort the results
    sorted1 = qsortBasedOn(list(dict1.keys()), list(dict1.values()))[0]
    # Then control based on FDR
    res1 = []
    alpha = float(alpha)
    #res1 = [sorted1[i] for i in range(len(sorted1)) if dict1[sorted1[i]] <= alpha/float(tests-i)]
    for i in range(len(sorted1)):
        if dict1[sorted1[i]] <= alpha*(float(i+1)/float(tests)):
            res1.append(sorted1[i])
        else:
            break
    return res1


def make_mirna_dicts(mirna_path):
    # 0. Create a dictionary to convert the miRNAs to there respective ids
    miRNAIDs = {}
    miRNAIDs_rev = {}
    with open(mirna_path, 'r') as inFile:
        while 1:
            inLine = inFile.readline()
            if not inLine:
                break
            splitUp = inLine.split(' ')
            if not splitUp[1] in miRNAIDs_rev:
                miRNAIDs_rev[splitUp[1]] = splitUp[0].lower()
            if not splitUp[0].lower() in miRNAIDs:
                miRNAIDs[splitUp[0].lower()] = splitUp[1]
            else:
                print('Uh oh! %s' % splitUp)

    return miRNAIDs, miRNAIDs_rev


def make_refseq2entrez(gene2refseq_path, gene2refseq_pkl_path):
    # 1. Read in gene2refseq mappings and make a dictionary
    if not os.path.exists(gene2refseq_pkl_path):
        with gzip.open(gene2refseq_path, 'r') as inFile:
            #inFile.readline() # skip header
            refSeq2entrez = {}
            while 1:
                line = inFile.readline()
                if not line:
                    break
                line = line.decode("utf-8")

                # Only add those that have the correct NCBI organism ID
                splitUp = line.strip().split('\t')
                if int(splitUp[0]) == 9606:
                    # Check that the nucleotide ID is not a '-' and that it has
                    # genomic coordiantes assocaited with it
                    if not splitUp[3] == '-':
                        tmp = splitUp[3].split('.')[0]
                        if not tmp in refSeq2entrez:
                            refSeq2entrez[deepcopy(tmp)] = int(splitUp[1])

        with open(gene2refseq_pkl_path, 'wb') as pklFile:
            pickle.dump(refSeq2entrez,pklFile)
    else:
        with open(gene2refseq_pkl_path, 'rb') as pklFile:
            refSeq2entrez = pickle.load(pklFile)

    print('length refseq2entrez = %d' % len(refSeq2entrez))
    return refSeq2entrez


def read_sequences(seq_path):
    # 2. Read in sequences
    with gzip.open(seq_path, 'r') as seqFile:
        seqLines = [line.decode("utf-8") for line in seqFile.readlines()]
        ids = [i.strip().split(',')[0].upper() for i in seqLines]
        sequences = [i.strip().split(',')[1] for i in seqLines]
        seqs = dict(zip(ids,sequences))
    return seqs


# these are globals used in the Weeder run
fastaFiles = None
weederPSSMs1 = None

def run_mirvestigator(seqs, refSeq2entrez, use_entrez, exp_dir='exp'):
    """
    2 Stages:
      1. run Weeder and generate the PSSMS for mirvestigator
      2. run miRvestigator
    """
    global fastaFiles, weederPSSMs1

    # Setup for multiprocessing
    mgr = Manager()
    fastaFiles = mgr.list()

    # For each cluster file in exp from Goodarzi et al.
    # Cluster files should have a header and be tab delimited to look like this:
    # Gene\tGroup\n
    # NM_000014\t52\n
    # <RefSeq_ID>\t<signature_id>\n
    # ...
    files = os.listdir(exp_dir)
    for file in files:
        # 3. Read in cluster file and convert to entrez ids
        with open(os.path.join(exp_dir, file), 'r') as inFile:
            dataset = file.strip().split('.')[0]
            inFile.readline()
            lines = inFile.readlines()
            clusters = defaultdict(list)
            for line in lines:
                gene, group = line.strip().split('\t')
                group = int(group)
                if use_entrez:
                    entrez = int(gene)
                    clusters[group].append(entrez)
                else:
                    if gene in refSeq2entrez:
                        clusters[group].append(refSeq2entrez[gene])

        # 5. Make a FASTA file & run weeder
        for cluster in clusters:
            print(cluster)
            # Get seqeunces
            clusterSeqs = {}
            for target in clusters[cluster]:
                if str(target) in seqs:
                    clusterSeqs[target] = seqs[str(target)]
                else:
                    print("Did not find seq for '%s'" % target)

            # Make FASTA file
            print('Fasta output...')
            fastaFiles.append('tmp/weeder/fasta/'+str(cluster)+'_'+str(dataset)+'.fasta')
            if not os.path.exists('tmp/weeder/fasta'):
                os.makedirs('tmp/weeder/fasta')

            with open('tmp/weeder/fasta/'+str(cluster)+'_'+str(dataset)+'.fasta','w') as fastaFile:
                for seq in clusterSeqs:
                    fastaFile.write('>'+str(seq)+'\n'+str(clusterSeqs[seq])+'\n')

    # Run this using all cores available
    weederPSSMs1 = mgr.list()
    print('Starting Weeder runs...')
    cpus = cpu_count()
    print('There are %d CPUs available.' % cpus)
    pool = Pool(processes=cpus)
    pool.map(runWeeder, range(len(fastaFiles)))
    print('Done with Weeder runs.')

    # Compare to miRDB using my program
    return [{'name': p.name,
             'nsites': p.nsites,
             'evalue': p.eValue,
             'genes': p.genes,
             'matrix': p.matrix}
            for p in list(weederPSSMs1)]

    # run mirvestigator here
    """
    m2m = miRvestigator(weederPSSMs1, seqs.values(), seedModel=[6,7,8], minor=True, p5=True,
                        p3=True, wobble=False, wobbleCut=0.25)
    with open('m2m' + '_' + str(dataset) + '.pkl', 'wb') as outFile:
        pickle.dump(m2m, outFile)
    """


db = None
dataset = None
datasetGenes = None
totalTargets = None
clusters = None
miRNATargetDict = None

def run_target_prediction_dbs(refSeq2entrez, use_entrez, exp_dir='exp',
                              pred_db_dir='TargetPredictionDatabases'):
    global db, dataset, datasetGenes, totalTargets, clusters, miRNATargetDict

    dbs = [d for d in os.listdir(pred_db_dir)
           if os.path.isdir(os.path.join(pred_db_dir, d))]
    # Now do PITA and TargetScan - iterate through both platforms
    for db in dbs:
        # Get ready for multiprocessor goodness
        mgr = Manager()
        cpus = cpu_count()

        # Load up db of miRNA ids
        ls2 = [x for x in os.listdir(os.path.join(pred_db_dir, db))
               if x.endswith('.csv')]

        # Load the predicted target genes for each miRNA from the files
        tmpDict = {}
        for f in ls2:
            miRNA = f.rstrip('.csv')
            with open(os.path.join(pred_db_dir, db, f), 'r') as inFile:
                tmpDict[miRNA.lower()] = [int(line.strip()) for line in inFile.readlines()
                                          if line.strip()]
        miRNATargetDict = mgr.dict(tmpDict)

        # Total background
        print('\n2')
        with open(os.path.join(pred_db_dir, db, db + '_ids_entrez.bkg'), 'r') as inFile:
            targetList = [int(x) for x in inFile.readlines() if x]
            tmp1 = targetList
            totalTargets = mgr.list(tmp1)

        # For each cluster file in expfiles from Goodarzi et al.
        files = os.listdir(exp_dir)
        for file in files:
            # 3. Read in cluster file and convert to entrez ids
            with open(os.path.join(exp_dir, file), 'r') as inFile:
                dataset = mgr.list([file.strip().split('.')[0]])
                print(dataset[0])
                inFile.readline()
                lines = inFile.readlines()
                tmpDict = defaultdict(list)
                genes = []
                for line in lines:
                    gene, group = line.strip().split('\t')
                    group = int(group)
                    if use_entrez:
                        entrez = int(gene)
                    else:
                        if gene in refSeq2entrez:
                            entrez = refSeq2entrez[gene]
                        else:
                            entrez = None

                    if entrez in targetList:
                        genes.append(entrez)
                        tmpDict[group].append(entrez)

            clusters = mgr.dict(tmpDict)
            datasetGenes = mgr.list(genes)

            # Iterate through clusters and compute p-value for each miRNA
            if not os.path.exists('miRNA_' + db):
                os.mkdir('miRNA_' + db)

            # Run this using all cores available
            print('Starting '+ dataset[0] + ' runs...')
            keys2 = clusters.keys()
            pool = Pool(processes=cpus)
            pool.map(clusterHypergeo,keys2)
            print('Done.')

        # 1. Get a list of all files in miRNA directory
        overlapFiles = os.listdir('miRNA_' + db)

        # 2. Read them all in and grab the top hits
        with open('miRNA/mergedResults_' + db + '.csv','w') as outFile:
            outFile.write('Dataset,Cluster,miRNA,q,m,n,k,p.value')
            enrichment = []
            for overlapFile in overlapFiles:
                inFile = open('miRNA_' + db + '/' + overlapFile, 'r')
                inFile.readline() # Get rid of header
                lines = [line.strip().split(',') for line in inFile.readlines()]
                miRNAs = [line[0].lstrip(db+'_') for line in lines]
                intSect = [line[1] for line in lines]
                miRNAPred = [line[2] for line in lines]
                allNum = [line[3] for line in lines]
                clustGenes = [line[4] for line in lines]
                pVals = [float(line[5]) for line in lines]
                inFile.close()
                min1 = float(1)
                curMiRNA = []
                daRest = []
                for i in range(len(miRNAs)):
                    if pVals[i] < min1 and int(intSect[i])>=1:
                        min1 = pVals[i]
                        tmpMiRNA = miRNAs[i].lower()
                        if tmpMiRNA[-3:]=='-5p':
                            tmpMiRNA = tmpMiRNA[:-3]
                        curMiRNA = [tmpMiRNA]
                        daRest = [intSect[i], miRNAPred[i], allNum[i], clustGenes[i]]
                    elif pVals[i]==min1 and int(intSect[i])>=1:
                        tmpMiRNA = miRNAs[i].lower()
                        if tmpMiRNA[-3:]=='-5p':
                            tmpMiRNA = tmpMiRNA[:-3]
                        curMiRNA.append(tmpMiRNA)
                tmp = overlapFile.rstrip('.csv').split('_')
                dataset = tmp[0]+'_'+tmp[1]+'_'+tmp[2]
                cluster = tmp[3]
                outFile.write('\n' + dataset + ',' + cluster + ',' + ' '.join(curMiRNA) +
                              ',' + ','.join(daRest) + ',' + str(min1))
                enrichment.append({'dataset':dataset,
                                   'cluster':cluster,
                                   'miRNA':curMiRNA,
                                   'q' : daRest[0],
                                   'm' : daRest[1],
                                   'n' : daRest[2],
                                   'k' : daRest[3],
                                   'pValue' : min1,
                                   'percTargets' : float(daRest[0]) / float(daRest[3]),
                                   'significant' : False})

        # Filter using benjamini-hochberg FDR <= 0.001, >=10% target genes in cluster
        bhDict = {}
        for clust in range(len(enrichment)):
            bhDict[enrichment[clust]['dataset']+'_'+enrichment[clust]['cluster']] = enrichment[clust]['pValue']
        significant = benjaminiHochberg(bhDict, tests=len(clusters), alpha=0.001)
        # Do filtering
        filtered = []
        for clust in range(len(enrichment)):
            if (enrichment[clust]['dataset']+'_'+enrichment[clust]['cluster'] in significant) and (float(enrichment[clust]['q'])/float(enrichment[clust]['k']) >= 0.1):
                enrichment[clust]['significant'] = True
                filtered.append(enrichment[clust])

        # Write out filtered results
        with open('filtered_' + db + '.csv','w') as outFile:
            outFile.write('Dataset,Signature,miRNA,Percent.Targets')
            tot = 0
            for clust in range(len(filtered)):
                outFile.write('\n'+filtered[clust]['dataset']+','+filtered[clust]['cluster']+','+miRNA+','+str(float(enrichment[clust]['q'])/float(enrichment[clust]['k'])))


def write_combined_report():
    # Get miRvestigator results
    miRNA_matches = {}
    with open(os.path.join('miRNA', 'scores.csv'),'r') as inFile:
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]
        for line in lines:
            if not line[1]=='NA':
                miRNA_mature_seq_ids = []
                for i in line[1].split('_'):
                    miRNA_mature_seq_ids += miRNAInDict(i.lower(),miRNAIDs)
                cluster_name = [i for i in line[0].split('_')]
                cluster_name = cluster_name[1]+'_'+cluster_name[2]+'_'+cluster_name[3]+'_'+cluster_name[0]
                miRNA_matches[cluster_name] = {'miRNA':line[1],'model':line[2],'mature_seq_ids':miRNA_mature_seq_ids}

    print('Loaded miRvestigator.')
    # Get PITA results
    with open(os.path.join('miRNA', 'mergedResults_PITA.csv'), 'r') as inFile:
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]

    pita_miRNA_matches = {}
    for line in lines:
        if not line[2]=='':
            miRNA_mature_seq_ids = []
            mirs = [i.lower().strip('pita_') for i in line[2].split(' ')]
            for i in mirs:
                miRNA_mature_seq_ids += miRNAInDict(i,miRNAIDs)
            if not line[0]+'_'+line[1] in miRNA_matches:
                miRNA_matches[line[0]+'_'+line[1]] = {'pita_miRNA':' '.join(mirs),'pita_perc_targets':str(float(line[3])/float(line[6])),'pita_pValue':line[7],'pita_mature_seq_ids':miRNA_mature_seq_ids}
            else:
                miRNA_matches[line[0]+'_'+line[1]]['pita_miRNA'] = ' '.join(mirs)
                miRNA_matches[line[0]+'_'+line[1]]['pita_perc_targets'] = str(float(line[3])/float(line[6]))
                miRNA_matches[line[0]+'_'+line[1]]['pita_pValue'] = line[7]
                miRNA_matches[line[0]+'_'+line[1]]['pita_mature_seq_ids'] = miRNA_mature_seq_ids
    print('Loaded PITA.')

    # Get TargetScan results
    with open(os.path.join('miRNA', 'mergedResults_TargetScan.csv'),'r') as inFile:
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]

    targetScan_miRNA_matches = {}
    for line in lines:
        if not line[2]=='':
            miRNA_mature_seq_ids = []
            mirs = [i.lower().strip('scan_') for i in line[2].split(' ')]
            for i in mirs:
                miRNA_mature_seq_ids += miRNAInDict(i.lower().strip('targetscan_'),miRNAIDs)
            if not line[0]+'_'+line[1] in miRNA_matches:
                miRNA_matches[line[0]+'_'+line[1]] = {'ts_miRNA':' '.join(mirs),'ts_perc_targets':str(float(line[3])/float(line[6])),'ts_pValue':line[7],'ts_mature_seq_ids':miRNA_mature_seq_ids}
            else:
                miRNA_matches[line[0]+'_'+line[1]]['ts_miRNA'] = ' '.join(mirs)
                miRNA_matches[line[0]+'_'+line[1]]['ts_perc_targets'] = str(float(line[3])/float(line[6]))
                miRNA_matches[line[0]+'_'+line[1]]['ts_pValue'] = line[7]
                miRNA_matches[line[0]+'_'+line[1]]['ts_mature_seq_ids'] = miRNA_mature_seq_ids
    print('Loaded TargetScan.')

    # Big list of all miRNAs for all clusters
    with open('combinedResults.csv','w') as outFile:
        outFile.write('Dataset,signature,miRvestigator.miRNA,miRvestigator.model,miRvestigator.mature_seq_ids,PITA.miRNA,PITA.percent_targets,PITA.P_Value,PITA.mature_seq_ids,TargetScan.miRNA,TargetScan.percent_targets,TargetScan.P_Value,TargetScan.mature_seq_ids')
        for i in miRNA_matches:
            splitUp = i.split('_')
            writeMe = '\n'+splitUp[0]+'_'+splitUp[1]+'_'+splitUp[2]+','+splitUp[3]
            if 'miRNA' in miRNA_matches[i]:
                writeMe += ',' + miRNA_matches[i]['miRNA'] + ','+miRNA_matches[i]['model']+',' + ' '.join(miRNA_matches[i]['mature_seq_ids'])
            else:
                writeMe += ',NA,NA,NA'
            if 'pita_miRNA' in miRNA_matches[i]:
                writeMe += ',' + miRNA_matches[i]['pita_miRNA'] + ',' + miRNA_matches[i]['pita_perc_targets'] + ',' + miRNA_matches[i]['pita_pValue'] + ',' + ' '.join(miRNA_matches[i]['pita_mature_seq_ids'])
            else:
                writeMe += ',NA,NA,NA,NA'
            if 'ts_miRNA' in miRNA_matches[i]:
                writeMe += ',' + miRNA_matches[i]['ts_miRNA'] + ',' + miRNA_matches[i]['ts_perc_targets'] + ',' + miRNA_matches[i]['ts_pValue'] + ',' + ' '.join(miRNA_matches[i]['ts_mature_seq_ids'])
            else:
                writeMe += ',NA,NA,NA,NA'
            outFile.write(writeMe)


DESCRIPTION = """FIRM - Framework for Inference of Regulation by miRNAs"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('-ue', '--use_entrez', action='store_true',
                        help="input file uses entrez IDs instead of RefSeq")
    args = parser.parse_args()

    use_entrez = args.use_entrez
    miRNAIDs, miRNAIDS_rev = make_mirna_dicts('common/hsa.mature.fa')
    refSeq2entrez = make_refseq2entrez('common/gene2refseq.gz', 'common/refSeq2entrez.pkl')
    seqs = read_sequences('common/p3utrSeqs_Homo_sapiens.csv.gz')

    # this is actually running weeder
    weeder_pssms = run_mirvestigator(seqs, refSeq2entrez, use_entrez)
    with open('pssms.json', 'w') as outfile:
        json.dump(weeder_pssms, outfile)

    with open('seqs.txt', 'w') as outfile:
        for seq in seqs.values():
            outfile.write("%s\n", seq)

    # TODO: run our actual weeder as external tool


    """
    run_target_prediction_dbs(refSeq2entrez, use_entrez)
    write_combined_report()
    """
