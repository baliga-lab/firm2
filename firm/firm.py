import argparse
from firm.pssm import pssm
import gzip, os, re, sys
from subprocess import *
from multiprocessing import Pool, cpu_count, Manager
from collections import defaultdict
import json


def mirnas_for(mirna, mirna_ids):
    return [mirnas for m_name, mirnas in mirna_ids.items()
            if mirna_name_matches(mirna, m_name)]


def mirna_name_matches(a, b):
    if a == b:
        return 1
    if len(a) < len(b):
        re1 = re.compile(a + '[a-z]$')
        if re1.match(b):
            return 1
    else:
        re1 = re.compile(b + '[a-z]$')
        if re1.match(a):
            return 1
    return 0


def check_weederlauncher():
    """Check whether weederlauncher exists on this system"""
    for d in os.environ["PATH"].split(os.pathsep):
        filepath = os.path.join(d, 'weederlauncher')
        if os.path.isfile(filepath) and os.access(filepath, os.X_OK):
            return True
    return False


# Run weeder and parse its output
# First weederTFBS -W 6 -e 1, then weederTFBS -W 8 -e 2, and finally adviser
def run_weeder(params):
    seqFile, tmpdir, freqfiles = params
    print(".", end="", file=sys.stderr, flush=True)
    weeder_pssms = []
    percTargets = 50
    revComp = False

    # First run weederTFBS for 6bp motifs
    weeder_args = ' ' + str(seqFile) + ' HS3P small T50 F%s' % freqfiles
    if revComp:
        weeder_args += ' -S'
    errOut = open(os.path.join(tmpdir, 'weeder', 'stderr.out'), 'w')
    weederProc = Popen("weederlauncher " + weeder_args, shell=True, stdout=PIPE, stderr=errOut)
    output = weederProc.communicate()

    # Now parse output from weeder
    PSSMs = []
    with open(str(seqFile)+'.wee','r') as output:
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
            matrix += [[ float(nums[0]) / float(colSum),
                         float(nums[1]) / float(colSum),
                         float(nums[2]) / float(colSum),
                         float(nums[3]) / float(colSum)]]
            col = outLines.pop(0)
        weeder_pssms.append(pssm(name=name,
                                 sites=instances,
                                 evalue=hitBp[len(matrix)][1],
                                 pssm=matrix,
                                 genes=redMotifs))
    return weeder_pssms

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


def cluster_hypergeo(params):
    global db, dataset_genes, total_targets, mirna_target_dict
    cluster, cluster_genes, dataset, dboutdir = params
    sys.stderr.write(".")
    sys.stderr.flush()
    outpath = os.path.join(dboutdir, '%s_%s.csv' % (dataset, str(cluster)))
    with open(outpath, 'w') as outfile:
        outfile.write('miRNA,Cluster.Targets,miRNA.Targets,Cluster.Genes,Total,P.Value\n')
        # k = overlap, N = potential target genes, n = miRNA targets, m = cluster genes
        # Take gene list and compute overlap with each miRNA
        allGenes = set(dataset_genes).intersection(set(total_targets))
        genes = set(cluster_genes).intersection(set(allGenes))
        writeMe = []
        keys1 = mirna_target_dict.keys()
        m1s = []
        q = []
        m = []
        n = []
        k = []
        for m1 in keys1:
            m1s.append(m1)
            miRNAGenes = set(mirna_target_dict[m1]).intersection(allGenes)
            q.append(len(set(miRNAGenes).intersection(genes)))
            m.append(len(miRNAGenes))
            n.append(len(allGenes) - len(miRNAGenes))
            k.append(len(genes))
        results = phyper(q, m, n, k)
        for i in range(len(m1s)):
            writeMe.append(str(m1s[i]) + ',' + str(q[i]) + ',' + str(m[i]) + ',' + str(n[i]) + ',' + str(k[i]) + ',' + str(results[i]))
        outfile.write('\n'.join(writeMe))


# Benjamini-Hochberg - takes a dictionary of { name: pValue, ... }
def benjamini_hochberg(dict1, tests, alpha=0.001):
    # create a list of names sorted by the p-values in the input dictionary
    sorted1 = [i[0] for i in sorted(list(dict1.items()), key=lambda t: t[1])]

    # Then control based on FDR (false discovery rate)
    res1 = []
    alpha = float(alpha)
    for i in range(len(sorted1)):
        if dict1[sorted1[i]] <= alpha*(float(i+1)/float(tests)):
            res1.append(sorted1[i])
        else:
            break
    return res1


def make_mirna_dict(mirna_path):
    """reads the mature.fa.gz and returns a mapping from miRNA names => miRNA ids.
    This was changed from the original FIRM so we can use updated mature.fa.gz
    files instead of a prefiltered, potentially outdated file"""
    mirna_ids = {}
    with gzip.open(mirna_path, 'r') as infile:
        while True:
            desc = infile.readline()
            if not desc:
                break
            infile.readline()  # skip the sequence line
            desc = desc.decode('utf-8')
            if desc.startswith('>hsa'):
                comps = desc.lstrip('>').split(' ')
                mirna_name = comps[0].lower()
                mirna_id = comps[1]

                if not mirna_name in mirna_ids:
                    mirna_ids[mirna_name] = mirna_id
                else:
                    raise Exception("already exists: '%s'" % mirna_name)

    return mirna_ids


def make_refseq2entrez(gene2refseq_path):
    # Read in gene2refseq mappings and make a dictionary
    with gzip.open(gene2refseq_path, 'r') as inFile:
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
                        refSeq2entrez[tmp] = int(splitUp[1])

        return refSeq2entrez


def read_sequences(seq_path):
    # 2. Read in sequences
    with gzip.open(seq_path, 'r') as seqFile:
        seqLines = [line.decode("utf-8") for line in seqFile.readlines()]
        ids = [i.strip().split(',')[0].upper() for i in seqLines]
        sequences = [i.strip().split(',')[1] for i in seqLines]
        seqs = dict(zip(ids,sequences))
    return seqs


def prepare_weeder_input(seqs, refSeq2entrez, use_entrez, exp_dir, tmpdir):
    # For each cluster file in exp from Goodarzi et al.
    # Cluster files should have a header and be tab delimited to look like this:
    # Gene\tGroup\n
    # NM_000014\t52\n
    # <RefSeq_ID>\t<signature_id>\n
    # ...
    fasta_files = []
    files = os.listdir(exp_dir)
    weeder_fasta_dir = os.path.join(tmpdir, 'weeder', 'fasta')
    if not os.path.exists(weeder_fasta_dir):
        os.makedirs(weeder_fasta_dir)

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
            # Get sequences
            cluster_seqs = {}
            for target in clusters[cluster]:
                if str(target) in seqs:
                    cluster_seqs[target] = seqs[str(target)]
                else:
                    print("Did not find seq for '%s' (cluster %d)" % (target, cluster),
                          file=sys.stderr)

            # Make FASTA file
            fname = "%d_%s.fasta" % (cluster, dataset)
            fpath = os.path.join(weeder_fasta_dir, fname)
            fasta_files.append(fpath)

            with open(fpath, 'w') as outfile:
                for seq_name, seq in cluster_seqs.items():
                    outfile.write('>%s\n' % seq_name)
                    outfile.write('%s\n' % seq)
    return fasta_files


def find_motifs(fasta_files, tmpdir, freqfiles):
    # Setup for multiprocessing
    # Run this using all cores available
    cpus = cpu_count()
    print('Starting Weeder runs (%d CPUs available)...' % cpus, file=sys.stderr)
    pool = Pool(processes=cpus)
    pssms_list = pool.map(run_weeder, [(f, tmpdir, freqfiles) for f in fasta_files])
    print('Done with Weeder runs.', file=sys.stderr)

    # Compare to miRDB using my program
    final_pssms = [pssm for pssms in pssms_list for pssm in pssms]
    return [{'name': p.name,
             'sites': p.sites,
             'evalue': p.evalue,
             'genes': p.genes,
             'matrix': p.matrix}
            for p in final_pssms]


db = None
dataset = None
dataset_genes = None
total_targets = None
mirna_target_dict = None

def run_target_prediction_dbs(refSeq2entrez, exp_dir,
                              outdir, tmpdir, pred_db_dir, use_entrez):
    global db, dataset, dataset_genes, total_targets, mirna_target_dict

    mgr = Manager()
    dbs = [d for d in os.listdir(pred_db_dir)
           if os.path.isdir(os.path.join(pred_db_dir, d))]
    # Now do PITA and TargetScan - iterate through both platforms
    for db in dbs:
        print("checking against prediction database: '%s'" % db, file=sys.stderr)
        dboutdir = os.path.join(tmpdir, 'miRNA_%s' % db)
        if not os.path.exists(dboutdir):
            os.mkdir(dboutdir)


        # Get ready for multiprocessor goodness
        cpus = cpu_count()

        # Load up db of miRNA ids
        ls2 = [x for x in os.listdir(os.path.join(pred_db_dir, db))
               if x.endswith('.csv')]

        # Load the predicted target genes for each miRNA from the files
        tmp_dict = {}
        for f in ls2:
            miRNA = f.rstrip('.csv')
            with open(os.path.join(pred_db_dir, db, f), 'r') as inFile:
                tmp_dict[miRNA.lower()] = [int(line.strip()) for line in inFile.readlines()
                                          if line.strip()]
        mirna_target_dict = mgr.dict(tmp_dict)

        # Total background
        with open(os.path.join(pred_db_dir, db, db + '_ids_entrez.bkg'), 'r') as inFile:
            target_list = [int(x) for x in inFile.readlines() if x]
            tmp1 = target_list
            total_targets = mgr.list(tmp1)

        # For each cluster file in expfiles from Goodarzi et al.
        files = os.listdir(exp_dir)
        for f in files:
            # 3. Read in cluster file and convert to entrez ids
            print("File '%s'" % f)
            with open(os.path.join(exp_dir, f), 'r') as inFile:
                dataset = f.strip().split('.')[0]
                print("Data set: '%s'..." % dataset, file=sys.stderr)
                inFile.readline()
                lines = inFile.readlines()
                clusters = defaultdict(list)
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

                    if entrez in target_list:
                        genes.append(entrez)
                        clusters[group].append(entrez)

            dataset_genes = mgr.list(genes)

            # Iterate through clusters and compute p-value for each miRNA
            # Run this using all cores available
            pool = Pool(processes=cpus)
            params = [(cluster, cluster_genes, dataset, dboutdir)
                      for cluster, cluster_genes in clusters.items()]
            pool.map(cluster_hypergeo, params)
            print('Done.', file=sys.stderr)

        # 1. Get a list of all resulting overlap files in miRNA directory
        overlapFiles = os.listdir(dboutdir)

        # 2. Read them all in and grab the top hits
        with open(os.path.join(outdir, 'mergedResults_%s.csv' % db), 'w') as outFile:
            outFile.write('Dataset,Cluster,miRNA,q,m,n,k,p.value')
            enrichment = []
            for overlapFile in overlapFiles:
                overlap_path = os.path.join(dboutdir, overlapFile)
                with open(overlap_path, 'r') as inFile:
                    inFile.readline() # Get rid of header
                    lines = [line.strip().split(',') for line in inFile.readlines()]
                    miRNAs = [line[0].lstrip(db+'_') for line in lines]
                    intSect = [line[1] for line in lines]
                    miRNAPred = [line[2] for line in lines]
                    allNum = [line[3] for line in lines]
                    clustGenes = [line[4] for line in lines]
                    pVals = [float(line[5]) for line in lines]

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
                comps = overlapFile.rstrip('.csv').split('_')
                dataset = '_'.join(comps[:-1])
                cluster = comps[-1]
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
        significant = benjamini_hochberg(bhDict, tests=len(clusters), alpha=0.001)
        # Do filtering
        filtered = []
        for clust in range(len(enrichment)):
            if (enrichment[clust]['dataset']+'_'+enrichment[clust]['cluster'] in significant) and (float(enrichment[clust]['q'])/float(enrichment[clust]['k']) >= 0.1):
                enrichment[clust]['significant'] = True
                filtered.append(enrichment[clust])

        # Write out filtered results
        fpath = os.path.join(outdir, 'filtered_%s.csv' % db)
        with open(fpath, 'w') as outFile:
            outFile.write('Dataset,Signature,miRNA,Percent.Targets')
            tot = 0
            for clust in range(len(filtered)):
                outFile.write('\n'+filtered[clust]['dataset']+','+filtered[clust]['cluster']+','+miRNA+','+str(float(enrichment[clust]['q'])/float(enrichment[clust]['k'])))


def write_combined_report(mirv_score_path, mirbase_path, outdir):
    mirna_ids = make_mirna_dict(mirbase_path)
    # Get miRvestigator results
    print("Retrieving miRvestigator results...", file=sys.stderr, end="", flush=True)
    miRNA_matches = {}
    with open(mirv_score_path,'r') as inFile:
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]
        for line in lines:
            if not line[1]=='NA':
                miRNA_mature_seq_ids = []
                for i in line[1].split('_'):
                    miRNA_mature_seq_ids += mirnas_for(i.lower(), mirna_ids)
                cluster_name = [i for i in line[0].split('_')]
                cluster_name = cluster_name[1]+'_'+cluster_name[2]+'_'+cluster_name[3]+'_'+cluster_name[0]
                miRNA_matches[cluster_name] = {'miRNA':line[1],'model':line[2],'mature_seq_ids':miRNA_mature_seq_ids}

    print('Done.', file=sys.stderr)

    # Get PITA results
    print("Retrieving PITA results...", file=sys.stderr, end="", flush=True)
    with open(os.path.join(outdir, 'mergedResults_PITA.csv'), 'r') as inFile:
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]

    for line in lines:
        if not line[2]=='':
            miRNA_mature_seq_ids = []
            mirs = [i.lower().strip('pita_') for i in line[2].split(' ')]
            for i in mirs:
                miRNA_mature_seq_ids += mirnas_for(i, mirna_ids)
            if not line[0]+'_'+line[1] in miRNA_matches:
                miRNA_matches[line[0]+'_'+line[1]] = {'pita_miRNA':' '.join(mirs),'pita_perc_targets':str(float(line[3])/float(line[6])),'pita_pValue':line[7],'pita_mature_seq_ids':miRNA_mature_seq_ids}
            else:
                miRNA_matches[line[0]+'_'+line[1]]['pita_miRNA'] = ' '.join(mirs)
                miRNA_matches[line[0]+'_'+line[1]]['pita_perc_targets'] = str(float(line[3])/float(line[6]))
                miRNA_matches[line[0]+'_'+line[1]]['pita_pValue'] = line[7]
                miRNA_matches[line[0]+'_'+line[1]]['pita_mature_seq_ids'] = miRNA_mature_seq_ids
    print('Done.', file=sys.stderr)

    # Get TargetScan results
    print("Retrieving Targetscan results...", file=sys.stderr, end="", flush=True)
    with open(os.path.join(outdir, 'mergedResults_TargetScan.csv'),'r') as inFile:
        inFile.readline() # get rid of header
        lines = [i.strip().split(',') for i in inFile.readlines()]

    for line in lines:
        if not line[2]=='':
            miRNA_mature_seq_ids = []
            mirs = [i.lower().strip('scan_') for i in line[2].split(' ')]
            for i in mirs:
                miRNA_mature_seq_ids += mirnas_for(i.lower().strip('targetscan_'), mirna_ids)
            if not line[0]+'_'+line[1] in miRNA_matches:
                miRNA_matches[line[0]+'_'+line[1]] = {'ts_miRNA':' '.join(mirs),'ts_perc_targets':str(float(line[3])/float(line[6])),'ts_pValue':line[7],'ts_mature_seq_ids':miRNA_mature_seq_ids}
            else:
                miRNA_matches[line[0]+'_'+line[1]]['ts_miRNA'] = ' '.join(mirs)
                miRNA_matches[line[0]+'_'+line[1]]['ts_perc_targets'] = str(float(line[3])/float(line[6]))
                miRNA_matches[line[0]+'_'+line[1]]['ts_pValue'] = line[7]
                miRNA_matches[line[0]+'_'+line[1]]['ts_mature_seq_ids'] = miRNA_mature_seq_ids
    print('Done.', file=sys.stderr)

    # Big list of all miRNAs for all clusters
    with open(os.path.join(outdir, 'combinedResults.csv'), 'w') as outFile:
        outFile.write('Dataset,signature,miRvestigator.miRNA,miRvestigator.model,miRvestigator.mature_seq_ids,PITA.miRNA,PITA.percent_targets,PITA.P_Value,PITA.mature_seq_ids,TargetScan.miRNA,TargetScan.percent_targets,TargetScan.P_Value,TargetScan.mature_seq_ids')
        for i in miRNA_matches:
            comps = i.split('_')
            dataset = '_'.join(comps[:-1])
            cluster = comps[-1]
            writeMe = '\n%s,%s' % (dataset, cluster)
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
