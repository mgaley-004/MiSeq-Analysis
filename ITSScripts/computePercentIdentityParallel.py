import sys
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import multiprocessing as mp
import random


results = []
 
 
def percIdentity(analign): 
    startalign = analign[-2] 
    endalign = analign[-1] 
    seq1 = analign[0][startalign:endalign] 
    seq2 = analign[1][startalign:endalign] 
    matches = sum(np.char.equal(list(seq1), list(seq2))) 
    gaps = (seq1+seq2).count('-') 
    score = matches/(len(seq1)-gaps) 
    return(score)
    
def hammingDistance(analign): 
    startalign = analign[-2] 
    endalign = analign[-1] 
    seq1 = analign[0][startalign:endalign] 
    seq2 = analign[1][startalign:endalign] 
    matches = sum(np.char.equal(list(seq1), list(seq2))) 
    score = matches/(len(seq1)) 
    return(score) 
     
def maxIdentity(alignments): 
    scores = [percIdentity(a) for a in alignments] 
    return(max(scores)) 
     
def getPerc(arecord, rdf): 
    tax = arecord.taxonomy 
    fasta = arecord.fasta 
    refseqs = rdf[rdf.taxonomy==tax].fasta
    if len(refseqs) > 3:
        refseqs = refseqs[:3]
    alignments = [pairwise2.align.localms(fasta, x, 2, -2, -2, -1) for x in refseqs]
    scores = [maxIdentity(a) for a in alignments]
    if len(scores) > 0:
        maxscore = max(scores)
        print(maxscore)
        return(maxscore)
    else:
        print("bad alignment")
        return(np.nan)
        
def getAltPerc(arecord, rdf): 
    taxid = arecord.taxid 
    fasta = arecord.fasta 
    refseq = rdf.loc[taxid,'fasta'] 
    alignments = pairwise2.align.localms(fasta, refseq, 2, -2, -2, -1)
    bestalign = alignments[np.argmax(np.array(alignments)[:,4])]
    score = percIdentity(bestalign)
    print(score)
    return(score)
    
def readTaxonomy(filename):
    f = open(filename)
    flines = f.read()
    f.close()
    frecords = flines.split("\n")[:-1]
    fids = [x.split("\t")[0] for x in frecords]
    ftax = [x.split("\t")[1] for x in frecords]
    return(pd.DataFrame({"id":fids, "taxonomy":ftax}))
    
def readFasta(filename):
    f = open(filename)
    flines = f.read()
    f.close()
    frecords = flines.split(">")[1:]
    fids = [x.split("\n")[0] for x in frecords]
    ffasta = [x.split("\n")[1] for x in frecords]
    return(pd.DataFrame({"id":fids, "fasta":ffasta}))
    
def calcPercDF(i, df, refdf):
    df['percentIdentity'] = df.apply(lambda x: getAltPerc(x, refdf), axis=1)
    return(i, df)
    
def collect_results(result):
    global results
    results.append(result)
    
if __name__ == "__main__":
    outname = str(sys.argv[1])
    #reftax = str(sys.argv[2])
    reffasta = str(sys.argv[2])
    seqtax = str(sys.argv[3])
    seqfasta = str(sys.argv[4])
    #seqdf = pd.read_csv(str(sys.argv[3]))
    
    #rftx = readTaxonomy(reftax)
    rffa = readFasta(reffasta)
    rffa.index = rffa.id
    sqtx = readTaxonomy(seqtax)
    sqfa = readFasta(seqfasta)
    
    #refdf = rftx.merge(rffa, how='inner', on='id')
    refdf = rffa
    seqdf = sqtx.merge(sqfa, how='inner', on='id')
    
    seqdf['taxid'] = seqdf.taxonomy.apply(lambda x: x.split("i__")[-1][:-1])
    unknown= "unknown;unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);"
    knowndf = seqdf[seqdf.taxonomy != unknown]
    
    cpus = mp.cpu_count()
    pool = mp.Pool(cpus)
    
    dfdivisions = int(np.round(len(knowndf)/cpus))
    indices = [x*dfdivisions for x in list(range(cpus))]
    
    dfs = [knowndf.iloc[x:x+dfdivisions,] for x in indices[:-1]]
    dfs.append(knowndf.iloc[indices[-1]:,])
    
    for i, df in enumerate(dfs):
        pool.apply_async(calcPercDF, args=(i, df, refdf), callback=collect_results)
        
    pool.close()
    pool.join()
    
    #getGenusDiv(x, spdfc, 50) for x in genera[:10]
    #[getGenusDiv(x, spdfc, 50) for x in genera[:10]]
    
    results.sort(key=lambda x: x[0])
    finaldfs = [r for i,r in results]
    
    bigdf = pd.concat(finaldfs)
    
    #knowndf['percentIdentity'] = knowndf.apply(lambda x: getAltPerc(x, refdf), axis=1)
    
    print("completed percent identity computation.")
    bigdf.to_csv(outname, index=False)