import sys
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import multiprocessing as mp
import random
import computePercentIdentityParallel as cpi


results = []

outname = str(sys.argv[1])
#reftax = str(sys.argv[2])
reffasta = str(sys.argv[2])
seqtax = str(sys.argv[3])
seqfasta = str(sys.argv[4])
cpus = int(sys.argv[5])
#seqdf = pd.read_csv(str(sys.argv[3]))

#rftx = readTaxonomy(reftax)
rffa = cpi.readFasta(reffasta)
rffa.index = rffa.id
sqtx = cpi.readTaxonomy(seqtax)
sqfa = cpi.readFasta(seqfasta)

#refdf = rftx.merge(rffa, how='inner', on='id')
refdf = rffa
seqdf = sqtx.merge(sqfa, how='inner', on='id')

seqdf['taxid'] = seqdf.taxonomy.apply(lambda x: x.split("i__")[-1][:-1])
unknown= "unknown;unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);unknown_unclassified(0);"
knowndf = seqdf[seqdf.taxonomy != unknown]

pool = mp.Pool(cpus)

def collect_results(result):
    global results
    results.append(result)

print("Partitioning sequences")
print(len(knowndf))
dfdivisions = int(np.round(len(knowndf)/cpus))
indices = [x*dfdivisions for x in list(range(cpus))]

dfs = [knowndf.iloc[x:x+dfdivisions,] for x in indices[:-1]]
dfs.append(knowndf.iloc[indices[-1]:,])
print(dfs[1].shape)

print("Initializing pool")
for i, df in enumerate(dfs):
    pool.apply_async(cpi.calcPercDF, args=(i, df, refdf), callback=collect_results)
    print("Adding partition " + str(i))

print("Closing Pool")
pool.close()
print("Running")
pool.join()

#getGenusDiv(x, spdfc, 50) for x in genera[:10]
#[getGenusDiv(x, spdfc, 50) for x in genera[:10]]

print("Sorting results")
results.sort(key=lambda x: x[0])
finaldfs = [r for i,r in results]

bigdf = pd.concat(finaldfs)

print("Finished")
#knowndf['percentIdentity'] = knowndf.apply(lambda x: getAltPerc(x, refdf), axis=1)

print("completed percent identity computation.")
bigdf.to_csv(outname, index=False)