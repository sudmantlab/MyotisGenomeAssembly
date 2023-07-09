         seed = -1
      seqfile = data/genes/allbats_fasta/alignments_concat_all/allBatGenes.phy
     treefile = data/allBats.GHOST.mcmcInputTree
      outfile = output/mcmctree/allBats_codon/out

        ndata = 21693
        ngene = 21693
      seqtype = 1
      usedata = 2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
        clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates

        model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        alpha = 0.5  * alpha for gamma rates at sites
        ncatG = 5    * No. categories in discrete gamma

    cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

        print = 1
*       burnin = 10000
*     sampfreq = 100
*      nsample = 10000

