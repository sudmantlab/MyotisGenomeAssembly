         seed = -1
      seqfile = alignments_noATGPlaceholder_batsOnly_universal.c3.phy
     treefile = 536Mammals-ChiropteraOnly-rooted-noBL-noSturnira_hondurensis-times.nwk
     mcmcfile = Runs_real_combined.mcmc
      outfile = Runs_real_combined.out

        ndata = 1
      seqtype = 0
      usedata = 2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
        clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates
      RootAge = '<0.66'

        model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        alpha = 0.5  * alpha for gamma rates at sites
        ncatG = 5    * No. categories in discrete gamma  

    cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?
      BDparas = 1 1 0.1

  rgene_gamma = 2 40 1
 sigma2_gamma = 1 10 1
        print = -1

