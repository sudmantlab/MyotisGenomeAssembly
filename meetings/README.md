


## 2023/03/13

- David will finish the annotations by the end of the week
    - we repair these w/ David's pipeline
    - these will be the final annotations
- There are several holes which are annotated in all but one genome. How to fix?
    - mask the annotations that exist, and blast them against the genomes to fill holes, you find missing genes that are homologs
    - David will run this
    - this will create a NEW set for evolutionary analyses (BLAT complement)
- He will share the macse alignments as well
- he will run HyPhy
- **we will run RER converge (John)**
- David will use a Velifer centric approach for genes
- Manny will build the distribution of trees
    - also begin looking at structural changes
- David will send the Cactus alignments
- Elise will upload the species tree


## 2023/03/27

- Sudmant lab is running RER converge
    - *action item*: John finish that and present in next couple of weeks
- Tree building
    - We will move forward with using GTR model
    - Try neighbour joining as well and look for concordance
    - *action item*: Build the trees using GTR model and compare w/ neighbor joining (Manny)
- Elise is working on calibrating flexsuite, machine learning selection tests
- Dragon variant calling is running, comparing this to the original GATK bwa pipeline, (vs dragmap/dragvariant, https://github.com/Illumina/DRAGMAP)
- David reran the Hyphy analysis of selection w/ the new annotations, strong excess of adaptation at VIPs, esp strong if you focus on things that look like strong orthologs with human VIPs, RNA vs DNA VIPs, adaptation DNA viruses dominates very strongly
    - *action item* David present 
- Can we compare the Myotis lineage to other bat lineages?
    - approach: use diamond + miniprot to find orthologs among other bats, and then make gene trees from this
    - make gene alignments from all of those
    - run HyPhy again
    - *action item*: Run the gene discovery in "other" genomes, bats + HG38 (Manny), align and HyPhy (David), look at some examples (Lucie et al)
- Sedeph finished, but Biser did not run
    - sedeph was run on all
    - *action item*: Visualize and present on the SDs (Manny + Peter), gene family expansion and loss, Elise run Biser
    
- mobile element distribution 
    - hold off until next time
    - Peter and Manny looked at rates
    - Saba and Elise interested in looking @ PKR locus


