Assembly of 9 Novel *Myotis* Genomes
====================================

Mammals exhibit a wide variety of maximum lifespans; cells from exceptionally long-lived species like bats demonstrate exceptional tolerance stresses including DNA damage, oxidation, and heat shock. Comparative studies between long- and -short-lived species link increased stress tolerance with the improved cellular homeostasis underlying their extended lifespans and exceptional health and vitality in old age. However, factors such as poor-quality genomes, enormous evolutionary distances, and uncontrolled environmental and technical factors reduce the power these studies have to identify genetic and phenotypic differences between long- and short-lived species. To overcome these barriers, we are creating cell lines, reference genomes, and population genetic data from a clade of nine closely-related bats, both long- and short-lived relative to body size. 


Collaborators
-------

* **Juan M Vazquez** <a itemprop="sameAs" content="https://orcid.org/0000-0001-8341-2390"
href="https://orcid.org/0000-0001-8341-2390" target="orcid.widget" rel="me noopener noreferrer" 
style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" 
style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> 
\- [docmanny](https://vazquez.bio)
* **Elise Lauterbur** 
\- [lauterbur](https://github.com/lauterbur)
* **David Enard**
\-
* **Peter Sudmant** <a itemprop="sameAs" content="https://orcid.org/0000-0002-9573-8248" 
href="https://orcid.org/0000-0002-9573-8248" target="orcid.widget" rel="me noopener noreferrer" 
style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" 
style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>
\- [petersudmant](https://github.com/petersudmant) 
\- [Lab Website](https://www.sudmantlab.org)


Repo Structure: 
-------

- `data`: Contains all initial data for assembly, including sequencing data;
- `envs`: Contains the Conda environment files defining the tools for assembly;
- `output`: Contains all the intermediate and final outputs of assembly;
- `rules`: Contains the Snakefiles with rules for assembling and QC-ing the genomes;
- `Snakefile`: The master Snakefile for assembling the genomes.

