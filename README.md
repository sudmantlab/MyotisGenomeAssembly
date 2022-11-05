Assembly of 9 Novel *Myotis* Genomes
====================================

Mammals exhibit a wide variety of maximum lifespans; cells from exceptionally long-lived species 
like bats demonstrate exceptional tolerance stresses including DNA damage, oxidation, and heat shock. 
Comparative studies between long- and -short-lived species link increased stress tolerance with the 
improved cellular homeostasis underlying their extended lifespans and exceptional health and 
vitality in old age. However, factors such as poor-quality genomes, enormous evolutionary distances, 
and uncontrolled environmental and technical factors reduce the power these studies have to identify 
genetic and phenotypic differences between long- and short-lived species. To overcome these barriers, 
we are creating cell lines, reference genomes, and population genetic data from a clade of nine 
closely-related bats, both long- and short-lived relative to body size. 


Collaborators
-------

* **Juan M Vazquez** <a itemprop="sameAs" content="https://orcid.org/0000-0001-8341-2390"
href="https://orcid.org/0000-0001-8341-2390" target="orcid.widget" rel="me noopener noreferrer" 
style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" 
style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> 
\- [docmanny](https://vazquez.bio)
* **M Elise Lauterbur** <a itemprop="sameAs" content="https://orcid.org/0000-0002-7362-3618"
href="https://orcid.org/0000-0002-7362-3618" target="orcid.widget" rel="me noopener noreferrer" 
style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" 
style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> 
\- [lauterbur](https://github.com/lauterbur)
* 
* **David Enard** <a itemprop="sameAs" content="https://orcid.org/0000-0002-7362-3618"
href="https://orcid.org/0000-0002-7362-3618" target="orcid.widget" rel="me noopener noreferrer" 
style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" 
style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> 
\- [davidpierreenard](https://github.com/DavidPierreEnard)
* **Peter Sudmant** <a itemprop="sameAs" content="https://orcid.org/0000-0003-2634-8016" 
href="https://orcid.org/0000-0003-2634-8016" target="orcid.widget" rel="me noopener noreferrer" 
style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" 
style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>
\- [petersudmant](https://github.com/petersudmant) 
\- [Lab Website](https://www.sudmantlab.org)


Repo Structure: 
-------

- `assembly/`: Contains the scripts, Snakefiles, and other parts necessary for assembling all genomes
- `analyses/`: Contains all of the analyses being done for the paper in their own folders. 
Under debate whether to house each analysis in its own folder in this repo or to have git submodules
linked to folders here.
- `data/`: PLACE FINAL PRODUCTS <1GB IN HERE. 
Serves as a central spot to place files and outputs that are dependencies for other projects. 
Note that this might change if we agree on a better way to track and share dependencies.
- `manuscript/`: Holds the final figures and texts for the paper.
- `tools.bib`: Bibliography file of tools and scripts used herein. Under debate if this should be one
file or a `tools.bib` file in each folder.


For more information, please find a detailed README in each folder. 
