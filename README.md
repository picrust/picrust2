# PICRUSt2

The current code and documentation for PICRUSt2 is found here. **Please note that PICRUSt2 is still actively being developed and tested**, so in the meantime we recommend you use [PICRUSt1](http://picrust.github.io/picrust/).

Documentation is found on this [repository's wiki](https://github.com/picrust/picrust2/wiki).

We are planning to release PICRUSt 2 by the end of 2018. It includes these and other improvements:
* Updated pre-calculated files with additional genomes and gene families.
* Allows output of MetaCyc ontology predictions that will be comparable with common shotgun metagenomics outputs.
* Straight-forward running of genome prediction, which will allow users to predict functions for all their study sequences (even those not found in a database).
* Addition of hidden-state prediction algorithms from the ```castor``` R package.

# Citations
PICRUSt2 wraps a number of tools to generate functional predictions from amplicon sequences. If you use PICRUSt2 you also need to cite these tools:

### For phylogenetic placement of reads:
* [PaPaRa](https://academic.oup.com/bioinformatics/article/27/15/2068/400617)
* [EPA-NG](https://www.biorxiv.org/content/early/2018/03/29/291658)
* [gappa](https://github.com/lczech/gappa), which is based on the [genesis library](https://github.com/lczech/genesis).

### For hidden state prediction:
* [castor](https://academic.oup.com/bioinformatics/article/34/6/1053/4582279)

### For pathway inference:
 * [MinPath](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000465) - A modified version of this tool from the HMP project is also packaged with PICRUSt2. This tool was releasd under the GNU General Public License. [Source code](http://omics.informatics.indiana.edu/MinPath/).
