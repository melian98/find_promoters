# find_promoters
The find promoters tool will find how many instances of a promoter sequence are upstream of a gene. This tool will output the results in a text file in tabular format with the genes as column headers and the promoter sequences as row names. 

The infividual cells will count how many times that promoter sequence was found in the gene. In order to run, the python code must be provided with a genes text file, a promoters text file, and a gff3 and fasta file (must be .fa). An example of the inputs as well as the output results text file are provided in this repository. 

If a gff3 and fasta file are not present when the tool is run, the tool will allow the user to download the required files from "ftp.ensemblgenomes.org" in which case the user must provide their username, password and the names and filepaths of the fasta and gff3 file that they wish to download.

The files present in this repository are as follows:

find_promoters.py : A .py file which when run, will execute the code and generate the results.txt file given the appropriate genes, promoters, fasta and gff3 files

genes.txt: A single column of properly formatted gene names which will be cross referenced with the gff3 file. The file in this repository serves as an example and can be replaced with any other properly formatted genes file in which case that file must be renamed to "genes.txt"

promoters.txt: A single column of properly formatted promoter sequences which will be used to search for the promoter sequence in the fasta file. Sequences found in square brackets (i.e AC[AG]) will be interpreted as being one or the other (i.e ACA and ACG). The file in this repository serves as an example and can be replaced with any other properly formatted promoter file in which case that file must be renamed to "promoters.txt"

Zea_mays.B73_RefGen_v4.48.chromosome.1.gff3 : a gff3 file provided as an example but can be replaced with any properly formatted .gff3 file

Zea_mays.B73_RefGen_v4.dna.chromosome.1.fa : a fasta file provided as an example but can be replaced with any properly formatted .fa file 
