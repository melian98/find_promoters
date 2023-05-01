# import necessary packages
import ftplib
import os
import gzip
import re
from Bio import SeqIO
import random
import pandas as pd
import sys

# initialize needed variables
yes_no_zip = ''  # confirm if the user has downloaded their needed gunzip files from ftp ensembl
i = 0  # counter for the while loops
gene_information = []  # contain necessary gene information
fasta = False
gff3 = False
DNA_sequence = ""
TSS = ""
chromosome = ""


# function to download gz files from the ftp server given the name of the desired file as well as the server
def get_ftp(filename, server):
    fname = filename
    # checks if the gunzip file is already present in the current directory in which case no download is required
    if not os.path.isfile(os.getcwd() + '/' + fname):
        # asks user for username, password and the filepath for the desired file
        username = input("please enter your username : ")
        password = input("please enter your password : ")
        filepath = input("please enter the filepath of the file you wish to download : ")
        # connect to the server
        ftp = ftplib.FTP(server)
        # logs in using the users information
        ftp.login(user=username, passwd=password)
        # sets the encoding to ensure it is in utf-8
        ftp.encoding = "utf-8"
        # navigates to the inputted filepath
        ftp.cwd(filepath)
        # writes the file onto the current working directory
        with open(fname, 'wb') as file:
            ftp.retrbinary(f"RETR {fname}", file.write)


# function which will check for unknown nucleotides (represented by N) and only keep the sequence after the last N
def check_for_n(seq):
    # confirms if N is in the sequence
    if 'N' in seq:
        # finds the index of the last instance of the letter N
        index = seq.rfind("N")
        # retains only the part of the sequence after the last N
        seq = seq[index + 1:]
    # returns the new sequence
    return seq


# function which will take the reverse complement of a given sequence
def rev_comp(seq):
    # replaces each nucleotide with their complement. These are not capitalized to ensure that they won't then be
    # replaced again
    rev = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    # capitalizes the new complementary sequence
    rev = rev.upper()
    # reverses the sequence so that the reverse complement is now obtained
    rev = rev[::-1]
    # returns the reverse complement sequence
    return rev

#check if the directory contains a fasta and gff3 file, then opens them
dir_list = os.listdir(".")
for diry in dir_list:
    if diry.endswith('.fa'):
        fasta = True
        chromosome = open("./" + diry, mode='rt')

    elif diry.endswith('.gff3'):
        gff3 = True
        TSS = open("./" + diry, mode='rb')


# if the gunzip file has either not been downloaded or has been downloaded but not extracted, it will be downloaded
# from the ftp ensembl website using the users username and password
if not fasta or not gff3:
    #inform the user that the files were not found and ask if they would like to download them from ftp ensemble
    download = input("Either the Fasta or the gff3 file were not present in the current working directory. "
                     "Would you like to download them from the ftp ensemble website (input y for yes or n for no) : ")
    # ask user for the gunzip file name for the fasta file
    if download == 'y':
        chromosome_gzip = input('please enter the name of the chromosome fasta file from ftp ensemble that you wish to '
                                'download : ')

        # call get_ftp to check if the file is downloaded and download it if it has not been downloaded
        get_ftp(chromosome_gzip, "ftp.ensemblgenomes.org")

        # ask user for the gunzip file name of the gff3 file
        TSS_gzip = input('please input the name of the gff3 file : ')

        # call get_ftp to check if the file is downloaded and download it if it has not been downloaded
        get_ftp(TSS_gzip, "ftp.ensemblgenomes.org")

        # extract and open the fasta file
        chromosome = gzip.open(chromosome_gzip, mode='rt')

        # extract and open the gff3 file
        TSS = gzip.open(TSS_gzip, mode='rb')
    #if user doesn't want to download the file from ftp ensemble, exit the program
    else:
        sys.exit("thank you for using the find_promoters tool")

# open and read in the genes text file
gene_list = open('genes.txt', 'r')
genelist_lines = gene_list.readlines()

# open the list of promoters and save them into a list using readlines
promoters = open('promoters.txt', 'r')
promoters = promoters.readlines()

# create a list containing the lines for both the gff3 and gene list files
TSS_lines = TSS.readlines()

# filter the gff3 file to remove all lines that don't contain a gene
while i < len(TSS_lines):
    # string will contain the line as the while loop loops through each line
    string = str(TSS_lines[i])
    # convert the gff3 line into a list delimiting by tabs
    string_list = re.split(r'\\t', string)
    # all lines containing gene information will be converted into a list longer than 4, only header
    # information will be filtered out
    if len(string_list) > 4:
        # string_list[2] contains the feature type (exon, CDS, gene, mRNA, etc.) only gene is wanted
        if str(string_list[2]) != 'gene':
            # if string_list[2] was anything other than gene, the given line is deleted from the TSS file
            del TSS_lines[i]
            # the increment of i is subtracted by one because by deleting the previous line, i+1 has turned into i
            i = i - 1
    else:
        # if string_list was shorter than 4, then this won't contain gene information and will be deleted
        del TSS_lines[i]
        # the increment of i is subtracted by one again because a deletion was made
        i = i - 1
    # i is incremented by one to continue the loop
    i = i + 1

# i is reset to zero for the next loop
i = 0

# while loop will delete all genes from gene_list that are not present in the gff3 file
while i < len(genelist_lines):
    # substring contains the gene name from gene_list converted to a string
    substring = str(genelist_lines[i])
    # removes all meta character such as \n\t and others
    substring = re.sub('[^A-Za-z0-9]+', '', substring)
    # delete = 1 will cause the current line in gene list to be deleted
    delete = 1
    # loop through TSS lines to see if the gene in gene_list is present
    for j in range(len(TSS_lines)):
        # string is set  to be the string of the current line of the gff3 file converted to a string
        string = str(TSS_lines[j])
        # if the gene in gene list isn't found in the gff3 file, delete is set to false (zero)
        if string.find(substring) != -1:
            delete = 0
            # break the loop because it's confirmed that this gene is in the gff3 file
            break
    # if the gene wasn't found in the gff3 file, delete will be equal to one
    if delete == 1:
        # delete the current genelist line
        del genelist_lines[i]
        # because a line was deleted, i+1 is now i, so it is reduced by one
        i = i - 1
    # increment i by one
    i = i + 1

# loop will populate the gene_information list with the name of the gene, its start and end position as well as if
# it is on the positive or negative strand
for i in range(len(genelist_lines)):
    # remove all special characters from the genes being tested
    genelist_lines[i] = re.sub('[^A-Za-z0-9]+', '', genelist_lines[i])
    # loop through the lines in the gff3 file and find the line that contains the gene being tested
    for j in range(len(TSS_lines)):
        if str(TSS_lines[j]).find(genelist_lines[i]) != -1:
            # split the gff3 line into a list and extract the desired properties
            string_list = re.split(r'\\t', str(TSS_lines[j]))
            # append gene information with the name, start and end position and the strand
            append_list = [str(genelist_lines[i]), string_list[3], string_list[4], string_list[6]]
            gene_information.append(append_list)
            # break from the second for loop since a match was found for this gene
            break

# parse the fasta file and save the nucleotide sequence into DNA_sequence
for seq_record in SeqIO.parse(chromosome, "fasta"):
    DNA_sequence = seq_record.seq

# loop will extract the nucleotide sequence either upstream of the start site or downstream of the end site as well as
# take the reverse complement if necessary and check for N's
for i in range(0, len(genelist_lines)):
    # check if the gene is on the negative strand
    if gene_information[i][3] == '-':
        # extract the nucleotides downstream of the end site denoted in the gff3 file
        sequence = DNA_sequence[int(gene_information[i][2]):int(gene_information[i][2]) + 400]
        # take the reverse complement of the sequence
        sequence = rev_comp(sequence)
    # if the gene is on the positive strand
    else:
        # extract the nucleotides upstream of the start sire
        sequence = DNA_sequence[int(gene_information[i][1]) - 400:int(gene_information[i][1])]
    # check the extracted sequence for nucleotides and keep nucleotides downstream of the last N
    sequence = check_for_n(sequence)
    # append gene information with the sequence
    gene_information[i].append(sequence)

column_names = genelist_lines
column_names.insert(0, "Promoters")
final_counts = pd.DataFrame(columns=column_names)

for i in range(0, len(promoters)):
    promoters[i] = re.sub('[\\n]', '', promoters[i])
    count = [promoters[i]]

    for j in range(0, len(gene_information)):
        count.append(len(re.findall(promoters[i], str(gene_information[j][4]))))

    new_row = pd.DataFrame({k: v for k, v in zip(column_names, count)}, index=[0])

    final_counts = pd.concat([final_counts.loc[:], new_row]).reset_index(drop=True)

path = r'.\results.txt'

if os.path.isfile(path):
    os.remove("results.txt")

with open(path, 'a') as f:
    df_string = final_counts.to_string(header=True, index=False)
    f.write(df_string)

f.close()

# tell the user that the code is done and the results.txt file can be accessed.
print('done. The results can be found in the results.txt file in the '
      'current working directory. Thank you for using the find_promoters tool.')
