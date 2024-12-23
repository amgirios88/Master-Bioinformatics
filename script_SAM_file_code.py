# Main task is to extract split reads from an alignment file (sam file) to identify the locations of introns in genes (in a tab-separated file)

# script structure: modules, functions, command line arguments

# usage: python3 Assignment1_3030116G.py mySamFile.sam myInputTable.txt

# ***************************************
# *  Modules to load
# ***************************************

import sys
import re


# ***************************************
# *  Processing SAM file
# ***************************************

# defining a function to process SAM files and analyze read alignments
def process_sam_file(sam_file):
    """
    Function that processes a SAM file to extract split reads (junctions) and count supporting reads for each junction into a dictionary. 
    Key for the dictionary is a tuple (chromosome, junction_start, junction_end) and the value is the count of supporting reads.
    """
    # initialising an empty dictionary to store junctions and their supporting read counts
    junctions = {}

    try:
        with open(sam_file, 'r') as sam:
            for line in sam:
                if line.startswith('@'): # skipping the header
                    continue

                # processing the SAM file line by line
                line = line.strip().split('\t') # splitting the line by tab to get data for each column
                rname = line[2]  # chromosome to which the read aligned (column 3 in SAM file)
                pos = int(line[3])  # position on the chromosome where the alignment starts (column 4 in SAM file)
                cigar = line[5]  # CIGAR string; string describing the alignment in cigar format (column 6 in SAM file)
                if line[-1].startswith("NH:i:"): # checking if the last column starts with NH:i: since NH column is an optional column
                    nh_column = line[-1]  # assigning NH column when present
                else:
                    continue # skipping reads that do not have NH: in the last column

                # checking if the read aligns only once (NH:i:1), to process them further
                if nh_column.split(':')[-1] != '1':
                    continue  # skipping reads that do not align exactly once
                
                if 'N' in cigar: # checking if the read is split (at least one N in the cigar string) to process them further
                    # processing the CIGAR string to extract the junctions
                    matches = re.findall(r'(\d+)([MIDNS])', cigar) # literal pattern is: one or more digits followed by one (and one only) of the characters M, I, D, N, or S

                    # checking if the read is split (at least one N in the cigar string)
                    current_pos = pos  # position on the reference genome
                    for length, operation in matches:
                        length = int(length)

                        if operation in ['M', 'D']:  # I and S operations do not affect the position, so we can ignore them
                            current_pos += length # adds the length of the match or deletion to the current position
                        elif operation == 'N':  # N means skipped region and indicates a junction
                            junction_start = current_pos # start position of the junction equals the current position
                            current_pos += length # move the current position to the end of the junction, adding the length of the skipped region
                            junction_end = current_pos # end position of the junction equals the current position

                            # key for the junction
                            key = (rname, junction_start, junction_end)

                            # incrementing the count for the junction to check how many reads support it
                            if key in junctions: # using the previously defined dictionary to store the junctions and their supporting read counts
                                junctions[key] += 1 # increment the count if the junction is already in the dictionary
                            else:
                                junctions[key] = 1 # add the junction to the dictionary with a count of 1 if it is not already present

    except FileNotFoundError:
        print(f"SAM file {sam_file} not found.")

    except Exception as e:
        print(f"An error occurred while processing SAM file: {e}")
    
    # added an additional error handling step, since the first time I forgot to close an } in the function processing tab-separated files, and returned an empty list that
    # was not iterable in the junctions_in_genes(). This prints a statement to confirm that the file was processed successfully (checks if junctions is not empty to print the statement)
    if junctions:
        print(f"SAM file was processed successfully.")

    return junctions
 

# ChatGPT 4o was used for the following tasks:
    # clarifying if I understood some tasks correctly (i,e, I asked "Are we counting number of reads supporting each junction for reads that only 
    #   align once?" inputting the instructions for each data row in the sam file, to make sure I was understanding the task correctly)
    # if line[-1].startswith("NH:i:"): ; checking if this method could be applied to the last item in the line tuple, which I was not completely sure about
    # if nh_column.split(':')[-1] != '1': ; was suggested by chatGPT when I asked if my code made sense, which was a longer version (and I liked the shorter version better):
        # nh_parts = nh_column.split(':')
        # nh_value = nh_parts[-1]
        # if nh_value != '1':
            # continue
    # if it made sense to use both exceptions in the same try block, which I was not sure about; I also checked class notes and Python documentation


# function for processing the genes file
def process_genes_table(file_name):
    """
    Function that processes a tab-separated file to extract gene_id, location, chromosome, start, end, and strand.
    """
    genes = []  # initialising list to store gene details, where each gene is a dictionary, we don't need to have unique identifiers for the genes 
                # because we will loop through all of them
    try:
        with open(file_name, 'r') as file:
            header = file.readline()  # this will skip the header

            # processing the tab-separated file line by line
            for line in file:
                columns = line.strip().split('\t')
                gene_id = columns[0]
                # transcript_id was not mentioned in the tasks, so I skipped it which would be columns[1]
                location = columns[2] 

                # this is a complex string, so we need to split firt the chromosome name (chrom) and the rest of the string
                chrom, rest = location.split(':')
                # then we split the "rest" string to get the start, which is separated from end with '..' and the strand is separated from end with '('
                start, end = rest.split('..') # end here still includes the strand, so we need to further separate it
                end, strand = end.split('(') # strand will still have the ')' at the end, so we need to remove it
                strand = strand.strip(')')
                # removing commas and converting position numbers to integers for start and end
                start = int(start.replace(',', ''))
                end = int(end.replace(',', '')) 

                # appending the gene details to the list previously initialised; each gene is a dictionary with the following keys: gene_id, chrom, start, end, strand
                genes.append({
                    "gene_id": gene_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand
                })
    except FileNotFoundError:
        print(f"Error: tab-separated file {file_name} not found.")
    
    except Exception as e:
        print(f"An error occurred while processing tab-separated file: {e}")
    
    # added an additional error handling step, since the first time I forgot to close an } in the function processing tab-separated files, and returned an empty list that
    # was not iterable in the junctions_in_genes(). This prints a statement to confirm that the file was processed successfully (checks if genes is not empty to print the statement)
    if genes:
        print(f"Tab-separated file was processed successfully.")

    return genes
   

# ChatGPT 4o was used for the following tasks:
    # I double checked if I needed the transcript_id in the genes list, because the column was mentioned in the structure but not in the tasks, so I just
    #   wanted to make sure I was right in not needing the transcript_id data, so I removed it
    # I also used ChatGPT to troubleshoot an error I was having "'NoneType' object is not iterable" when I tried to run the script, and it happened because
    #   I forgot to return the genes list at the end of the function, so it was returning None, and I was trying to iterate over None in the junctions_in_genes() function

# processing the junctions and genes to find introns
def junctions_in_genes(genes, junctions, output_file):
    """
    Finds junctions within gene boundaries and writes them to a tab-separated file. This function takes tree arguments: genes (list of dictionaries - output from 
    process_gene_table function), junctions (dictionary - output from process_sam_file function), and output_file (which will contain the junctions that are 
    within the boundaries of each gene).
    """
    try:
        with open(output_file, 'w') as output:
            output.write("gene_id\tjunction_start\tjunction_end\tcount\n")
            # iterating trough the dictionary for each gene in the genes list
            for gene in genes:
                gene_id = gene["gene_id"] # these access the mentioned key in the dictionary and assigns it to the named variable
                chrom = gene["chrom"]
                start = gene["start"]
                end = gene["end"]

                # finding all junctions within the boundaries of the gene
                # this code loops through the junctions dictionary, identifies juntions that match the current gene chromosome and are within the gene boundaries, and 
                # adds them to the matching_junctions list
                # conditions to filter the junctions: if same chromosome, and junction start is within the gene boundaries (<= greater or equal to start and >= smaller or 
                # equal to end) and junction end is within the gene boundaries (<= greater or equal to start and >= smaller or equal to end)
                matching_junctions = []
                matching_found = False
                for (j_chrom, junction_start, junction_end), count in junctions.items():
                    if j_chrom == chrom and start <= junction_start <= end and start <= junction_end <= end:
                        matching_junctions.append((junction_start, junction_end, count))
                        matching_found = True


                # writing the junctions for the current gene; the corresponding gene_id will be added because of the nested iteration, since each gene is processed separately 
                # the matching_junctions list will be iterated for each gene
                for junction_start, junction_end, count in matching_junctions:
                    # will write a single line for each junction
                    output.write(f"{gene_id}\t{junction_start}\t{junction_end}\t{count}\n") 
                    # writing the gene_id, junction start, junction end, and the count of reads supporting the junction
                
                # adding an empty row after each gene before continuing to the next gene
                if matching_found == True:
                    output.write("\n")
                else:
                    continue
                # the idea for adding a boolean variable to check if junctions were found to then add or not an empty row was
                # provided by Kathryn during a lab class

    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")
    print("File written successfully.") 
    return output_file

# ChatGPT 4o was used for the following tasks:
    # I first created a list comprehension to filter the junctions that are within the gene boundaries, but I was not sure if it was the best approach. This is the list comprehension:
                 # matching_junctions = [
                    # (junction_start, junction_end, count) # this is the tuple that will be added to the list
                    # for (j_chrom, junction_start, junction_end), count in junctions.items() # loop that goes through the junctions dictionary, .items() accesses the key-value pairs
                    # the first part matches the key structure, and count is the associated value to each key
                    # if j_chrom == chrom and start <= junction_start <= end and start <= junction_end <= end # conditions to filter the junctions
                # ]
    # However, since I was not 100% sure that I wrote it okay, I changed it for a for loop that appends the junctions that match the conditions to the matching_junctions list, which I was more comfortable with
    # and asked chatGPT if the function would work as expected anyway, and it was confirmed that it should work as expected
    # I also used it to debug this functions, since I was having an issue with the header not being written, and it was because my code to write the header was inside the loop which I did not realise

# I also used this website (https://blogs.glowscotland.org.uk/sh/ahscomputingpython/higher/writing-data-to-a-text-file/#:~:text=The%20file.,%5C”%20should%20also%20be%20written.)
# as support to write the output to a file

# code to deal with the command line arguments
if len(sys.argv) != 3: # checking if the number of arguments is correct, if not, the script will print the usage and exit; input should only be the sam file and the gene table
    print("Usage: python3 myScript.py mySamFile.sam myInputTable.txt") # printing usage if the number of arguments is incorrect
    sys.exit(1) # and exiting the script

# assigning the command line arguments to variables
sam_file = sys.argv[1]
gene_file = sys.argv[2]
output_file = "3030116G.txt" # output file name with my student number (hardcoded)

try:
    # processing the SAM file to get junctions using the function defined above (process_sam_file) and storing the results in the junctions dictionary
    junctions = process_sam_file(sam_file)

    # processing the tab-separated gene table to get gene details using the function defined above (process_gene_table) and storing the results in the genes list
    genes = process_genes_table(gene_file)

    # finding junctions within the boundaries of genes using the function defined above (junctions_in_genes) and writing the results to a tab-separated file hardcoded as 3030116G.txt
    junctions_in_genes(genes, junctions, output_file)
    print(f"Results written to {output_file}")

except Exception as e:
    print(f"An unexpected error occurred running the script: {e}")

# ChatGPT 4o was used for the following tasks:
    # I asked if the code structure was correct, or if I should place the block dealing with command line arguments at the beginning of the script, 
    # and it was confirmed that it was okay where it was
    # I also asked if the way I wrote the output file name was hardcoding, just to double check and be sure, and it was confirmed that it was okay
    # I also used it to double check the "write good code" part of the task, to make sure I was following the best practices
