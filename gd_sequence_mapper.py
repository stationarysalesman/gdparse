__author__ = 'tyler'
"""Genomediff Parser. Copyright Tyler Camp.

This Python script takes genomediff files as input and outputs graphical information about
the mutations identified in the files."""

from Bio import SeqIO
import os
import re
from GenomeDiffSequenceMap import GenomeDiffSequenceMap


"""Pares genomediff files for statistical information about sample mutations.

This function defines input/output directories  used to gather
count-based data on mutations.

Precondition: Correctly formatted Genomediff files with tab-separated fields."""


def parse_files_cds(cat_map, input_dir, output_dir, plasmid_dir):
    """Parse information from .gd files into a map organized by CDS."""

    # Define string constants
    err_no_plasmid = "Error: no plasmid file found: "

    # Enumerated lists
    cfp_genes = ("Bba_E0020")
    yfp_genes = ("Bba_E0030", "Bba_K592101", "Bba_K864100")
    bfp_genes = ("Bba_K592100")

    # Map genes to category (this can be moved in future versions)
    gene_map = dict()
    for item in cfp_genes:
        gene_map[item.lower()] = "CFP"
    for item in yfp_genes:
        gene_map[item.lower()] = "YFP"
    for item in bfp_genes:
        gene_map[item.lower()] = "BFP"

    for dirName, subdirList, fileList in os.walk(input_dir):
        print subdirList
        for gdFile in fileList:

            with open(input_dir+gdFile, "r") as data:
                # Obtain a SeqRecord containing all info from Genbank file
                first_line = data.readline()
                second_line = data.readline()
                if not(second_line): # no mutations
                    return
                ref_seq_name = re.split("\t", second_line)[3]
                current_record = SeqIO.read(plasmid_dir+ref_seq_name+".gb", "genbank")
                if not(current_record):
                    print err_no_plasmid + ref_seq_name + "\n"
                    continue

                current_record_features = current_record.features

                # Gather information common to all mutations in this .gd file

                # Filter by strand (no duplicate features)
                top_strand_features = filter(item.strand == "1" for item in current_record_features)

                # Get type of CDS
                cds = filter(feat.type == "CDS" for feat in top_strand_features)
                cds_id = cds.qualifiers['label']
                cds_category = gene_map[cds_id]

                # Use this for all further count updates
                mutation_map = cat_map[cds_category]

                for line in data:
                    # Update total count
                    mutation_map.update_count()

                    split_line = re.split("\t", line)

                    # Common fields
                    mut_type = split_line[0]
                    ref_seq = split_line[3] # Plasmid sequence name

                    # Unique fields
                    position = ""
                    size = ""
                    repeat_seq = ""
                    repeat_length = ""
                    repeat_ref_num = ""
                    repeat_new_copies = ""
                    repeat_name = ""
                    strand = ""
                    # duplication_size = ""

                    # Update count based on mutation type
                    mutation_map.update_type_count(mut_type)

                    if ("MC" in mut_type): # Missing coverage follows unique format
                        #code to handle missing coverage to appear in later versions
                        pass

                    position = split_line[4]

                    # Update count based on feature (precondition: mutation has position)
                    # TODO: use FeatureLocation.start and .end to form sequence and check if mutation is positioned within
                    # Determine category

                    if ("SNP" in mut_type):
                         new_seq = split_line[5]


    return


"""Map individual reference sequences onto their respective category.


@:arg plasmid_dir: the directory containing reference sequence files
@:arg type: string specifying how to map reference files"""


def get_genbank_info(handle):

    # Map reference sequences onto CDS categories
    genbank_record = SeqIO.read(handle, "genbank")

    return
"""This Python script provides an interactive interface to a program that uses GenomeDiffSequenceMap objects.

Users are first asked how they wish to categorize sequences. Defaults are provided. Then the input and output
directories are specified. Then the analysis begins."""
def main():
    # Define string constants
    opening_message = ("gd_sequence_mapper - An interactive tool for processing Genomediff files\n",
                       "Copyright Tyler Camp. All rights reserved.\n")
    query_categorize = "Please specify how you want to categorize your sequences:\n"
    category_options = ("1. By coding sequence type (uses iGEM Spring 2015 CDS list)\n",
                        "2. By promoter strength (uses iGEM Spring 2015 promoter list)\n"
                        "3. By RBS strength (uses iGEM Spring 2015 CDS list)")
    err_bad_category = "Error: please try again.\n"
    err_input_dir_noexist = "Error: input directory does not exist.\n"
    err_input_dir_noaccess = "Error: input directory not accessible. Please check permissions and try again.\n"
    err_output_dir_fail = "Error: failed to create output directory. Please check permissions and try again.\n"
    err_output_dir_exists = "Error: output directory exists. Overwrite contents? [Y/n]\n"
    err_yesno = "Error: enter y or n.\n"
    # This dictionary maps categories to their GenomeDiffSequenceMap objects.
    cat_map = dict()

    # Category defaults (add more in later version)
    cds_categories = ("YFP", "BFP", "CFP")


    # Begin workflow
    print opening_message
    print query_categorize
    done = False
    while not(done):
        cat_num = input(("1. By coding sequence type (uses iGEM Spring 2015 CDS list)\n" +
                         "2. By promoter strength\n" +
                         "3. By RBS strength\n"))
        if (cat_num < 1 or cat_num > 3):
            print err_bad_category
            continue
        elif (cat_num == 1):
            print "yay"
            done = True
        else:
            print err_bad_category

    done = False
    while not(done):
        plasmid_dir = input("Please specify the plasmid directory path.\n")
        if not(os.access(plasmid_dir, os.F_OK)):
            print err_input_dir_noexist
            continue
        if not(os.access(plasmid_dir, os.EX_OK)):
            print err_input_dir_noaccess
            continue
        else:
            done = True

    done = False
    while not(done):
        user_input_dir = input("Please enter the input directory path.\n")
        # Check existence and access

        if not(os.access(user_input_dir, os.F_OK)):
            print err_input_dir_noexist
            continue
        if not(os.access(user_input_dir, os.EX_OK)):
            print err_input_dir_noaccess
            continue
        else:
            done = True

    done = False
    while not(done):
        user_output_dir = input("Please enter the output directory path.\n")
        if (os.access(user_output_dir, os.F_OK)):
            ans_done = False
            while not(ans_done):
                ans = input("Error: output directory exists. Overwrite contents? [Y/n]: ")
                if "y" in ans.lower():
                    ans_done = True
                    done = True
                    continue
                elif "n" in ans.lower():
                    ans_done = True
                    continue
                else:
                    print err_yesno
                    continue

    # Create the maps organized by category.

    # Categorize by CDS
    if (cat_num == 1):
        # Initialize each dictionary containing
        print "Initializing mutation tables...\r"
        for category in cds_categories:
            cat_map[category] = GenomeDiffSequenceMap()
        print "done."
        print "Parsing genomediff files.\n"
        parse_files_cds(cat_map, user_input_dir, user_output_dir, plasmid_dir)


    return


main()














