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


def parse_files(cat_map, input_dir, output_dir):

    for dirName, subdirList, fileList in os.walk(input_dir):
        print subdirList
        for gdFile in fileList:
            with open(input_dir+gdFile, "r") as data:
                for line in data:
                    print line

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
        user_input_dir = input("Please enter the input directory path.\n")
        # Check existence and access
        print "Entered", str(user_input_dir)
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
        parse_files(cat_map, user_input_dir, user_output_dir)


    return


main()














