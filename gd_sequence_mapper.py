__author__ = 'tyler'
"""Genomediff Parser. Copyright Tyler Camp.

This Python script takes genomediff files as input and outputs graphical information about
the mutations identified in the files.

This file is part of gdparse.

    gdparse is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gdparse is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with gdparse.  If not, see <http://www.gnu.org/licenses/>."""

from Bio import SeqIO
import os
import re
from collections import Counter
from GenomeDiffSequenceMap import GenomeDiffSequenceMap

"""get_category(): Return a key to use in the mapping structure that exists in the caller.


This function takes input in the form of an open file handle (GenomeDiff file), the current template, features, and a
categorization number which determines the logic for determining the key string to return."""
def get_category(data, current_record, top_strand_features, categorization_number):

     # Error strings
    err_no_cds = "Error: no cds found in current template."

    # Enumerated lists
    cfp_genes = ("bba_e0020",)
    yfp_genes = ("bba_e0030", "bba_k592101", "bba_k864100")
    bfp_genes = ("bba_k592100",)
    cds_tuple = tuple(tuple(cfp_genes) + tuple(yfp_genes) + tuple(bfp_genes))

    identified_category = '' # we will use this as key in map in calling function
    if categorization_number == 1 or categorization_number == 2:
        # Map genes to category (this can be moved in future versions)
        gene_map = dict()
        for item in cfp_genes:
            gene_map[item.lower()] = "CFP"
        for item in yfp_genes:
            gene_map[item.lower()] = "YFP"
        for item in bfp_genes:
            gene_map[item.lower()] = "BFP"

        # Get type of CDS
        cds = filter(lambda feat: feat.type == "CDS" and (feat.qualifiers['label'][0]).lower() in cds_tuple, top_strand_features)
        if not cds:
            print err_no_cds
            return 'None'
        cds = cds[0] # get object from list
        cds_id = cds.qualifiers['label'][0].lower()
        if categorization_number ==1:
            identified_category = gene_map[cds_id]
        identified_category = cds_id

    return identified_category
"""parse_file_data(): parses information contained within .gd files.


This function encapsulates the map update functionality of the program to simplify the map update process.
@input: data from a GenomeDiff file
@output: GenomeDiffSequenceMap object
"""


def parse_file_data(data, mutation_map, features, category):
    # Get cutoff
    cutoff = None
    curr_cds = filter(lambda feat: category in (feat.qualifiers['label'][0]).lower(), features)[0]
    if curr_cds:
        cutoff = int(curr_cds.location.end) - 300
        print "cutting off at", str(cutoff)
    for line in data:

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


        if (mut_type == "MC"): # Missing coverage follows unique format
            #code to handle missing coverage to appear in later versions
            pass

        position = int(split_line[4])
        # Ignore mutations after the cutoff
        if cutoff and (position > cutoff):
            continue
        # Update count based on feature
        containing_feature_type = None
        containing_feature_label = None
        containing_feature = filter(lambda feat: position in feat.location, features)
        if (containing_feature):
            containing_feature_type = containing_feature[0].type
            containing_feature_label = containing_feature[0].qualifiers['label'][0]
            if 'BBa_K608002' in containing_feature_label:
                print "found problem in", data, "at nt pos", str(position)
        else: # not within an annotation; we are currently counting these
            #containing_feature_type = 'None'
            containing_feature_type = "None"
            containing_feature_label = "None"
            print "Found mutation outside of annotations."
            print "IN:", str(data)
            print "TYPE:", mut_type
            print "AT:", str(position)
        mutation_map.update_feature_map(containing_feature_label)
        """Add a mapping of this label to its type for output."""
        mutation_map.update_label_type_map(containing_feature_label, containing_feature_type)

        # Update other counts based on detailed mappings
        mutation_map.update_type_feat_map(mut_type, containing_feature_label)
        mutation_map.update_feat_type_map(containing_feature_label, mut_type)

        # Update count based on mutation type
        mutation_map.update_type_map(mut_type)
        # Update total count
        mutation_map.update_count()

    return mutation_map

"""Parse genomediff files for statistical information about sample mutations.

This function defines input/output directories  used to gather
count-based data on mutations.

Precondition: Correctly formatted Genomediff files with tab-separated fields."""


def parse_files_cds(cat_map, categorization_number, input_dir, output_dir, plasmid_dir):
    """Parse information from .gd files into a map organized by CDS."""

    # Define string constants
    err_no_plasmid = "Error: no plasmid file found: "
    err_no_category = "Error: sample category not defined."

    print "Scanning input directory..."
    file_count = 0
    for dirName, subdirList, fileList in os.walk(input_dir):
        for f in fileList:
            file_count +=1
    print "Found", file_count, "files."
    print "Processing files:"
    for dirName, subdirList, fileList in os.walk(input_dir):
        for gdFile in fileList:
            with open(input_dir+gdFile, "r") as data:
                # Obtain a SeqRecord containing all info from Genbank file
                first_line = data.readline()
                second_line = data.readline()
                if not(second_line): # no mutations
                    continue
                data.seek(18) # return to beginning of second line
                ref_seq_name = (re.split("\t", second_line)[3]).lower()
                current_record = SeqIO.read(plasmid_dir+ref_seq_name+".gb", "genbank")
                if not(current_record):
                    print err_no_plasmid + ref_seq_name + "\n"
                    continue
                # Filter by strand (no duplicate features)
                top_strand_features = filter(lambda item: item.strand == 1, current_record.features)
                # Determine category key in category map
                category = get_category(data, current_record, top_strand_features, categorization_number)
                if not cat_map[category]:
                    print err_no_category
                    continue
                temp_map = parse_file_data(data, cat_map[category], top_strand_features, category)
                
                if (temp_map):
                    cat_map[category] = temp_map




    return cat_map


def parse_file_labels(cat_map, input_dir, output_dir, plasmid_dir):
    """Parse samples based on labels."""
    print "Scanning input directory..."
    file_count = 0
    for dirName, subdirList, fileList in os.walk(input_dir):
        for f in fileList:
            file_count +=1
    print "Found", file_count, "files."
    print "Processing files:"
    for dirName, subdirList, fileList in os.walk(input_dir):
        for gdFile in fileList:
            with open(input_dir+gdFile, "r") as data:
                print gdFile, "...\r"
                # Obtain a SeqRecord containing all info from Genbank file
                first_line = data.readline()
                second_line = data.readline()
                if not(second_line): # no mutations
                    continue
                data.seek(18) # return to beginning of second line
                ref_seq_name = (re.split("\t", second_line)[3]).lower()
                current_record = SeqIO.read(plasmid_dir+ref_seq_name+".gb", "genbank")
                if not(current_record):
                    print "nope: " + ref_seq_name + "\n"
                    continue
                # Filter by strand (no duplicate features)
                top_strand_features = filter(lambda item: item.strand == 1, current_record.features)
                temp_map = parse_file_data(data, cat_map['all'], top_strand_features, 'CDS')
                if temp_map:
                    cat_map['all'] = temp_map

    return cat_map
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

    # Enumerated lists
    cfp_genes = ("bba_e0020",)
    yfp_genes = ("bba_e0030", "bba_k592101", "bba_k864100")
    bfp_genes = ("bba_k592100",)
    cds_tuple = tuple(tuple(cfp_genes) + tuple(yfp_genes) + tuple(bfp_genes))

    # Defaults for user input
    INPUT_DIR_DEFAULT = "genomediff/"
    OUTPUT_DIR_DEFAULT = "output/"
    PLASMID_DIR_DEFAULT = "plasmids/"


    # Begin workflow
    print opening_message
    print query_categorize
    done = False
    while not(done):
        cat_num = input(("1. By coding sequence type (uses iGEM Spring 2015 CDS list)\n" +
                         "2. By specific CDS (iGEM Spring 2015 CDS list)\n" +
                         "3. By label\n" +
                         "4. By RBS strength\n"))
        if (cat_num < 1 or cat_num > 4):
            print err_bad_category
            continue
        elif(cat_num in range(1, 5)):
            done = True
        else:
            print err_bad_category

    done = False
    user_plasmid_dir = ""
    while not(done):
        user_input = raw_input("Please specify the plasmid directory path.\n1. Default directory "+
                               "(local)\n2. Custom\n")
        if (user_input == "1"):
            user_plasmid_dir = PLASMID_DIR_DEFAULT
            done = True
            continue

        if not(os.access(user_plasmid_dir, os.F_OK)):
            print err_input_dir_noexist
            continue
        if not(os.access(user_plasmid_dir, os.EX_OK)):
            print err_input_dir_noaccess
            continue
        else:
            done = True

    done = False
    user_input_dir = ""
    while not(done):
        user_input = raw_input("Please enter the input directory path.\n1. Default directory (local)\n2. " +
                               "Custom\n")
        if (user_input == "1"):
            user_input_dir = INPUT_DIR_DEFAULT
            done = True
            continue
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
    user_output_dir = ""
    while not(done):
        user_input = raw_input("Please enter the output directory path.\n1. Default directory (local)\n2. " +
                               "Custom\n")
        if (user_input == "1"):
            user_output_dir = OUTPUT_DIR_DEFAULT
            done = True
            continue

        if (os.access(user_output_dir, os.F_OK)):
            ans_done = False
            while not(ans_done):
                ans = raw_input("Error: output directory exists. Overwrite contents? [Y/n]: ")
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
    categorization_number = 0
    # Categorize by CDS
    if (cat_num == 1):
        # Initialize each dictionary containing
        for category in cds_tuple:
            cat_map[category] = GenomeDiffSequenceMap()
        categorization_number = 1

    # Categorize by specific CDS
    if (cat_num == 2):
        for cds in cds_tuple:
            cat_map[cds] = GenomeDiffSequenceMap()
        categorization_number = 2

    # Categorize by label
    if (cat_num == 3):
        cat_map['all'] = GenomeDiffSequenceMap()
        categorization_number = 3

    new_map = None
    print "Parsing genomediff files.\n"
    if categorization_number < 3:
        new_map = parse_files_cds(cat_map, categorization_number, user_input_dir, user_output_dir,
                                  user_plasmid_dir)
    elif categorization_number == 3:
        new_map = parse_file_labels(cat_map, user_input_dir, user_output_dir, user_plasmid_dir)

    done = False
    user_input = ""
    while not(done):
        user_input = raw_input("Choose output type.\n1. CSV: Organize by mutation type\n2. " +
                               "CSV: Organize by label\n")
        if (user_input == "1"):
            output_muttype_csv()
            done = True
            continue
        if (user_input == "2"):
            output_label_csv(new_map)
            done = True
            continue


def output_label_csv(new_map):

    master_feat_type_map = dict()
    master_counter = Counter(master_feat_type_map)
    x = 1
    for k in new_map.keys():
        with open("output"+str(x)+".csv", "a") as of:
            header= ",MOB,INS,DEL,SNP,TOTAL\n"
            of.write(header)
            for feat in new_map[k].feat_type_map.keys():
                curr_mob_count = 0
                curr_ins_count = 0
                curr_del_count = 0
                curr_snp_count = 0
                curr_feat = new_map[k].feat_type_map[feat]
                if 'MOB' in curr_feat:
                    curr_mob_count = curr_feat['MOB']
                if 'INS' in curr_feat:
                    curr_ins_count = curr_feat['INS']
                if 'DEL' in curr_feat:
                    curr_del_count = curr_feat['DEL']
                if 'SNP' in curr_feat:
                    curr_snp_count = curr_feat['SNP']
                data_write = (str(feat) + "," + str(curr_mob_count) + "," + str(curr_ins_count) + "," +
                              str(curr_del_count) + "," + str(curr_snp_count) + "\n")
                of.write(data_write)
        x += 1
    return

def output_muttype_csv(new_map):

    total_count = 0
    snp_total = 0
    mob_total = 0
    ins_total = 0
    deletions_total = 0
    with open("output.csv", "a") as of:
        header= ",MOB,INS,DEL,SNP,TOTAL\n"
        of.write(header)
        for k in new_map.keys():
            data_lst = new_map[k].output_type_csv()
            print new_map[k].feat_type_map
            # Update totals
            mob_total += data_lst[0]
            ins_total += data_lst[1]
            deletions_total += data_lst[2]
            snp_total += data_lst[3]
            curr_total = data_lst[0] + data_lst[1] + data_lst[2] + data_lst[3]
            total_count += curr_total
            data_write = (str(k) + "," + str(data_lst[0]) + "," + str(data_lst[1]) + "," + str(data_lst[2]) +
            "," + str(data_lst[3]) + "," + str(curr_total) + "\n")
            of.write(data_write)
        totals_row = ("TOTALS," + str(mob_total) + "," + str(ins_total) + "," + str(deletions_total) + "," +
        str(snp_total) + "," + str(total_count) + "\n")
        
        of.write(totals_row)
   # Output mutations organized by feature
   

    return

if __name__ == "__main__":
    main()














