__author__ = 'tyler'

"""gd_sequence map: A class for storing information about mutations found in Genomediff files in predefined mappings.

This class is designed to simplify analysis and processing of Genomediff files. It stores the counts of various kinds
of mutations, organized in different dictionaries according to certain conditions; an example would be the type_map,
which organizes mutations by type.

This class is designed to work with the gd_sequence_mapper script. Details about that script can be found there or in
the documentation."""


class GenomeDiffSequenceMap:

    # More should be added to MUT_TYPES to accommodate all .gd files
    MUT_TYPES = ("INS", "DEL", "SNP", "MOB")

    # More may be added to accomodate additional features
    FEAT_TYPES = ("promoter", "RBS", "CDS", "misc_feature", "origin", "tactag", "tactagag")
    """__init__: instantiate all instance attributes"""
    def __init__(self):
        # Define mappings.
        self.type_map = dict()
        self.feat_map = dict()
        self.type_feat_map = dict()
        self.feat_type_map = dict()

        # Initialize all counts to 0.
        for mut_type in self.MUT_TYPES:
            self.type_map[mut_type] = 0
        for feat in self.FEAT_TYPES:
            self.feat_map[feat] = 0
        for mut_type in self.MUT_TYPES:
            for feat in self.FEAT_TYPES:
                self.type_feat_map[mut_type] = dict()
                self.type_feat_map[mut_type][feat] = 0
        for feat in self.FEAT_TYPES:
            for mut_type in self.MUT_TYPES:
                self.feat_type_map[feat] = dict()
                self.feat_type_map[feat][mut_type] = 0
        return

