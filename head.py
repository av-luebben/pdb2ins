head = """
    ########################################################################
    #                              PDB2INS                                 #
    #                  by Anna V. Luebben (Linux_x86_64@2016-09-08)        #
    ########################################################################

    Reads a PDB file and generates an .ins file for SHELXL.
    The PDB file is assumed to conform to the Protein Data Bank notes
    'Atomic Coordinate and Bibliographic Entry Format Description Version
    3.30'.

    Usage:
        pdb2ins <filename|@pdbcode> [options]

    <filename|@pdbcode> Exact name of a file in PDB format or
                        a legal pdb code prefixed by '@'.

    Enter '--help' for complete options list.
    Enter 'q' or 'exit' to exit the program.

    """