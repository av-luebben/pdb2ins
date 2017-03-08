head = """
    ########################################################################
    #                              PDB2INS                                 #
    #                  by Anna V. Luebben (Windows_AMD64@2017-03-06)       #
    ########################################################################

    Reads a PDB file and generates an .ins file for SHELXL.
    The PDB file is assumed to conform to the Protein Data Bank notes
    'Atomic Coordinate and Bibliographic Entry Format Description Version
    3.30'. For remarks and problems please contact aluebbe@gwdg.de.

    Usage:
        pdb2ins <filename|@pdbcode> [options]

    <filename|@pdbcode> Exact name of a file in PDB format or
                        a legal pdb code prefixed by '@'.

    Enter '--help' for complete options list.
    Enter 'q' or 'exit' to exit the program.

    """