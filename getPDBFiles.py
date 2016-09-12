__author__ = 'anna'
"""
This is an addition to pdb2ins.py, the pdb2ins program.
It is inspired by an excerpt from the program pymol (dali.py), fetching a pdb file after the pdb code is given.

"""

# dali_file = "dali.txt"
#
# pdb_dir = "dali_pdb"
#
# max_pairs = 10


def fetchPDB(pdbCode, options, force=False):
    """
    Given a PDB code this function retrieves the PDB file from the PDB database directly.
    :param pdbCode: A valid pdb code
    :param force: If force is true an existing file with the 'pdbfile' name is overwritten by the new one.
    :return:
    """
    import urllib
    import gzip
    import os
    import string

    # print pdbCode
    if options['o']:
        pdbFile = ''.join(str(options['o']).split('.')[:-1]) + '.pdb'
    else:
        pdbFile = pdbCode + '.pdb'
    #    pdbFile = pdbCode.lower() + '.pdb'
    remoteCode = string.upper(pdbCode)
    # if not os.path.exists(pdb_dir):
    #     os.mkdir(pdb_dir)
    if not os.path.exists(pdbFile) or force:
        try:
            filename = urllib.urlretrieve(
                'http://www.rcsb.org/pdb/cgi/export.cgi/' +
                remoteCode + '.pdb.gz?format=PDB&pdbId=' +
                remoteCode + '&compression=gz')[0]
        except:
            print "Warning: {} not found.\n".format(pdbCode)
        else:
            if os.path.getsize(filename) > 0:  # If 0, pdb code was invalid
                try:
                    w = open(pdbFile, 'w')
                    w.write(gzip.open(filename).read())
                    w.close()
                    print "INFO: Fetched PDB file for: {}".format(pdbCode)
                except IOError:
                    try:
                        os.remove(pdbFile)
                    except OSError:
                        print ' *** ERROR: Filename or pdb code not valid. *** '
                        exit()
            else:
                print "Warning: {} not valid.\n".format(pdbCode)
            os.remove(filename)
    return pdbFile


def fetchPDBredo(pdbCode, options, force=False):
    """
    Given a PDB code this function retrieves the PDB file from the PDBredo database directly.
    :param pdbCode: A valid pdb code
    :param force: If force is true an existing file with the 'pdbfile' name is overwritten by the new one.
    :return: pdbfile
    """
    import urllib
    import gzip
    import os
    import string

    # print 'INFO: Fetching PDB redo file.'
    if options['o']:
        pdbFile = ''.join(str(options['o']).split('.')[:-1]) + '.pdb'
    else:
        pdbFile = pdbCode.lower() + '.pdb'
    remoteCode = string.lower(pdbCode)
    # if not os.path.exists(pdb_dir):
    #     os.mkdir(pdb_dir)
    if not os.path.exists(pdbFile) or force:
        try:
            filename = urllib.urlretrieve(
                'http://www.cmbi.ru.nl/pdb_redo/' +
                remoteCode[1:3] + '/' +
                remoteCode + '/' + remoteCode + '_besttls.pdb')[0]
        except:
            print "Warning: {} not found.\n".format(pdbCode)
        else:
            if os.path.getsize(filename) > 0:  # If 0, then pdb code was invalid
                try:
                    w = open(pdbFile, 'w')
                    w.write(open(filename).read())
                    w.close()
                    print "INFO: Fetched pdb file {} from PDB_REDO (best TLS).".format(pdbCode)
                except IOError:
                    os.remove(pdbFile)
            else:
                print "Warning: {} not valid.\n".format(pdbCode)
            os.remove(filename)
    return pdbFile