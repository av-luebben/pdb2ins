__author__ = 'anna'
"""
This is an addition to pdb2ins.py, the pdb2ins program.
It is inspired by an excerpt from the program pymol (dali.py), fetching a pdb file after the pdb code is given.
Now changed to run with current RCSB PDB RESTful web interface service.
"""

# dali_file = "dali.txt"
#
# pdb_dir = "dali_pdb"
#
# max_pairs = 10


def fetchPDB(pdbCode, options, force=False):
    """
    Given a PDB code this function retrieves the PDB file from the PDB database directly.
    URL changed in Juni 2017 to current version, before :
    old url 'http://www.rcsb.org/pdb/cgi/export.cgi/' + remoteCode + '.pdb.gz?format=PDB&pdbId=' + remoteCode +
    '&compression=gz'
    Also pdb file is now named in the scheme pdbcode + _a.pdb to differentiate from shelxl pdb file.
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
        pdbFile = pdbCode + '_a.pdb'
        if os.path.exists(pdbFile):
            os.remove(pdbFile)
    #    pdbFile = pdbCode.lower() + '.pdb'
    remoteCode = string.upper(pdbCode)
    # if not os.path.exists(pdb_dir):
    #     os.mkdir(pdb_dir)

    if not os.path.exists(pdbFile) or force:  # new url: https://files.rcsb.org/download/4ZXB.pdb.gz
        try:
            filename = urllib.urlretrieve(
                'http://files.rcsb.org/download/' + remoteCode + '.pdb.gz')[0]
            # exit()
        except IOError:
            print 'ERROR: The RCSB PDB could not be reached.'
        # except Exception as ex:
        #     template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        #     message = template.format(type(ex).__name__, ex.args)
        #     print message
        # except:
        #     print "Warning: {} not found.\n".format(pdbCode)
        else:
            if os.path.getsize(filename) > 0:  # If 0, pdb code was invalid
                # print 'file is there', os.path.isfile(filename)
                try:
                    w = open(pdbFile, 'w')
                    # print open(filename).read()
                    w.write(gzip.open(filename).read())
                    w.close()
                    print "INFO: Fetched PDB file for code {} and saved as {}.".format(pdbCode, pdbFile)
                except IOError:
                    # print 'exception raised'
                    try:
                        os.remove(pdbFile)
                        print('ERROR: The PDB code you entered is not valid: {}'.format(pdbCode))
                        return None
                    except OSError:
                        print ' *** ERROR: Filename or pdb code not valid. *** '
                        exit()
            else:
                # print 'filesize', os.path.getsize(filename)
                print "WARNING: {} not valid. Download not successful. \n".format(pdbCode)
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
    # import gzip
    import os
    import string

    # print 'INFO: Fetching PDB redo file.'
    if options['o']:
        pdbFile = ''.join(str(options['o']).split('.')[:-1]) + '.pdb'
    else:
        pdbFile = pdbCode.lower() + '_a.pdb'
    remoteCode = string.lower(pdbCode)
    # if not os.path.exists(pdb_dir):
    #     os.mkdir(pdb_dir)
    if not os.path.exists(pdbFile) or force:
        try:
            filename = urllib.urlretrieve(
                'http://www.cmbi.ru.nl/pdb_redo/' +
                remoteCode[1:3] + '/' +
                remoteCode + '/' + remoteCode + '_besttls.pdb')[0]
        # except Exception as ex:
        #     template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        #     message = template.format(type(ex).__name__, ex.args)
        #     print message
        except IOError:
            print "WARNING: {} not found.\n".format(pdbCode)
        else:
            if os.path.getsize(filename) > 0:  # If 0, then pdb code was invalid
                try:
                    w = open(pdbFile, 'w')
                    w.write(open(filename).read())
                    w.close()
                    print "INFO: Fetched pdb file {} from PDB_REDO (best TLS) and saved as {}.".format(pdbCode, pdbFile)
                except IOError:
                    os.remove(pdbFile)
            else:
                print "WARNING: {} not valid.\n".format(pdbCode)
            os.remove(filename)
    return pdbFile