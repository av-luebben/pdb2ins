__author__ = 'anna'
"""
first project pdb2ins
by Anna Vera Luebben
start February 2015
version 2016/2 (October)

Read pdb file and generate .ins file for SHELXL.
The pdb file is assumed to conform to the Protein Data Bank notes
'Atomic Coordinate and Bibliographic Entry Format Description Version 3.30'.
"""
import os
import sys
import time
from sys import exit

import numpy as np

import pdb2hkl
import transformations
from LigandsInstructions import ligandRestraints
from ResiInstructions import instructions
from cmd import CommandlineParser
from head import head
from spagsydata import getSymmCards, testSpaceGroup

options = None
# Padding fills up the residue number to four digits (from the left) with zeros.
# padding = '0>4.0'
padding = '<4.0'

buildin_raw_Input = raw_input

######################for water model analysis########################
specAtom = 'OH'
specResi = 'TYR'
specResiFlag = 'wat'  # 'wat' or 'non'
# in the function 'writeSpecifiedAtomFile' an AtomName needed for gitty is specified. Needs to be changed for other Resi
######################################################################


def dummy(_):
    """
    Dummy
    :param _: Dummy
    :return: None
    """
    pass


class OutputProxy(object):
    """
    Object replacing sys.stdout to reroute output from stdout to a callback.
    """
    def __init__(self, cb):
        self.cb = cb

    def write(self, *args):
        """
        Interface to the 'print' statement.
        :param args: list of strings
        :return: None
        """
        self.cb(args)


def setSlaveMode(opt, cb=None):
    """
    Sets pdb2ins script to SlaveMode thereby using an dictionary object created by cmd.CommanineParser.__call__()
    instead of using sys.argv to build the dictionary itself and intercepting all output to sys.stdout and rerouting
    it to the callback object 'cb'
    :param opt: Dictionary type defining all options as output by cmd.CommandlineParser.__call__().
    :param cb: Callable that gets called with the string arguments that are usually written to sys.stdout.
    :return: None
    """
    if not cb:
        cb = dummy
    global options
    options = opt
    sys.stdout = OutputProxy(cb)


def raw_input(*args, **kwargs):
    """
    Enables the user to terminate the program by typing either 'q' or 'exit' during raw_input.
    :param args:
    :param kwargs:
    :return:
    """
    inputString = buildin_raw_Input(*args, **kwargs)
    if inputString.lower() == 'q' or inputString.lower() == 'exit':
        print '*** PDB2INS has been terminated ***'
        exit()
    return inputString


class Data(object):
    def __init__(self):
        self.strings = None
        self.IO = IO(options)
        self.hklf = None
        self.hklfile = None
        self.atomContainer = AtomContainer()
        self.header = Header()
        self.askHKL()
        # print 'now starting makeHKLfile()'
        self.makeHKLfile()
        # print 'finished with making hkl file.'
        self.IO.read()
        self.hasHAtoms = False
        self.neutronData = False
        self.readContent()
        self.printWarnings()
        self.header.extractWavelength()
        self.atomContainer.extractAllElements()
        self.dealWithHAtoms()
        self.atomContainer.getResidueList()
        # self.header.extractCell()
        if specAtom and specResi:
            # print 'spec atom found'
            self.atomContainer.findSpecifiedAtom(specResi, specAtom, None)
            self.atomContainer.writeSpecifiedAtomFile()
            self.atomContainer.findAllNonAAAtoms()
            self.atomContainer.findAllWaterAtoms()
            self.atomContainer.findNonAANearSpecAtom(specResiFlag)
            self.atomContainer.makeOMITrecords()
        self.header.abbreviateSpaceGroup()
        self.atomContainer.createSequenceDict()
        self.header.handleResidueSequence(self.atomContainer.getSequenceDict())
        # self.header.extractScale()
        self.header.makeGeneralRefinementInstructions()
        if 'HOH' in self.atomContainer.getOtherResiSet():
            self.askWaterOccupancy()

        self.joinstrings()
        self.IO.writeFile(self)

    def askHKL(self):
        """
        The first question the user is asked after starting pdb2ins. (pdb2ins will need the format of the hkl file
        therefore it is necessary to run pdb2hkl subroutine first.)
        'yes' will lead down the pdb2hkl path. Default answer is 'no'
        :return:
        """
        if not options['i'] and not options['b']:
            while True:
                doHKL = raw_input('\nCreate .hkl file from structure factor file (cif) or PDB code? (y or n) '
                                  '[N]: ')
                if not doHKL:
                    break
                if doHKL == 'Y' or doHKL == 'y':
                    options['b'] = True
                    break
                if doHKL == 'N' or doHKL == 'n':
                    break
        if options['b'] and not options['filename']:  # check if indentation correct!
            self.askHKLfilename()

    def askHKLfilename(self):
        """
        Choose a file or enter a pdbcode prefixed with '@'. The prefix will signal the program to download the -sf.cif
        file from the RCSB PDB.
        When the file should be downloaded, the function checks if the file is already available locally.
        :return:
        """
        while True:
                self.hklfile = raw_input("\nEnter name of structure factor file to read (To download a PDB "
                                         "file enter \'@<PDBCODE>\'): ")#.upper()
                if not os.path.isfile(self.hklfile) and not self.hklfile.startswith('@'):
                    newstring = str(self.hklfile[:-4].upper())+str(self.hklfile[-4:])
                    if not os.path.isfile(newstring) and not os.path.isfile(self.hklfile.lower()):
                        # print self.hklfile.upper(), newstring, self.hklfile.lower()
                        print 'INFO: File \'{}\' not found.'.format(self.hklfile)
                    if os.path.isfile(newstring):
                        self.hklfile = newstring
                        break
                    if os.path.isfile(self.hklfile.lower()):
                        self.hklfile = self.hklfile.lower()
                        break
                    if not self.hklfile.endswith('.pdb'):
                        self.hklfile += '.pdb'
                        if os.path.isfile(self.hklfile):
                            print 'INFO: Using file \'{}\' instead.'.format(self.hklfile)
                            break
                if self.hklfile.startswith('@'):
                    if len(self.hklfile) == 5:
                        options['filename'] = self.hklfile
                        break
                    else:
                        print 'ERROR: Given pdb code is not correct! Please check.'
                        pass
                else:
                    break
        options['d'] = self.hklfile

    def makeHKLfile(self):
        """
        options 'b':tells us a hkl file should be created.
        options 'd': a input file in .cif format containing structure factors is given to create an .hkl
        options 'o': an output filename for the .ins file is specified and will be used for the .hkl file also.
        options 'filename': if starting with an '@', this code is used to download the *-sf.cif file from the RCSB PDB.
        This function tries to generate a input and output filename for the structure factor files.
        Those and other information is given as optionsForPdb2hkl to the subroutine pdb2hkl.
        The subroutine pdb2hkl is started.
        :return:
        """
        # print ' first line makeHKLfile.'
        # print 'these are the Options', options
        if options['b']:
            # print 'options b found.'
            if options['d']:
                filename = options['d']
            else:
                filename = options['filename']
            outfile = None
            if options['o']:
                insOutfilename = str(options['o'])
                outfile = ''.join(str(insOutfilename).split('.')[:-1]) + '.hkl'
                # print '2', outfile
            elif str(filename).startswith('@'):
                insfilename = options['filename']
                # print '1', insfilename
                # if insfilename:
                if not str(insfilename).startswith('@'):
                    outfile = ''.join(str(insfilename).split('.')[:-1]) + '.hkl'
                else:
                    outfile = str(insfilename[1:]) + '.hkl'
                    # print '3', outfile
                # else:
                #     outfile = str(insfilename[1:]) + '.hkl'
                # print '4', outfile
            if not outfile:  # if not str(filename).startswith('@'):
                if '-sf' in filename:
                    outfile = ''.join(str(filename).split('-')[:-1]) + '.hkl'
                else:
                    outfile = ''.join(str(filename).split('.')[:-1]) + '.hkl'  # "".join(f.split('.')[:-1]) + '.ins'
            # print '5', outfile
            if options['i']:
                i = True
            else:
                i = False
            # print outfile
            optionsForPdb2hkl = {'filename': filename, 'i': i, 'o': outfile}
            # print 'these are the options for pdb2hkl', optionsForPdb2hkl
            try:
                print 'INFO: Starting pdb2hkl.'
                try:
                    self.header.hklf = pdb2hkl.run(optionsForPdb2hkl)
                except SystemExit:
                    # print 'error occured with options 'b', System Exit! '
                    pass
                # self.hklf = pdb2hkl.Data.getHKLF()
                # print self.header.hklf
            except SystemError:
                pass

    def getHKLF(self):
        return self.header.hklf

    def buildInstructions(self):
        """
        The list of all residues present in the pdb file is used to generated the general refinement instructions.
        First it is checked if HOH (water) is present and a CONN, ISOR line is added.
        Next RTAB, HFIX and restraints are fetched for all natural amino acids.
        The necessary instructions for aa are given in the ResiInstructions file.
        For a selection of ligands the restraints can be added from LigandInstructions subroutine.
        At last all residues without restraints (i.e. ligands, metals) are listed.
        :return: string
        """
        self.atomContainer.makeRestraintsForTermini(self.header.getResiDict())
        if 'HOH' in self.atomContainer.getOtherResiSet():
            waterInstructions = ['ISOR_HOH 0.1 $O  !water atoms are restraint to near isotropic behavior',
                                 'CONN_HOH 0 O  !generation of connectivity table fine-tunning']
            # waterInstructions = ['ISOR_HOH 0.1 $O', 'CONN_HOH 0 $O']
        else:
            waterInstructions = ["REM ISOR and CONN 0 recommended on adding water"]

        blockNames = ['RTAB', 'HFIX', 'Restraints']
        missingRestraints = set()
        instructionStrings = [waterInstructions]
        instructionStrings.append(['\nREM Restraints and HFIX for terminal residues as follows:\n'])
        instructionStrings.append(self.atomContainer.getRestraintsForTermini())
        for i, instruction in enumerate(instructions):
            instructionStrings.append(['\nREM {} instruction block\n'.format(blockNames[i])])
            for resiName in ['All']+self.atomContainer.getResidueList() + self.atomContainer.getOtherResiSet():
                try:
                    resiInstruction = instruction[resiName]
                    instructionStrings.append(resiInstruction)
                except KeyError:
                    if i == 2 and resiName not in ['HOH', 'All']:
                        missingRestraints.add(resiName)
                    pass
        foundRestraints = set()
        for i in missingRestraints:
            try:
                ligandRestraintsFound = ligandRestraints[i]
                instructionStrings.append(['\n'])
                instructionStrings.append(['REM Restraints for ligand {}:\n'.format(i)])
                instructionStrings.append(ligandRestraintsFound)
                foundRestraints.add(i)
                instructionStrings.append(['\n'])
            except KeyError:
                continue
        missingRestraints -= foundRestraints
        if missingRestraints:
            if len(missingRestraints) <= 4:
                instructionStrings.append(['\nREM Restraints missing for the following residues: ' + ", ".join(
                    [str(i) for i in list(missingRestraints)])])
            else:
                instructionStrings.append(['\nREM Restraints missing for the following residues: \n',
                                           'REM ' + ', ' .join([str(i) for i in list(missingRestraints)])])
            instructionStrings.append(['\n'])
            text = '\nINFO: Following ligands/residues ' \
                   'have no restraints:\n      ' + ', ' .join([str(i) for i in list(missingRestraints)])
            print text + '\n      Please remember to manually add restraints for this residues.\n'

        import itertools
        return '\n'.join(list(itertools.chain.from_iterable(instructionStrings)))

    def readContent(self):
        """
        reads the first 6 letters in the pdb file and sorts the line into the appropriate class for further use.
        The line ANISOU is only given to atom container if useAnisou is TRUE.
        :return: line
        """
        useAnisou = True
        alreadyasked = False
        # try:
        for line in self.IO.dataf:
            if line[0] == "#":
                continue
            if line[:6] == 'EXPDTA':
                self.isCrystData(line)
            if line[:6] == "ATOM  " or line[:6] == "HETATM":
                self.atomContainer.extractAtom(line)
            if line[:6] == "HET   ":
                self.atomContainer.extractHetRecord(line)
            if line[:6] == "ANISOU" and useAnisou:
                if not alreadyasked:
                    useAnisou = self.askAnisou()
                    alreadyasked = True
                if useAnisou:
                    self.atomContainer.extractAtomAnisou(line)
            else:
                self.header.interpretLine(line)
        # except TypeError:
        #     print line
        #     print '\nERROR: File is not a PDB file.\n *** PDB2INS is terminated without writing an .ins file. ***'
        #     exit()

    def printWarnings(self):
        """
        All Warnings that can occur multiple times during reading the ATOM lines of the pdb file should only be printed
        once. Therefore this function is called after reading all lines and prints the messages.
        :return:
        """
        # if self.atomContainer.negResiNumber:
        #     print '\n*** WARNING: Negative residue numbers found in file. SHELXL might not be able to handle.***\n'
        if self.atomContainer.overlongResiNum:
            print '\n*** WARNING: One or more residues have a residue number larger than 10 000. Please check!***\n'
        if self.atomContainer.wrongResiName:
            print '\n*** WARNING: The files contain residue names starting with a number. \n' \
                  '    Older SHELXL versions can only handle residues names starting with a letter.\n'
        if self.atomContainer.waterOffsetWarning:
            print '\nINFO: Water residue numbers were changed to handle residues with insertion codes.\n'
        if self.atomContainer.resiNUmberCollision:
            print '*** Warning: The residues with a number larger than 1000 could collide with applied offset \n' \
                  '    for insertion code residues. Please check!\n'

    def dealWithHAtoms(self):
        """
        This functions is called to find out whether X-Ray diffraction data was given with H atoms (eg. from PHENIX pdb
        files). When H atoms are found in the pdb file, the function askAboutHAtoms is called.
        :return:
        """
        # print self.atomContainer.elementDict.keys()
        # for key in self.atomContainer.elementDict.keys():
        #     print key
        if self.atomContainer.elementDict['H'] > 0 and not self.neutronData:
            self.hasHAtoms = True
        if self.hasHAtoms:
            self.askAboutHAtoms()

    def askAboutHAtoms(self):
        """
        This function is called when H atoms where found in an pdb file containing x-ray diffraction data.
        The user is informed that shelxl can use the HFIX command to generate H atoms and it is recommended to erase
        them.
        It should be noted that pdb files containing some of their natural aa residues in PART instruction will have
        trouble with HFIX in shelxl.
        :return:
        """
        print ('\nINFO: This pdb file contains Hydrogen atoms from X-ray diffraction data. \n'
               'It is recommended to delete all Hydrogen atoms now and use HFIX in shelxl \n'
               'to place them again. This program automatically creates the necessary \n'
               'HFIX instructions for natural amino acids.')
               # 'It should be noted that is not recommended to use the HFIX instructions after disorder has been
               # modeled')
        if not options['e']:
            if not options['i']:
                reply = raw_input('Delete all Hydrogen atoms in .ins file? (y or n) [Y]: ')
                if reply == 'N' or reply == 'n':
                    self.atomContainer.keepHAtoms = True
                elif reply == 'Y' or reply == 'y' or not reply:
                    self.atomContainer.keepHAtoms = False
                else:
                    self.atomContainer.keepHAtoms = False
            else:
                print 'INFO: Hydrogen atoms will not be transferred to .ins file.'
                self.atomContainer.keepHAtoms = False
        else:  # This part should only be run if options 'e' is True.
            self.atomContainer.keepHAtoms = True

    def isCrystData(self, line):
        """
        takes the pdb line starting with EXPDTA as an input (which specifies the experiment) and searches for
        'X-RAY DIFFRACTION' in the string. if the string does not contain this phrase, the user is given the line and
        asked whether the program should continue. Default answer is 'no', which will terminate the program.
        :param line: string
        :return:
        """
        self.neutronData = False
        if 'NEUTRON' in line:
            self.neutronData = True
            self.atomContainer.neut = True
            print '\nINFO: This file contains NEUTRON diffraction data. ' \
                  'The .ins file will now contain the NEUT instruction.\n' \
                  'Please be aware that the restraints for ligands may not be suitable for neutron data.'
        if 'X-RAY DIFFRACTION' not in line and 'NEUTRON' not in line:
            print 'This pdb file contains the following experimental data:'
            line1 = line[6:]
            print line1
            reply = raw_input('\nThis pdb file might not contain X-RAY diffraction data.\n'
                              'Important information necessary to create an .ins file might be missing.\n'
                              'Missing data can cause the program to terminate inadvertently.\n'
                              'Continue anyway? (y or n) [N]: ')
            # print reply
            if reply == 'N' or reply == 'n' or not reply:
                print ' *** PDB2INS is terminated without writing an .ins file. ***'
                exit()

    def askAnisou(self):
        """
        -only called when anisotropic data is present.
        Checks options for an entry under 'i'. If the option is given, the boolean useAnisou is set to TRUE.
        If no option was given (interactive modus), the user is prompted whether the anisotropic data should be
        converted to isotropic data.
        default answer is 'yes'. 'yes' sets the useAnisou to FALSE, 'no' to TRUE.
        :return: useAnisou (boolean) or options['a']
        """
        if options['a']:
            return options['a']
        else:
            if not options['i']:
                reply = raw_input('\nThe pdb file contains anisotropic atom data. \nConvert anisotropic atoms to '
                                  'isotropic? (y or n) [Y]: ')
                if reply == 'Y' or reply == 'y' or not reply:
                    useAnisou = False
                elif reply == 'N' or reply == 'n':
                    useAnisou = True
                else:
                    useAnisou = False
                    print " *** ERROR: Invalid response \'{}\'. " \
                          "Anisotropic atoms will be converted to isotropic as default. *** ".format(reply)
                return useAnisou
            else:
                useAnisou = options['a']
                if useAnisou:
                    answer = 'No'
                else:
                    answer = 'Yes'
                print 'INFO: The pdb file contains anisotropic data. Convert to isotropic? (y or n) [Y]: ' \
                      '\n     PDB2INS used default answer "{}".'.format(answer)
                return options['a']

    def askWaterOccupancy(self):
        """
        if the interactive mode is True, the user is asked whether the water occupancy should be reset to unity.
        default is to reset occupancy by calling the resetOccupancy function in class Atom.
        :return:
        """
        if not options['i']:
            reply = raw_input("\nReset water occupancy to unity? (y or n) [Y]: ")
            if not reply or reply == 'y' or reply == 'Y':
                resetOccupancy = True
            elif reply == 'n' or reply == 'N':
                resetOccupancy = False
            else:
                print 'Invalid response. Water occupancy reset to unity as default.'
                resetOccupancy = True
        else:
            print 'INFO: PDB2INS reset water occupancy to unity as default.'
            resetOccupancy = True
        if resetOccupancy:
            Atom.resetOccupancy('HOH', 1)

    def joinstrings(self):
        """
        creates strings for the lines needed in the header of .ins files.
        appends all strings created in atom container.
        :return: a joined string of all strings needed for .ins files
        """
        self.strings = []
        self.strings.append("TITL converted from file {}\n".format(self.IO.workfile))
        self.strings.append("CELL {:7.5f} {} \n".format(self.header.getWavelength(),
                                                        '{:7.3f} {:7.3f} {:7.3f} {:7.2f} {:7.2f} {:7.2f}'
                                                        .format(*self.header.getCell())))
        self.strings.append("ZERR      {} {} \n\n".format(self.header.getZvalue(),
                                                          '{:7.3f} {:7.3f} {:7.3f} {:7.2f} {:7.2f} {:7.2f}'
                                                          .format(*self.header.getCellError())))
        self.strings.append("REM Space group {} \n\n".format(self.header.getSpaceGroup()))
        self.strings.append("LATT {} \n".format(self.header.getLattice()))
        #self.header.validateSpaceGroup(self.header.getAbbrSpaceGroup())
        self.strings += getSymmCards(self.header.getAbbrSpaceGroup())
        # try:
        #     self.strings += getSymmCards(self.header.getAbbrSpaceGroup())
        # except KeyError:
        #     print 'Space group {} is not valid.'.format(self.header.getAbbrSpaceGroup())
        #     exit()
        if self.neutronData:
            self.strings.append("NEUT \n")
        self.strings.append("SFAC  {} \n".format(('{} '*len(self.atomContainer.getElementList()))
                                                 .format(*self.atomContainer.getElementList())))
        if self.header.makeDISP:
            self.header.makeDISPinstruction(self.atomContainer.getElementList())
            self.strings += self.header.getDISPinstructions()
        self.strings.append("UNIT  {} \n\n".format(('{} '*len(self.atomContainer.getElementList()))
                                                   .format(*[self.atomContainer.elementDict[key] * self.header.getZvalue()
                                                             for key in self.atomContainer.getElementList()])))
        atomContainerString = self.atomContainer.asShelxString(self.header.getCell())
        self.strings.append(self.header.getGeneralRefinementInstructions())
        self.strings.append("\n")
        self.strings.append(self.header.getRIGUInstructions(self.atomContainer.getElementList2()))
        self.strings.append(self.header.getDELUInstruction(self.atomContainer.getElementList2()))
        self.strings.append(self.header.getSIMUInstruction(self.atomContainer.getElementList2()))
        self.strings.append(self.header.getGeneralRefinementInstructions3())
        if self.atomContainer.incompleteResiString:
            self.strings.append("\n")
            self.strings.append("REM HFIX 0 instructions were added for incomplete residues.\n")
            self.strings += self.atomContainer.incompleteResiString
        self.strings.append("\n")
        self.strings.append(self.buildInstructions())
        self.strings.append("\n")
        if self.atomContainer.insertionCodeString:
            self.strings.append('REM Instructions for residues with insertion code.\n')
            self.strings.append("\n")
            for i in self.atomContainer.insertionCodeString:
                self.strings.append(i + '\n')
            self.strings.append("\n")
        self.atomContainer.findSSBonds()
        if self.atomContainer.ssBonds:
            self.strings.append("REM Instructions for disulfide-bridges:\n")
            self.strings += self.atomContainer.getSSBonds()
            self.strings.append("\n\n")
        if self.atomContainer.omitList:
            self.strings.append("REM Instructions for OMIT atoms: \n")
            self.strings += self.atomContainer.getOmitAtoms()
        self.strings.append(self.atomContainer.asShelxString(self.header.getCell()))
        # self.strings.append(self.atomContainer.asShelxString(self.header.getCell()))
        self.strings.append("\n\n")
        self.strings.append(self.header.getHklf()+"\nEND")


class IO(object):

    def __init__(self, options, *args, **kwargs):
        self.workfile = None
        self.dataf = None
        self.outputFilename = None
        self.usePDBredo = False
        self.options = options

    def askPDBredo(self):
        """
        User is promted whether the RCSB or PDB-REDO server should be used to download the .pdb file. Options are
        enumerated 1 and 2.
        :return:
        """
        # if not self.options['GUI']:
        if not self.options['i'] and not self.options['r']:
            while True:
                pdbRedo = raw_input('\nDownload PDB file from RCSB Protein Data Base (1) or PDB_REDO database '
                                    '(2)? [1]: ')
                if not pdbRedo or pdbRedo == '1':
                    self.usePDBredo = False
                    break
                elif pdbRedo == '2':
                    self.usePDBredo = True
                    break
                else:
                    pass
        if self.options['i'] and not self.options['r']:
            self.usePDBredo = False

    def askFilename(self):
        """
        Asks the user for a name of the pdb file that needs to be converted in the form 'name.format'.
        It is also possible to give the PDB code of the desired file in the form '@PDBCODE'.
        In this case the function fetch_pdb in fromDali is called upon to retrieve the code from the pdb database.
        Some pdb files to main:
        3LOH (big file!)
        1AZI (relatively small)
        3YMI, 1CTJ (space group R 3)
        2QXX (space group R 3 2)
        :return: self.workfile
        """
        self.workfile = self.options['filename']
        # print self.workfile, 'this is the workfile'
        if self.workfile:  # when cmd parser is used, the filename should already be specified, skipping user input
            try:
                # print '1'
                if not os.path.isfile(self.workfile) and '@' not in self.workfile:
                    # print '2'
                    newstring = str(self.workfile[:-4].upper())+str(self.workfile[-4:])
                    if not os.path.isfile(newstring) and not os.path.isfile(self.workfile.lower()):
                        # print '3'
                        print 'INFO: File {} not found.'.format(self.workfile)
                        if not self.workfile.endswith('.pdb'):
                            self.workfile += '.pdb'
                            if os.path.isfile(self.workfile):
                                print 'INFO: Using file \'{}\' instead.'.format(self.workfile)
                            else:
                                print ' *** Error: Given file name not valid. *** '
                                self.workfile = None
                    if os.path.isfile(newstring):
                        # print '4'
                        self.workfile = newstring
                    if os.path.isfile(self.workfile.lower()):
                        # print '5'
                        self.workfile = self.workfile.lower()
                # if self.workfile.startswith('@'):
                #     trystring = str(self.workfile[-4:]) + '_a.pdb'  # if the file was already loaded in the GUI
                #     if os.path.isfile(trystring):
                #         print trystring
                #         self.workfile = trystring
            except TypeError:
                self.workfile = None
        if not self.workfile:  # in interactive mode without cmd options, the user is asked for the filename
            if not self.options['d']:
                while True:
                    self.workfile = raw_input("\nEnter name of PDB file to read (To download a PDB "
                                              "file enter \'@<PDBCODE>\'): ")  # .upper()
                    if not os.path.isfile(self.workfile) and not self.workfile.startswith('@'):
                        newstring = str(self.workfile[:-4].upper())+str(self.workfile[-4:])
                        if not os.path.isfile(newstring) and not os.path.isfile(self.workfile.lower()):
                            # print self.workfile.upper(), newstring, self.workfile.lower()
                            print 'INFO: File \'{}\' not found.'.format(self.workfile)
                        if os.path.isfile(newstring):
                            self.workfile = newstring
                            break
                        if os.path.isfile(self.workfile.lower()):
                            self.workfile = self.workfile.lower()
                            break
                        if not self.workfile.endswith('.pdb'):
                            self.workfile += '.pdb'
                            if os.path.isfile(self.workfile):
                                print 'INFO: Using file \'{}\' instead.'.format(self.workfile)
                                break
                    else:
                        break
            else:  # here the possibility is handled, that the user is in interactive mode and created a .hkl already
                hklfilename = options['d']  # the filename of the sf file is taken and an input filename suggested
                if hklfilename.startswith('@'):
                    possiblePdbFilename = str(hklfilename)
                # else:
                #     possiblePdbFilename = ''.join(str(hklfilename).split('.')[:-1]) + '_a.pdb'
                while True:
                    self.workfile = raw_input("\nEnter name of PDB file to read (To download a PDB "
                                              "file enter \'@<PDBCODE>\')[{}]: ".format(possiblePdbFilename))  #.upper()
                    if not self.workfile:
                        self.workfile = possiblePdbFilename
                    if not os.path.isfile(self.workfile) and not self.workfile.startswith('@'):
                        newstring = str(self.workfile[:-4].upper())+str(self.workfile[-4:])
                        if not os.path.isfile(newstring) and not os.path.isfile(self.workfile.lower()):
                            # print self.workfile.upper(), newstring, self.workfile.lower()
                            print 'INFO: File \'{}\' not found.'.format(self.workfile)
                        if os.path.isfile(newstring):
                            self.workfile = newstring
                            break
                        if os.path.isfile(self.workfile.lower()):
                            self.workfile = self.workfile.lower()
                            break
                        if not self.workfile.endswith('.pdb'):
                            self.workfile += '.pdb'
                            if os.path.isfile(self.workfile):
                                print 'INFO: Using file \'{}\' instead.'.format(self.workfile)
                                break
                    else:
                        break
        else:
            # print '6'
            self.workfile = self.workfile
        if self.workfile.startswith('@'):  # if the user was asked for a filename, it is transferred to options
            # print '7'
            if not self.options['filename']:
                self.options['filename'] = self.workfile  # now a correct output filename can be created
            if self.options['r']:
                self.usePDBredo = True
            else:
                self.askPDBredo()
            if self.usePDBredo:
                self.options['r'] = True
        if self.workfile.startswith('@') and self.options['r']:
            from getPDBFiles import fetchPDBredo
            self.workfile = fetchPDBredo(self.workfile[1:], self.options)
        elif self.workfile.startswith('@'):  # here elif when the if statement before is used!
            # print '8'
            from getPDBFiles import fetchPDB
            # print self.workfile[1:]
            self.workfile = fetchPDB(self.workfile[1:], self.options)
            # print "INFO: Fetching PDB file for entry {}.".format(self.workfile)
        else:
            pass
            # self.workfile = '3LOH.pdb'  # nur zu Testzwecken

    def read(self):
        """
        Takes the pdb file given from askFilename and opens it.
        :return: Lines from pdb file as data.f
        """
        self.askFilename()
        try:
            workf = open(self.workfile, 'r')
            self.dataf = workf.readlines()
            # print "reading workfile \'{}\'.".format(self.workfile)
            print 'INFO: File {} successfully opened.'.format(self.workfile)
        except IOError:
            self.askFilename()

    def askOutputFilename(self):
        """
        Asks the user for a name for the output file .ins.
        :return: output filename
        """
        # defaultName = os.path.splitext(self.workfile)[0] + '.ins'
        try:
            if '_' in self.workfile and self.options['filename'].startswith('@'):
                defaultName = os.path.splitext(self.workfile)[0].split('_')[0] + '.ins'
            else:
                defaultName = os.path.splitext(self.workfile)[0] + '.ins'
        except AttributeError:
            if '_' in self.workfile:
                defaultName = os.path.splitext(self.workfile)[0].split('_')[0] + '.ins'
            else:
                defaultName = os.path.splitext(self.workfile)[0] + '.ins'
        if not self.options['i']:
            self.outputFilename = raw_input("\nEnter .ins filename to be created [{}]: ".format(defaultName))
            if not self.outputFilename:
                self.outputFilename = defaultName
            elif not self.outputFilename.endswith('.ins'):
                self.outputFilename += '.ins'
        else:
            self.outputFilename = self.options['o']
            if not self.outputFilename:
                self.outputFilename = defaultName
            elif not self.outputFilename.endswith('.ins'):
                self.outputFilename += '.ins'

    def writeFile(self, data):
        """
        writes the joined strings from DATA in an output file either named by the user in askOutputFilename or calls it
        newfile.ins.
        :param data: Gives joined string of all data needed for the .ins file
        :return: A new file in .ins def getTerminalResidues(self):format.
        """
        self.askOutputFilename()
        if not self.outputFilename:
            print "INFO: File newfile.ins created."
            f = open("newfile.ins", 'w')
            f.write(''.join(data.strings))
        else:
            print "INFO: File {} created.".format(self.outputFilename)
            f = open("{}".format(self.outputFilename), 'w')
            f.write(''.join(data.strings))
        f.close()
        global start_time
        t = (time.time() - start_time)
        print 'INFO: +++ PDB2INS finished running after {:3.4f} seconds. +++ '.format(t)


class Header(object):

    def __init__(self):
        self.cell = None
        self.crystLine = None
        self.cellError = []
        self.remark200Lines = []
        self.scaleLine = []
        self.wavelength = None
        self.lattice = None
        self.sequenceLines = []
        self.ssBondLines = []
        self.resiDict = {}
        self.generalRefinement = []
        self.generalRefinement2 = []
        self.generalRefinement3 = []
        self.generalRefinement2Order = None
        self.ssBondList = []
        self.userCell = None
        self.hklf = None
        self.zValue = None
        self.spaceGroup = None
        self.shortenedSpaceGroup = None
        self.kisselfile = None
        self.kisselContent = None
        self.dispInstructions = []
        self.makeDISP = False
        self.sequenceDict = {}

    def interpretLine(self, line):
        """
        Takes all lines from the input file and sorts the according to their first six letters.
        :param line: From the input PDB file.
        :return: crystLine to extract all crystallographic data and remark200Line to extract wavelength.
        """
        if line[:6] == "CRYST1":
            # print line
            self.crystLine = line
            self.askHklf()
            self.extractCell()
            self.extractSpaceGroup()
            self.extractZvalue()
            self.estimateCellError()
        if line[:10] == "REMARK 200":
            self.remark200Lines.append(line)
        if line[:6] == "SEQRES":
            self.sequenceLines.append(line)
        if line[:5] == "SCALE":
            self.scaleLine.append(line)
        if line[:6] == 'NUMMDL':
            self.extractNumberOfModels(line)

    def extractNumberOfModels(self, line):
        """
        If the line of the PDB file starts with NUMMDL, the number of models in the file should follow.
        Is this number greater than one. pdb2ins will terminate. Pdb2ins as well as shelxl cannot handle multiple
        models.
        :param line:
        :return:
        """
        l = line[6:].strip()
        try:
            numberOfModels = int(l)
        except ValueError:
            numberOfModels = 2
            pass
        if numberOfModels > 1:
            print 'ERROR: This PDB file contains more than one model. PDB2INS cannot handle multiple models.\n' \
                  '*** PDB2INS is terminated without writing an .ins file. ***'
            exit()

    def extractCell(self):
        """
        If the entry describes a structure determined by a technique other than X-ray crystallography,
         CRYST1 contains a = b = c = 1.0, alpha = beta = gamma = 90 degrees, space group = P 1, and Z = 1.
         """
        if options['c']:
            try:
                arg = options['c'].split(',')
                arg = [float(i) for i in arg if i]
                cell_a = arg[0]
                cell_b = arg[1]
                cell_c = arg[2]
                alpha = arg[3]
                beta = arg[4]
                gamma = arg[5]
                self.cell = [cell_a, cell_b, cell_c, alpha, beta, gamma]
                print 'INFO: Cell is set to {}.'.format(self.cell)
            except ValueError:
                self.cell = None
            except IndexError:
                self.cell = None
        if not self.cell:
            try:
                cell_a = float(self.crystLine[6:15])
                cell_b = float(self.crystLine[15:24])
                cell_c = float(self.crystLine[24:33])
                alpha = float(self.crystLine[33:40])
                beta = float(self.crystLine[40:47])
                gamma = float(self.crystLine[47:54])
                self.cell = [cell_a, cell_b, cell_c, alpha, beta, gamma]
                # print self.cell
                if options['c']:
                    print 'INFO: Cell is set to {}.'.format(self.cell)
            except ValueError:
                # print cell_a, cell_b, cell_c
                self.cell = None
        if self.cell:
            if (cell_a + cell_b + cell_c) < 1:
                self.cell = None
            if cell_a <= 15.00 or not 20 < alpha < 160:
                print "\nINFO: Warning: Cell may not be correct! Please check."
                # print self.cell, float(cell_a)
                #if options['i']:
                #    print ' ** Error: No cell found. ** '
                #    exit()
                reply = raw_input("Cell found: ({a} {b} {c} {alpha} {beta} {gamma}). This may not be correct. \n"
                                  "Please enter correct cell [{a} {b} {c} "
                                  "{alpha} {beta} {gamma}]:".format(a=float(cell_a), b=float(cell_b), c=float(cell_c),
                                                                    alpha=float(alpha), beta=float(beta),
                                                                    gamma=float(gamma)))
                if not reply:
                    pass
                else:
                    try:
                        arg = reply.split(' ')
                        arg = [float(i) for i in arg if i]
                        cell_a = arg[0]
                        cell_b = arg[1]
                        cell_c = arg[2]
                        alpha = arg[3]
                        beta = arg[4]
                        gamma = arg[5]
                        self.cell = [cell_a, cell_b, cell_c, alpha, beta, gamma]
                    except IndexError:
                        self.cell = None
                        print 'Cell given is not valid: ', reply
                    except TypeError:
                        self.cell = None
        else:
            while not self.cell:
                self.askCell()

    def askCell(self):
        if not options['c']:
            reply = raw_input('Enter cell (a b c alpha beta gamma):')
            try:
                arg = reply.split(' ')
                arg = [float(i) for i in arg if i]
                cell_a = arg[0]
                cell_b = arg[1]
                cell_c = arg[2]
                alpha = arg[3]
                beta = arg[4]
                gamma = arg[5]
            except IndexError:
                return False
            except TypeError:
                return False
            if cell_a <= 2.00 or not 20 < alpha < 160:
                print "INFO: Warning: Cell may not be correct! Please check."
                #if options['i']:
                #    print ' ** Error: No cell found. ** '
                #    exit()
            self.cell = [cell_a, cell_b, cell_c, alpha, beta, gamma]
            # print self.cell
        else:
            self.cell = options['c']
            # print 'INFO: cell is set to {}.'.format(self.cell)

    def askHklf(self):
        """
        Sets the HKLF code. Default value is 4.
        :return:
        """
        # print options['i'], options['h']
        # if options['b']:
        #     self.hklf = Data.getHKLF()
        if not self.hklf:
            if not options['i'] and not options['h']:
                reply = raw_input('\nEnter HKLF code (3 for F, 4 for F-squared) [4]:')
                if not reply:
                    self.hklf = '4'
                elif reply == '4' or reply == '3':
                    self.hklf = reply
            elif options['h']:
                self.hklf = options['h']
            else:
                self.hklf = 3
                print 'INFO: HKLF is set to default: HKLF{}. Please check if HKLF is correct.'.format(self.hklf)
        # else:
        #    self.hklf = '4'
        #    print " ** ERROR: Input HKLF code {} not valid. " \
        #          "HKLF code was set to value '4' for F-squared. ** ".format(reply)

    def getHklf(self):
        return 'HKLF {}\n'.format(self.hklf)

    def getCell(self):
        # self.userCell = options['c']
        # if self.userCell:
        #     return self.userCell
        # else:
        return self.cell

    def estimateCellError(self):
        """
        The cell error is calculated.
        The sides of the cell get an error of 0.1% and
        the cell angles get a deviation of 0.05 if the angles are not within an range of 0.01 degrees to 90 degrees.
        :return:
        """
        for i in self.cell[:3]:
            self.cellError.append(i * 0.001)
        for i in self.cell[3:]:
            if 89.99 < i < 90.01:
                self.cellError.append(0.00)
            else:
                self.cellError.append(0.05)

    def getCellError(self):
        return self.cellError

    def extractSpaceGroup(self):
        """
        The full International Tables Hermann-Mauguin symbol is used, e.g., P 1 21 1 instead of P 21.
        The Hermann-Mauguin space group symbol is given without parenthesis, e.g., P 43 21 2.
        Please note that the screw axis is described as a two digit number.
        Even when no PDB entry today is fulfilling this criteria,
        the main whether the space group is starting with an 'R" even if it is rhombohedral obverse on hexagonal axes
        remains in the program.
        """
        doNotReplace1 = ['P 1', 'A 1', 'B 1', 'C 1', 'I 1', 'F 1']
        if not options['s']:
            try:
                self.spaceGroup = self.crystLine[55:66].strip('\n')
                if self.spaceGroup.strip() not in doNotReplace1:
                    self.spaceGroup = self.spaceGroup.replace(' 1 ', '').replace(' 1', '').lstrip()
                if self.spaceGroup[1] == "R":
                    x = self.cell[6] - self.cell[5]
                    if x >= 20:
                        self.spaceGroup[1] = "H"
                self.abbreviateSpaceGroup()
                if not self.validateSpaceGroup(self.getAbbrSpaceGroup()):
                    self.askSpaceGroup()
            except IndexError:
                self.askSpaceGroup()

    def getSpaceGroup(self):
        if self.spaceGroup:
            return self.spaceGroup
        else:
            return ''

    def askSpaceGroup(self):
        """
        Asks the user for input of a space group when called. The entered return is stripped of white space and all
        single one's are replaced when the return is not one of the space groups in doNotReplace1. Also the space group
        is checked when it starts with 'R': The cell parameters should fit, otherwise the first letter is replaced by H.
        Last, the space group is abbreviated by calling the function abbreviateSpaceGroup to call on the function
        validateSpaceGroup using the abbreviated function as param. Only a correct space group, validated, can break
        the while loop.
        :return:
        """
        doNotReplace1 = ['P 1', 'A 1', 'B 1', 'C 1', 'I 1', 'F 1']
        while True:
            answer = raw_input("\nSpace group not found or incorrect in .pdb file. Space group given: {} \n"
                                        "Please enter the correct space group:".format(self.getSpaceGroup()))
            # print self.spaceGroup.strip(), type(self.spaceGroup)

            if answer.strip() == '':
                self.spaceGroup = self.getSpaceGroup()
            else:
                if answer.strip() not in doNotReplace1:
                    self.spaceGroup = answer.replace(' 1 ', '').lstrip()
                if answer[1] == "R":
                    x = self.cell[6] - self.cell[5]
                    if x >= 20:
                        answer[1] = "H"
                    self.spaceGroup = answer
                self.abbreviateSpaceGroup()
                if self.validateSpaceGroup(self.getAbbrSpaceGroup()):
                    break

    def abbreviateSpaceGroup(self):
        """
        The space group is stripped of all white space. Also, when program crashes at this point, the pdb file does not
        contain the correct keywords or is empty.
        :return:
        """
        try:
            self.shortenedSpaceGroup = self.spaceGroup.replace(" ", "")
        except AttributeError:
            print 'ERROR: File does not contain X-ray diffraction data or not in the correct format. '
            exit()
        # print 'This is the space group: ', self.shortenedSpaceGroup

    def getAbbrSpaceGroup(self):
        short = self.shortenedSpaceGroup.rstrip('\r\n')
        return short

    def validateSpaceGroup(self, abbrSpaceGroup):
        """
        takes the abbreviated space group and returns either True or False, depending on the return of testSpaceGroup
        function.
        :param abbrSpaceGroup: string, abbreviated space group
        :return: Boolean
        """
        if testSpaceGroup(abbrSpaceGroup):
            return True
        else:
            return False

    def estimateLattice(self):
        """
        Takes the space group extracted from the pdb file and extracts the lattice:
        The first letter in space group notation gives the lattice system, which is extracted and translated into a
        number.
        The number will be printed after the :LATT: command in the .ins file.
        :return: Lattice as a number.
        """
        self.getSpaceGroup()
        latticeDict = {'P': '-1', 'I': '-2', 'H': '-3', 'F': '-4', 'A': '-5', 'B': '-6', 'C': '-7'}
        try:
            self.lattice = latticeDict[self.spaceGroup[0]]
        except KeyError:
            self.lattice = '-1'

    def getLattice(self):
        if not self.lattice:
            self.estimateLattice()
        return self.lattice

    def extractZvalue(self):
        """
        The Z value is the number of polymeric chains in a unit cell. In the case of heteropolymers,
        Z is the number of occurrences of the most populous chain.
        """
        validZ = ['0.5', '1', '2', '3', '4', '5', '6', '8', '9', '10', '11', '12', '16', '24', '36']
        # print 'this is zvalue',  self.crystLine[66:70].strip()
        # if options['z']:
        #     z = options['z']
        #     try:
        #         self.zValue = int(z)
        #     except ValueError:
        #         self.zValue = None
        #     except TypeError:
        #         self.zValue = None
        if not options['z']:
            try:
                self.zValue = int(self.crystLine[66:70])
            except KeyError:
                print 'No Z value found in PDB file.'
                self.zValue = self.askZvalue()
                # self.zValue = raw_input("\nPlease enter Z (number of molecules per cell):")
                # if not self.zValue:
                #     print ' *** Error: No valid Z value given. *** '
                #     self.zValue = None
            except ValueError:
                print 'No Z value was found in PDB file.'
                self.zValue = self.askZvalue()
            #     self.zValue = raw_input("Please enter Z (number of molecules per cell):")
            #     if not self.zValue:
            #         print ' *** Error: No valid Z value given. *** '
            #         self.zValue = None
            # if not options['i'] and not self.zValue:
            #     newZ = raw_input("Please enter Z (number of molecules per cell) [{}]:".format(self.zValue))
            #     if not newZ:
            #         self.zValue = 1
            #     # if newZ in validZ:
            #     else:
            #         self.zValue = int(newZ)
            # else:
            #     pass
        else:
            z = options['z']
            try:
                self.zValue = int(z)
                print 'INFO: Z value is set to {}'.format(z)
            except ValueError:
                self.zValue = 1
                print 'INFO: Z value was set to 1. Please check!'
            except TypeError:
                self.zValue = 1
                print 'INFO: Z value was set to 1. Please check!'

    def getZvalue(self):
        return int(self.zValue)

    def askZvalue(self):
        """
        this function is called if no z value could be extracted from the pdb file. the function is called as long as no
        z value has been specified.
        :return:
        """
        if options['i']:
            self.zValue = 1
            print 'INFO: Z value was set to 1. Please check!'
            return self.zValue
        else:
            while not self.zValue:
                newZ = raw_input("\nPlease enter Z (number of molecules per cell):")
                if not newZ:
                    pass
                else:
                    try:
                        self.zValue = int(newZ)
                        return self.zValue
                    except ValueError:
                        pass
                    except TypeError:
                        pass

    def extractWavelength(self):
        """
        searches the remark200 lines of the pdb file for the wavelength. If the wavelength is found, 
        the user is asked if the wavelength is correct and has a chance to correct it.
        If a wavelength could not ge found, the user is asked to enter a wavelength.
        :return: wavelength in Angstrom.
        """
        if options['w']:
            try:
                self.wavelength = float(options['w'])
            except TypeError:
                self.wavelength = None
            except ValueError:
                self.wavelength = None
        for line in self.remark200Lines:
            if not self.wavelength and "WAVELENGTH" in line:
                try:
                    wavelengthExtracted = float(line.split(':')[-1].lstrip(' ').split(' ')[0].partition(';')[0])
                    print '\nINFO: Wavelength found in pdb file: {}'.format(wavelengthExtracted)
                except ValueError:
                    print '\nAttention: PDB2INS could not find a valid wavelength in the pdb file.\n' \
                          'Information found:\n', line
                    continue
                if not options['i']:
                    reply = raw_input("Is the wavelength {} correct? (y or n) [Y]: ".format(wavelengthExtracted)).lower()
                    if reply == 'y' or reply == 'Y':
                        self.wavelength = wavelengthExtracted
                        return
                    if reply == 'N' or reply == 'n':
                        self.wavelength = None
                    if not reply:
                        self.wavelength = wavelengthExtracted
                else:
                    self.wavelength = wavelengthExtracted
                # else:
                #     self.wavelength = options['w']
        while not self.wavelength:
            # # this part is just for the automated pdb main, please remove after main!
            # if options['i']:
            #     self.wavelength = 1.54178
            # else:
            #     # this is the end of the change for the automated pdb main. Please remove after main!
            try:
                reply = raw_input("\nNo wavelength found in file. Enter wavelength in Angstroms [1.54178]: ")
                if not reply:
                    self.wavelength = 1.54178
                    # print "Wavelength set to default value 1.54178 Angstroms."
                elif 0.2 < float(reply) < 6.0:
                    self.wavelength = float(reply)
                else:
                    self.wavelength = None
            except ValueError:
                self.wavelength = 1.54178
                print " *** ERROR: Wavelength set to default value 1.54178 Angstroms. *** "
        self.checkWavelength()

    def checkWavelength(self):
        """
        checks the wavelength given if shelxl knows it. Shelxl recognizes kalpha radiation from copper, molybdenum,
        silver, gallium and Indium (list in this order).
        Otherwise DISP instructions must be made and the function makeKissel is called to read in the data necessary.
        :return:
        """
        # print 'starting check wavelength'
        if not any([abs(self.wavelength - shelxWL) <= 0.001 for shelxWL in [1.541867, 0.710747, 0.5608, 1.34139,
                                                                            0.513590]]):
            self.makeDISP = True
            # print 'the wavelength is not a shelxl recognized wavelength', self.makeDISP
            self.makeKissel()

    def getWavelength(self):
        return self.wavelength

    def extractScale(self):
        """
        reads the scale line from the pdb file and extracts the scale instructions in these lines to numpy arrays.
        :return: 2 np.arrays (scale, u)
        """
        for line in self.scaleLine:
            if line[:6] == "SCALE1":
                s11 = float(line[11:21])
                s12 = float(line[21:31])
                s13 = float(line[31:41])
                u1 = float(line[46:56])
                # print s11, s12, s13, u1
                x = np.array([s11, s12, s13])
            if line[:6] == "SCALE2":
                s21 = float(line[11:21])
                s22 = float(line[21:31])
                s23 = float(line[31:41])
                u2 = float(line[46:56])
                # print s21, s22, s23, u2
                y = np.array([s21, s22, s23])
            if line[:6] == "SCALE3":
                s31 = float(line[11:21])
                s32 = float(line[21:31])
                s33 = float(line[31:41])
                u3 = float(line[46:56])
                # print s31, s32, s33, u3
                z = np.array([s31, s32, s33])
        u = np.array([u1, u2, u3])
        scale = np.concatenate((x, y, z)).reshape(3, 3)
        # print scale, u
        return scale, u

    def handleResidueSequence(self, sequenceDict):
        """

        :return:
        """
        if self.sequenceLines:
            self.extractResiSequence()
        self.sequenceDict = sequenceDict

    def extractResiSequence(self):
        """
        This function extracts the aa sequence as given in the SEQREF record of the pdb file. The aa sequence is
        stored by chain as list in  dictionary 'resiDict'.
        :return: dict, aa sequence list sorted by chain
        """
        chainIDold = None
        aastring = []
        for line in self.sequenceLines:
            chainIDnew = line[11]
            if not chainIDnew == chainIDold and chainIDold:
                self.resiDict[chainIDold] = aastring
                chainIDold = chainIDnew
                aastring = []
            for i in xrange(13):
                i *= 4
                i += 19
                chunk = line[i:i+3].strip()
                if chunk:
                    aastring.append(chunk)
            chainIDold = chainIDnew
        self.resiDict[chainIDold] = aastring
        # for key in self.resiDict.keys():
        #     print 'done chain', key, len(self.resiDict[key]), self.resiDict[key]

    def getResiDict(self):
        if self.resiDict:
            return self.resiDict
        else:
            return self.sequenceDict

    def printResiSequence(self):
        for key in self.resiDict.keys():
            print 'chain', key, self.resiDict[key]

    def makeGeneralRefinementInstructions(self):
        """
        General refinement instructions are implemented here.
        In three parts (General refinement, general Refinement2, general Refinement3) the refinement instruction, which
        are present at the begin of the .ind file are listed.
        general Refinement is used for all fix instructions from DEFS to WPDB.
            CGLS is set to 0 to calculate fcf but nor refine for solvent analysis.
        generalRefinement2 is used to make DELU and SIMU instructions. Herefore the the elementList2 from AtomContainer
        is used to give a list of all element in the pdb, swhich are relevant for these instructions.
        generalRefinement3 is used to give the fixed instructions from BUMP to MORE.
        :return: 3 lists with general refinement instructions.
        """
        self.generalRefinement = ['DEFS 0.02 0.1 0.01 0.04', 'CGLS 0 -1', 'SHEL 999 0.1', 'FMAP 2', 'PLAN 200 2.3',
                                  'LIST 6', 'WPDB 2 \n']
        self.generalRefinement2 = {'C': '$C_*', 'N': '$N_*', 'O': '$O_*', 'S': '$S_*', 'P': '$P_*'}
        self.generalRefinement2Order = ['C', 'N', 'O', 'S', 'P']
        self.generalRefinement3 = ['XNPD 0.001 \n', 'REM BUMP', 'SWAT \n',
                                   'REM Remove MERG 4 instruction if Friedel opposites should not be merged.',
                                   'MERG 4 \n',
                                   'REM MORE 0 would reduce output if not required for diagnostic purposes. \n']

    def getRIGUInstructions(self, elementlist):
        """
        The RIGU instruction should replace DELu and SIMU instruction.
        [Thorn, Dittrich & Sheldrick, Acta Cryst. A68 (2012) 448-451]
        :param elementList: AtomContainer generated List of elements of interest for this instruction (elementList2)
        :return: string
        """
        return 'RIGU  !Apply enhanced rigid body restraints\n'
        # return 'RIGU {}\n'.format(' '.join([self.generalRefinement2[key] for key in self.generalRefinement2Order
        #                                     if key in elementlist]))

    def getDELUInstruction(self, elementlist):
        """
        The DELU instructions is generated from the ElementList2 and the generalRefinement2 in combination with
        generealRefinement2Order.
        :param elementlist: AtomContainer generated List of elements of interest for this instruction (elementList2)
        :return: string, DELU instruction
        """
        return 'DELU {}\n'.format(' '.join([self.generalRefinement2[key] for key in self.generalRefinement2Order
                                            if key in elementlist]))

    def getSIMUInstruction(self, elementlist):
        """
        The SIMU instructions is generated from the ElementList2 and the generalRefinement2 in combination with
        generealRefinement2Order.
        :param elementlist: AtomContainer generated List of elements of interest for this instruction (elementList2)
        :return: string, DELU instruction
        """
        return 'SIMU 0.1 {}\n'.format(' '.join([self.generalRefinement2[key] for key in self.generalRefinement2Order
                                                if key in elementlist]))

    def getGeneralRefinementInstructions(self):
        """
        The generalRefinement list is turned into a string for the .ins file
        :return: string
        """
        return '\n'.join(self.generalRefinement)

    def getGeneralRefinementInstructions3(self):
        return '\n'.join(self.generalRefinement3)

    def makeKissel(self):
        """"
        Imports the file kissel.py, which contains a dictionary. Keys = Elements and value = list of lists with f', f''
        and mu for different wavelength.
        From those tables the correct parameters for the given wavelength can be calculated.
        """
        # kisselwave = open('kisselwave.txt', 'r')
        # self.kisselContent = {element[:10].split()[1].rstrip('\n'): KisselData(element.split('\n')[1:])
        #                       for element in kisselwave.read().split('#S')[1:]}
        import kissel
        self.kisselContent = kissel.kisselData

    def makeDISPinstruction(self, elementList):
        """
        take the elementList and the wavelength to search the dictionary in kissel.py dictionary kisselData.
        KisselData has as key the element and as value a list of lists, where each list is a line from kisselwave.txt.
        Each line has wavelength, f', f'' and mu.
        The function searches for the two lines in the file for each element present in elementList in which the
        wavelength lies in between. Than a linear interpolation is used to find the correct f', f'' and mu for the
        wavelength.
        :param elementList: list of all elements present in pdb
        :return:
        """
        # print 'make Disp instructions has been called'
        wl = float(self.wavelength)
        for i in elementList:
            try:
                previous, next = self.kisselContent[str(i.capitalize())][wl]
                delta = previous[0] - next[0]
                f = (next[0] - wl) / delta
                fprime = abs(f) * (previous[1] - next[1]) + next[1]
                    # (1 - f) * next[1] + f * previous[1]
                # f * next[1] + f * previous [1]
                fdoubleprime = abs(f) * (previous[2] - next[2]) + next[2]
                mu = abs(f) * (previous[3] - next[3]) + next[3]
                self.dispInstructions.append('DISP ${:2} {:10.5f} {:10.5f} {:12.5f}\n'.format(i, fprime,
                                                                                              fdoubleprime, mu))
            except KeyError:
                pass

        # previous, next = self.kisselContent['O']['1.2345']

    def getDISPinstructions(self):
        return self.dispInstructions

# class KisselData(object):
#     def __init__(self, dataList):
#         self.dataList = [[float(i) for i in line.split()] for line in dataList]
#
#     def __getitem__(self, item):
#         previous = None
#         for i, line in enumerate(self.dataList[1:]):
#             previous = line
#             if line[0] < item:
#                 next = self.dataList[i+2]
#                 return previous, next


class AtomContainer(object):
    waterNames = ["HOH", "WAT", "OH2", "H2O"]
    naturalAA = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met',
                 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']

    def __init__(self):
        # print 'Initilized AtomContainer'
        self.atomDict = {}
        self.ResiList = []
        self.atomAnisouDict = {}
        self.elementDict = {}
        self.elementList = []
        self.elementList2 = []
        self.chainIDSet = set()
        self.residueSet = set()
        self.otherResiSet = set()
        self.ssBondList = []
        self.otherResiList = []
        self.chains = {}
        self.nTerminalNAtoms = {}
        self.waterResiNumbers = []
        self.cTerminalCAtoms = {}
        self.oldResiNum = None
        self.terminusRestraints = []
        self.cysSAtomsList = []
        self.nonNaturalAAAtomList = []
        self.waterAtomList = []
        self.specifiedAtomList = []
        self.removeAtomsList = []
        self.omitList = []
        self.resicounter = None
        self.hetIDlist = []
        self.ssBonds = False
        # self.overlongChains = {}
        self.oldResiName = None
        self.ligandChain = False
        self.incompleteResiString = []
        # self.overlongChainList = []
        # self.printWarning = True
        # self.overlongResiSeqNumbers = {}
        self.chainIDbefore = None
        # self.resicounter2 = None
        self.resiNumberBefore = None
        self.resiNameBefore = None
        self.chainIDdict = {}
        self.sequenceDict = {}
        self.chainCounter = 0
        self.firstNterm = True
        self.neut = False
        self.keepHAtoms = False
        self.resiNumberDict = {}
        self.resiOffset = 0
        self.resiSeqOffsetDict = {}
        self.insertionCodeDict = {}
        self.insCodeOccCounter = 0
        self.resiTuple = ()
        # self.insertionCodeBefore = None
        self.insertionCodeBefore = False
        self.insertionCodeString = []
        self.resiNameDict = {}
        # self.waterCounter = 0
        self.negResiNumber = False
        self.overlongResiNum = False
        self.wrongResiName = False
        self.altLocDict = {}
        self.wateroffset = None
        self.waterOffsetWarning = False
        self.resiNUmberCollision = False

    def extractAtom(self, line):
        """
        An Atom is created from the lines of the PDB file with 'Atom' mark
        sorts them with their atom name into a dictionary.
        * if a resinumber keeps but the resiname changes and an insertioncode is given, the resioffset is increased and
        with it also the residue number of this residue and all following it.
        * if the residue number contains letters all residues are renumbered that come afterwards.
        [This  following only applies to the old shelxl version]
        If the residue number is greater than 999 shelxl cannot work with the resulting file. (Shelxl does not recognize
        chain IDs therefore every chain gets an offset off a thousand.) Should it occur that a chain has more than 999
        residues, a new chain will be created with the remaining residues. Right now this happens automatically for all
        water residues later on and here for not bound residues listed in HET record.
        Important: An error should occur if the chain is neither! Now the chain just starts again with resi#2!!!
        [The former part only applies to the old shelxl version]
        * MAJOR CHANGES added with revision 173, now insertion Codes are handled differently. Still some issues with
        hole chains with insertion codes!
        :param line: one line with atom data from pdb file.
        :return: dictionary entry with atom name as key and the corresponding atoms line from the pdb file as value.
        """
        # print 'you are in extract atom'
        changedoffset = False
        newAtom = Atom(line)
        # print 'atom object created'
        self.chainIDSet.add(newAtom.getChainID())
        atomName = newAtom.getAtomName()
        self.atomDict[atomName] = newAtom
        chainIDnew = newAtom.getChainID()  # This line must be removed when chainIDs are renumbered with code below.

        # This part renumbers the chainIDs in sequence, relict from shelxl version with only large letter chainIDs
        # chainID = newAtom.getChainID()
        # chainIDnew = self.getChainID(chainID)
        # # print chainID, chainIDnew
        # newAtom.setChainID(chainIDnew)

        resiNumber = newAtom.getResiSeqNum().strip()
        resiSeqOffset = None
        # here residue names not starting with a letter are renamed
        resiName = newAtom.getResidueName().strip()
        # print chainIDnew, resiNumber, resiName

        # the following part handles residue names starting with a number
        # if resiName:
        #     if not resiName[0].isalpha():
        #         self.changeResiName(resiName)
        #         # resiNameNew = self.changeResiName(resiName)
        # else:
        #     print ' *** ERROR: The file has not the expected format. Data type: residue name is missing. ***'
        #     exit()

        if not resiName or not chainIDnew:
            print ' *** ERROR: The file has not the expected format. Data type: residue name is missing. ***'
            exit()

        # the following part handles residue names consisting only of numbers
        try:
            x = int(resiName)
            self.changeResiName(resiName)
            resiName = newAtom.getResidueName()
            # resiNameNew = self.changeResiName(resiName)
        except ValueError:
            pass

        # the following part handles negative residue seq numbers at the beginning of the chain
        if int(resiNumber) < 0:
            self.negResiNumber = True  # Warning is printed after parsing the pdb file.
        # try:
        #     i = self.resiSeqOffsetDict[chainID]
        # except KeyError:
        #     i = None
        #     self.resiSeqOffsetDict[chainID] = 0
        #     pass
        # if int(resiNumber) < 0:
        #     if i and abs(int(resiNumber)) > i:
        #         print "\nAttention: Residue numbers for chain ID {} might not be correct. " \
        #               "Please check!".format(chainID)
        #         resiSeqOffset = abs(int(resiNumber))
        #         self.resiSeqOffsetDict[chainID] = resiSeqOffset
        #     elif not i:
        #         resiSeqOffset = abs(int(resiNumber))
        #         self.resiSeqOffsetDict[chainID] = resiSeqOffset
        #     else:
        #         pass
        #         # if abs(int(resiSeqOffset)) > int(self.resiSeqOffsetDict[chainID]):
        # try:
        #     i = self.resiSeqOffsetDict[chainID]
        #     newAtom.setResiSeqNum(int(resiNumber) + i)
        # except KeyError:
        #     pass
        # resiNumber = newAtom.getResiSeqNum().strip()
        # resiName = newAtom.getResidueName()

        # In the next part insertion codes are handled, assuming complete residues are inserted.
        insertionCode = newAtom.getInsertionCode()
        # try:
        #     self.resiNumberDict[chainID].append(resiNumber)
        # except KeyError:
        #     self.resiNumberDict[chainID] = resiNumber
        if chainIDnew != self.chainIDbefore:
            self.insCodeOccCounter = 0
            self.resiOffset = 0
        resiNumber = int(resiNumber)
        resiNumberNew = self.resiOffset + resiNumber

        # handling insertion codes with an offset of a multiple of 1000:
        if insertionCode:
            resiNumberNew = self.handleInsertionCode(insertionCode, resiNumber, resiNumberNew)
        newAtom.setResiSeqNum(resiNumberNew)  # it is very important to set the new value for the variable in the object
        resiNumber = int(newAtom.getResiSeqNum())
        # print 'insertion code handled'

        # all residues with insertion code need restraints to bind them to their nearest residue
        if resiNumberNew != self.resiNumberBefore and not resiName == 'HOH':
            if (not insertionCode and self.insertionCodeBefore) or insertionCode:
                try:
                    delta = int(resiNumberNew)-int(self.resiNumberBefore)
                    # print delta, chainIDnew, resiNumberNew, self.chainIDbefore, self.resiNumberBefore
                    if abs(delta) >= 20 and self.chainIDbefore == chainIDnew:
                        self.insertionCodeRestraints(self.chainIDbefore, self.resiNumberBefore, chainIDnew,
                                                     resiNumberNew)
                except TypeError:
                    pass

        # print 'insertion Code restraints'

        # insertionCodeOffset = 0
        # j = 0
        # if insertionCode:
        #     # print resiNumber, chainIDnew, self.insertionCodeDict
        #     # if not resiNumber == self.resiNumberBefore:
        #         # print 'residue number changed, has insertion code'
        #     try:
        #         insertionCodeOffset = self.insertionCodeDict[chainIDnew]
        #     except KeyError:
        #         self.insertionCodeDict[chainIDnew] = 0
        #     if resiNumber == self.resiNumberBefore:
        #         # if insertionCode:
        #         print "Found Residue insertion code.", chainIDnew, resiNumber, resiName
        #         if not insertionCode == self.insertionCodeBefore or not self.resiNameBefore == resiName:
        #             # if not self.resiNumberBefore:
        #             # print 'the offset was changed here, 1'
        #             self.insertionCodeDict[chainIDnew] += 1
        #             self.resiOffset += 1
        #             changedoffset = True
        #             self.insertionCodeBefore = insertionCode
        #         # elif not resiName == self.resiNameBefore:
        #         #     # print 'the offset was changed here, 2'
        #         #     self.insertionCodeDict[chainIDnew] += 1
        #         #     self.resiOffset += 1
        #         #     changedoffset = True
        #         #     self.insertionCodeBefore = insertionCode
        #     if not resiNumber == self.resiNumberBefore:
        #         print 'resinum changes while insertioncode', resiNumber, self.resiNumberBefore, self.insCodeOccCounter
        #         try:
        #             i = int(self.resiNumberBefore) + 1
        #         except TypeError:
        #             i = 1
        #         k = i + 1
        #         if i > int(resiNumber):  # when the resinumber before is greater than the new resi number
        #             # print 'resinumber before greater than resinumber'
        #             j = self.insCodeOccCounter
        #             if j < i:
        #                 # print 'i > j ', i, j
        #                 difference = abs(i - int(self.resiNumberBefore))
        #                 self.insCodeOccCounter += difference
        #                 # print 'iiiiii', self.insCodeOccCounter, i, j, difference
        #             if j >= i:
        #                 self.insCodeOccCounter += 1
        #             #self.insCodeOccCounter += difference
        #         elif k < int(resiNumber):  # when there is a jump in resinumbers in combination with insertion code
        #             difference2 = abs(int(self.resiNumberBefore) - int(resiNumber))
        #             if self.insCodeOccCounter >= (int(resiNumber)):
        #                 pass
        #             else:
        #                 self.insCodeOccCounter += difference2
        #         # if i == j:
        #         #     self.insCodeOccCounter += 1
        #         else:  # i have no idea when this part should happen
        #             self.insCodeOccCounter += 1
        #         # print 'yyy', resiNumber, self.resiNumberBefore, self.insCodeOccCounter
        # if chainIDnew != self.chainIDold:  # new chain id means the insCodeOccCounter should be reset
        #     self.insCodeOccCounter = 0
        # # print self.insCodeOccCounter, resiNumber, self.resiNumberBefore
        # if not insertionCode and self.insertionCodeBefore:  # change from insertion code to no ins code, same chain
        #     if not resiNumber == self.resiNumberBefore:
        #         self.insCodeOccCounter += 1
        #     try:
        #         # self.insCodeOccCounter += 1
        #         self.insertionCodeDict[chainIDnew] += self.insCodeOccCounter
        #         self.insCodeOccCounter = 0
        #     except KeyError:
        #         self.insertionCodeDict[chainIDnew] = 0
        # if not insertionCode == self.insertionCodeBefore:  # change between different insertions codes in one chain
        #     try:
        #         self.insertionCodeDict[chainIDnew] += self.insCodeOccCounter
        #         self.insCodeOccCounter = 0
        #     except KeyError:
        #         self.insertionCodeDict[chainIDnew] = self.insCodeOccCounter
        #         self.insCodeOccCounter = 0
        #     self.resiOffset += 1
        #     self.insertionCodeDict[chainIDnew] += 1
        #     changedoffset = True
        #     self.insertionCodeBefore = insertionCode
        # elif not self.resiNameBefore or not self.resiNumberBefore:
        #     pass
        # else:
        #     pass
        # if not changedoffset:
        #     if resiNumber == self.resiNumberBefore or not self.resiNumberBefore:
        #         pass
        #     else:
        #         if not resiNumber.isdigit():
        #             self.insertionCodeDict[chainID] += 1
        #             self.resiOffset += 1
        #         else:
        #             pass
        # print self.insertionCodeDict
        # print 'zzzz', resiNumber

        # In the following part water residues are separated to a new chain.
        # if resiName == 'HOH':
        #     self.waterCounter += 1
        #     resiNumberNew = self.waterCounter
        #     newAtom.setResiSeqNum(resiNumberNew)

        # This part is older!
        # print 'I changed the resiNumber of a water residue. ', resiNumberNew
        # try:
        #     #newAtom.setResiSeqNum(int(''.join([c for c in resiNumber if c.isdigit()])) +
        #                                       self.insertionCodeDict[chainIDnew])
        #     newAtom.setResiSeqNum(int(resiNumber) + self.insertionCodeDict[chainIDnew])
        # except KeyError:
        #     newAtom.setResiSeqNum(int(''.join([c for c in resiNumber if c.isdigit()])) + self.resiOffset)
        #     print (int(''.join([c for c in resiNumber if c.isdigit()])) + self.resiOffset)
        #     print newAtom.getResiSeqNum()

        # if not self.chainIDold:
        #     self.resicounter2 = 1
        # if chainIDnew == self.chainIDold:
        #     if resiNumber == self.resiNumberBefore:
        #         pass
        #     else:
        #         self.resicounter2 += 1
        # print chainID, resiNumberNew, newAtom.getAltLoc(), insertionCode
        # print self.insCodeOccCounter, resiNumberNew

        # at this part water residues are renumbered if their numbers would collide with the insertion code offset
        if self.insCodeOccCounter and int(resiNumberNew) >= 1000:
            if resiName == 'HOH' and not self.wateroffset:
                self.waterOffsetWarning = True
                self.wateroffset = self.insCodeOccCounter + 1
                newAtom.setResiSeqNum(resiNumberNew + self.wateroffset * 1000)
                resiNumberNew = newAtom.getResiSeqNum()
            elif resiName == 'HOH' and self.wateroffset:
                self.waterOffsetWarning = True
                newAtom.setResiSeqNum(resiNumberNew + self.wateroffset * 1000)
                resiNumberNew = newAtom.getResiSeqNum()
                pass
            elif resiName != 'HOH':
                self.resiNUmberCollision = True
            else:
                print 'Error handling insertion codes encountering chains with more than 1000 residues.'
                exit()

        # the new atom is added to the chain and residue classes.
        try:
            # print 'i am in try', chainIDnew, resiNumberNew
            success = self.chains[chainIDnew][resiNumberNew].append(newAtom)
            # print 'this was created: ', chainIDnew, resiNumberNew, success
            if not success:
                # print 'not success'
                if newAtom.getAltLoc():
                    success = self.chains[chainIDnew][resiNumberNew].append(newAtom, True)
                    # print 'new success', success
                    if not success:
                        print '\nERROR: Problem while handling insertion codes starting with residue: ', chainIDnew, \
                            ':', resiNumber, resiName
                        exit()
                else:
                    # print 'else'
                    self.resiOffset += 1
                    resiNumberNew = self.resiOffset + resiNumber
                    # print resiNumber, resiNumberNew
                    newAtom.setResiSeqNum(resiNumberNew)  # never forget to set the new residue number for the object!!!
                    # print 'bla', chainIDnew, resiNumberNew, newAtom, self.chains
                    success = self.chains[chainIDnew][resiNumberNew].append(newAtom)
                    # print 'No, this was created: ', chainIDnew, resiNumberNew
                    if not success:
                        # print line
                        # for key, value in self.chains.items():
                        #     print key, value
                        #     for i, j in value.items():
                        #         print i, j.getResiName(), j.getResiSeqNum()
                        print '\nERROR: Problem while handling insertion codes starting with residue: ', chainIDnew, \
                            ':', resiNumber, resiName
                        exit()
        except NoResidueError:  # creates the residue if it was not there before, happens with first atom of new residue
            self.chains[chainIDnew][resiNumberNew] = Residue([newAtom])
            # print 'No, no: This was created: ', chainIDnew, resiNumberNew
        except KeyError:  # creates a new chain and new residue object
            newChain = Chain()
            self.chains[chainIDnew] = newChain
            newChain[resiNumberNew] = Residue([newAtom])
            # print 'You are all wrong, this was created: ', chainIDnew, resiNumberNew
        # self.insertionCodeBefore = insertionCode
        if insertionCode:
            self.insertionCodeBefore = True
        else:
            self.insertionCodeBefore = False
        self.resiNumberBefore = resiNumberNew
        self.resiNameBefore = resiName
        self.chainIDbefore = chainIDnew
        # print 'finished create atom'

    def changeResiName(self, resiName):
        """
        this function is called when a residue has a residue name starting with a number. SHELXL cannot handle residues
        starting with a number, the name must start with a letter. It is possible to rename the residue with the same
        name, without the program complaining. There are restraints available for ligands starting with a letter!
        :return:
        """
        try:
            resiNameNew = self.resiNameDict[resiName]
        except KeyError:
            # print resiName
            if not options['i']:
                while True:
                    resiNameNew = raw_input('\n *** WARNING: The residue {} has a name SHELXL cannot handle! ***\n'
                                            '     You can rename the residue now or keep the original name.\n'
                                            '     Please enter a new 3 digit residue name containing '
                                            'a letter [{}]: '.format(resiName, resiName))
                    if not resiNameNew:
                        resiNameNew = resiName
                    if len(resiName) <= 3:
                        print 'INFO: Residue {} successfully renamed to {}.'.format(resiName, resiNameNew)
                        break
            else:
                resiNameNew = resiName
                self.wrongResiName = True
                print '\nWARNING: The residue {} has a name SHELXL cannot handle! Please rename.\n'.format(resiName)
            self.resiNameDict[resiName] = resiNameNew
        # return resiNameNew

    def handleInsertionCode(self, insertionCode, resiNumber, resiNumberNew):
        """
        All atoms with an insertion code are handled here. The residue with the insertion code gets an offset. If the
        residue number becomes larger than 10000, a warning is printed. The insCodeOccCounter is needed for the water
        residue offset, should water residues collide with the new insertion code offset.
        :param insertionCode:
        :param resiNumber:
        :param resiNumberNew:
        :return:
        """
        try:
            a = self.insertionCodeDict[insertionCode]
        except KeyError:
            self.insCodeOccCounter += 1
            #  print len(self.insertionCodeDict)
            if len(self.insertionCodeDict) >= 9:  # offset should not be larger than 9000
                a = 9000 + (len(self.insertionCodeDict) - 8) * 100
                if a > 10000:
                    self.overlongResiNum = True  # Warning is printed after parsing the pdb file.
            else:
                a = (len(self.insertionCodeDict) + 1) * 1000
            self.insertionCodeDict[insertionCode] = a
            print '\nINFO: Residues with the insertion Code {} have now an offset of {} added to the residue ' \
                  'number!'.format(insertionCode, a)
        if resiNumberNew:
            resiNumberNew += a
        else:
            resiNumberNew = a + resiNumber
        return resiNumberNew

    def insertionCodeRestraints(self, chainIDbefore, resiNumberbefore, chainIDnew, resiNumberNew):
        """
        All inserted residues have to get the appropriate DFIX, DANG and FLAT restraints, since the general one's only
        work for residues with consecutive numbers.
        :param chainIDbefore:
        :param resiNumberbefore:
        :param chainIDnew:
        :param resiNumberNew:
        :return:
        """
        # print chainIDbefore, chainIDnew, resiNumberbefore, resiNumberNew
        self.insertionCodeString.append('DFIX 1.329 C_{}:{} N_{}:{}'.format(chainIDbefore, resiNumberbefore, chainIDnew,
                                                                            resiNumberNew))
        self.insertionCodeString.append('DANG 2.425 CA_{}:{} N_{}:{}'.format(chainIDbefore, resiNumberbefore,
                                                                             chainIDnew, resiNumberNew))
        self.insertionCodeString.append('DANG 2.250 O_{}:{} N_{}:{}'.format(chainIDbefore, resiNumberbefore, chainIDnew,
                                                                            resiNumberNew))
        self.insertionCodeString.append('DANG 2.435 C_{}:{} CA_{}:{}'.format(chainIDbefore, resiNumberbefore,
                                                                             chainIDnew, resiNumberNew))
        self.insertionCodeString.append('FLAT 2.0 O_{a}:{b} N_{c}:{d} '
                                        'C_{a}:{b} CA_{c}:{d}'.format(a=chainIDbefore, b=resiNumberbefore, c=chainIDnew,
                                                                      d=resiNumberNew))

    def getChainID(self, chainIDbefore):
        """
        The chainID  is checked if it is present in the chainID dictionary. if  not it is added to the chainIDdict.
        should the chainID be not a character, the chain ID is changed. Is this function necessary or plausible?
        :param chainIDbefore:
        :return:
        """
        try:
            chainIDnew = self.chainIDdict[chainIDbefore]
        except KeyError:
            self.chainCounter += 1
            if self.chainCounter >= 27:
                chainIDnew = self.chainCounter + 6
            else:
                chainIDnew = self.chainCounter
            self.chainIDdict[chainIDbefore] = chainIDnew
        return chr(chainIDnew + 64)

    # This function is not needed anymore, SHELXL now recognizes chainIDS.
    # def findOverlongChains(self):
    #     """
    #     this function looks at all chains and adds all chains longer than 999 residues to the
    #     dictionary overlongChains.
    #     :return:
    #     """
    #     for chainID, chain in self.chains.items():
    #         if len(chain) > 999:
    #             self.overlongChains[chainID] = chain
    #     return self.overlongChains

    def createSequenceDict(self):
        """
        Should the pdf file not contain a seq ref, the sequence should be extracted and saved to a dictionary by Chain.
        dictionary SequenceDict, key = ChainID, value = List of Residues, 3 letter code
        :return:
        """
        for chainID, chain in self.chains.items():
            chain.finalize()
            residueList = chain.getResiList()
            self.sequenceDict[chainID] = residueList
        # for key in self.sequenceDict.keys():
        #     print 'chain', key, 'residues', self.sequenceDict[key]

    def getSequenceDict(self):
        return self.sequenceDict

    def getTerminalResidues(self):
        """
        extracts all residues which are both the first(last) and natural aa residue in a chain.
        """
        for chain in self.chains.values():
            chain.finalize()
        nTerminals = [chain.getNTermResi() for chain in self.chains.values()]
        cTerminals = [chain.getCTermResi() for chain in self.chains.values()]
        return nTerminals, cTerminals

    def isBoundToSomethingDiff(self):
        """
        takes list of terminal natural aa residues from getTerminalResidues and gives each possible n terminal resi to
        the _nIsNotBound function and same for c terminal residue to _cIsNotBound function. If those functions return
        the residue, it is saved to a list. all the residues in the list are truely terminal residues and therefore need
        other restraints.
        :return:
        """
        n, c = self.getTerminalResidues()
        notBoundNResis = []
        notBoundCResis = []
        for resi in n:
            if not resi:
                continue
            unBoundResi = self._nIsNotBound(resi)
            if unBoundResi:
                notBoundNResis.append(unBoundResi)
        for resi in c:
            if not resi:
                continue
            unBoundResi = self._cIsNotBound(resi)
            if unBoundResi:
                notBoundCResis.append(unBoundResi)
        for i in notBoundCResis:
            for atom in i:
                if atom.getPDBAtomName() == 'O':
                    atom.rename('OT1')
        return notBoundNResis, notBoundCResis

    def _nIsNotBound(self, resi):
        """
        If the nitrogen has no other atom in direct distance, it is assumed that is the terminal nitrogen atom.
        :param resi:
        :return:
        """
        for atom in resi:
            if atom.getPDBAtomName() == 'N':
                atom1 = atom
                for atom2 in self.getOtherResiList():
                    if self.getAtomDistance(atom1.getAtomCoord(), atom2.getAtomCoord()) < 2.0:
                        return False
        return resi

    def _cIsNotBound(self, resi):
        """
        the IsBoundToSomethingDifferent function calls this function and gives a residue as an argument. Each atom of
        the residue is called and if the pdbatomname is 'C' the distance of this atom to all atoms of not natural aa
        residues is checked. if there is a close distance, the c atom is assumed coordinated and cannot be a terminal c
        atom. the function returns false. If no other atom close to C is found, the residue is assumed unbound and
        returned.
        (If the residue has an atom named OT1 it is assumed to be a terminal residue.)
        :param resi: Residue
        :return: Residue or None
        """
        for atom in resi:
            if atom.getPDBAtomName() == 'C':
                atom1 = atom
                for atom2 in self.getOtherResiList():
                    if self.getAtomDistance(atom1.getAtomCoord(), atom2.getAtomCoord()) < 3.0:
                        return False
        return resi

    def makeRestraintsForTermini(self, resiDict):
        """
        In this function the unbound terminal residues are called and split into N or C terminal residues. For each
        option the appropriate 'makeRestraints' function is called.
        :return:
        """
        self.terminusRestraints = []
        nTermini, cTermini = self.isBoundToSomethingDiff()
        for i in nTermini:
            self.makeRestraintsForNTerm(i)
        self.terminusRestraints.append('\n')
        for j in cTermini:
            # print 'make restraints for c termini for :', j
            self.makeRestraintsForCTerm(j, resiDict)

    def getRestraintsForTermini(self):
        return self.terminusRestraints

    def makeRestraintsForNTerm(self, i):
        """
        Here the N terminal nitrogen should get a 'REM HFIX 33' attribute. Is the residue the first residue in the
        chain, it is assumed to be the right n terminus. If not, the HFIX is added with preceeding command and
        :param i: residue
        :return: thats the question!
        """
        for atom in i:
            if atom.getPDBAtomName() == 'N':
                if self.firstNterm:
                    self.firstNterm = False
                    self.terminusRestraints.append('REM Remove "REM" before the following to generate HFIX for N '
                                                   'terminus.')
                chainID = atom.getChainID()
                num = atom.getResiSeqNum().strip()
                if num == '1':
                    self.terminusRestraints.append('REM HFIX 33 N_{}:{:{}f}'.format(chainID, float(num), padding))
                else:
                    self.terminusRestraints.append('\nREM Chain: {} residue: {} could be the N terminal residue.'.
                                                   format(chainID, num))
                    self.terminusRestraints.append('REM HFIX 33 N_{}:{:{}f}\n'.format(chainID, float(num), padding))

    def makeRestraintsForCTerm(self, j, resiDict):
        """
        Here it gets a little complicated. The c terminal residue should no longer have an oxygen named 'O'. This oxygen
        should carry the name 'OT1'. Also an oxygen named OT2 should be present. If the c terminus that was found is not
        the true c terminus (for example when the protein is not modeled any further but there should be more residues),
        nothing should be changed. But if it is the true c terminus, restraints like 'DFIX() 1.249 C OT1 C OT2' and
        others should be added. Is it possible to find out that the program has found the true c terminus and add the
        restraints without asking the user?
        :param j:
        :return:
        """
        chainID = j[0].getChainID()
        if chainID == 'XXX' or chainID == 'ZZZ':
            print 'ZZZ'
            return
        if j[0].getResidueName() == 'HOH':
            print 'HOH'
            return
        num = j[0].getResiSeqNum()
        nameList = [atom.getPDBAtomName() for atom in j]
        # for i in nameList:
        #     print i
        try:
            #  int(num) == len(resiDict[chainID]) or ()
            if 'OT1' in nameList and 'OT2' in nameList:
                self.terminusRestraints.append('DFIX 1.249 C_{0}:{1:{2}f} OT1_{0}:{1:{2}f} C_{0}:{1:{2}f} '
                                               'OT2_{0}:{1:{2}f}'.format(chainID, float(num), padding))
                self.terminusRestraints.append('DANG 2.379 CA_{0}:{1:{2}f} OT1_{0}:{1:{2}f} CA_{0}:{1:{2}f} '
                                               'OT2_{0}:{1:{2}f}'.format(chainID, float(num), padding))
                self.terminusRestraints.append('DANG 2.194 OT1_{0}:{1:{2}f} '
                                               'OT2_{0}:{1:{2}f}\n'.format(chainID, float(num), padding))
                # renaming is not possible anymore at this point, since this function is called after the shelxstring
                # of all atoms is already written. The change of the position for writing the shelxstring was changed
                # for some reason, but i can't remember it. Now all O atoms called 'OXT' are automatically renamed
                # to 'OT1' right when the atom object is created.
                # Addition:
                # Now the shelxstring is written again before it being added to the string that goes into the .ins file.
                # So all terminal oxygen atoms should be named appropriately and restrains should only be applied
                # without REM when both oxygen atoms (OT1, OT2) are present.
                # for atom in j:
                #     if atom.getPDBAtomName() == 'OXT':
                #         self.renameCtermOatom(atom)
            else:
                self.terminusRestraints.append('REM Chain: {} residue: {} could be the C terminal residue.'.
                                               format(chainID, num))
                self.terminusRestraints.append('REM Remove "REM" before the following to generate restraints '
                                               'for C terminus.')
                self.terminusRestraints.append('REM DFIX 1.249 C_{0}:{1:{2}f} OT1_{0}:{1:{2}f} C_{0}:{1:{2}f} '
                                               'OT2_{0}:{1:{2}f}'.format(chainID, float(num), padding))
                self.terminusRestraints.append('REM DANG 2.379 CA_{0}:{1:{2}f} OT1_{0}:{1:{2}f} CA_{0}:{1:{2}f} '
                                               'OT2_{0}:{1:{2}f}'.format(chainID, float(num), padding))
                self.terminusRestraints.append('REM DANG 2.194 OT1_{0}:{1:{2}f} '
                                               'OT2_{0}:{1:{2}f}\n'.format(chainID, float(num), padding))

        except KeyError:
            pass

    def renameCtermOatom(self, atom):
        atom.rename('OT1')

    def extractHetRecord(self, line):
        self.hetIDlist.append(line[7:10])

    def extractAtomAnisou(self, line):
        """
        An Atom is created from the lines of the PDB file with 'Anisou' mark
        sorts them with their atom name into a dictionary.
        :param line: one line with anisou data from pdb file.
        :return: dictionary entry with atom name as key and the corresponding atoms anisou line from the pdb file
         as value.
        """
        newAtomAnisou = Atom(line)
        try:
            self.atomDict[newAtomAnisou.getAtomName()].overrideADP(newAtomAnisou.getRawADP())
        except KeyError:
            print ' *** Error: Found ANISOU record for unknown atom. *** '
            exit(1)

    def extractAllElements(self):
        """
        The Element of the Atom is iterated over in element Dict and the occurrences are summed up.
        The Elements are sorted by relevance into the order C, H, N, O, S, others.
        Also the elements relevant for the DELU and SIMU instruction are found in elementList2.
        :return: molecular formula possible, dict of al elements in pdb file with number of occurrences.
        """
        for thisAtom in self.atomDict.values():
            newElement = thisAtom.getAtomElement()
            try:
                self.elementDict[newElement] += 1
            except KeyError:
                self.elementDict[newElement] = 1
        self.elementList = [element for element in self.elementDict.keys() if element in ['C', 'H']]
        self.elementList = sorted(self.elementList, reverse=False)
        if 'H' not in self.elementList:
            self.elementList.append('H')
            self.elementDict['H'] = 0
        otherElementList = [element for element in self.elementDict.keys() if element not in ['C', 'H']]
        otherElementList = sorted(otherElementList, reverse=False)
        self.elementList += otherElementList
        self.elementList2 = [element for element in self.elementDict.keys() if element in ['C', 'N', 'O', 'S', 'P']]

    def getElementList(self):
        return self.elementList

    def getElementList2(self):
        return self.elementList2

    def getResidueList(self):
        """
        This Function should produce a list of all natural amino acids residues that are present in the pdb file.
        Afterwards this list can be used to get the necessary restrains from ResiIns.
        Also all not natural amino acid residues are gathered into otherResiSet.
        :return: alphabetically sorted list of unique residues.
        """
        if not self.residueSet and not self.otherResiSet:
            for atom in self.atomDict.values():
                resi = atom.getResidueName().capitalize()
                if resi in AtomContainer.naturalAA:
                    self.residueSet.add(resi)
                else:
                    self.otherResiSet.add(resi.upper())
                    self.otherResiList.append(atom)
        return sorted(self.residueSet)

    def getOtherResiSet(self):
        """
        Returns a list with all residues that are not one of twenty natural amino acids.
        :return:list
        """
        return list(self.otherResiSet)

    def getOtherResiList(self):
        """
        Returns a list of all residues not including natural AA and present in the pdb file.
        """
        return self.otherResiList

    # the following part is used for solvent model analysis.

    def findAllNonAAAtoms(self):
        """
        All non natural amoni acid residue atoms are found and saved as object to a list.
        :return: list
        """
        # print 'find all non aa'
        naturalAA = [x.upper() for x in ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu',
                                         'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']]
        for atom in self.atomDict.values():
            resi = atom.getResidueName().upper()
            atomelement = atom.getAtomElement()
            # flag = atom.isComplete()[0]
            atomname = atom.getPDBAtomName()
            if resi not in naturalAA:
                # print resi, atomname
                self.nonNaturalAAAtomList.append(atom)

    def getNonNaturalAAAtoms(self):
        return self.nonNaturalAAAtomList

    def findAllWaterAtoms(self):
        """
        All atoms in water residues are saved to a list as an object.
        :return:
        """
        # print 'find all water atoms'
        for atom in self.atomDict.values():
            resi = atom.getResidueName().upper()
            atomname = atom.getAtomName()
            if resi == 'HOH':
                # print resi, atomname
                self.waterAtomList.append(atom)

    def getWaterAtomList(self):
        return self.waterAtomList

    def findSpecifiedAtom(self, type, spAtom, flag):
        """
        With given a Type of residue e.g. 'tyr' and an atom name e.g. 'O2' this function should find all residues with
        this name and perhaps an additional flag (none right now) and list
        them/return them. Atom object is writen to specifiedAtomList
        :param type:
        :param flag:
        :return: list
        """
        # print 'in findSpecAtom'
        residuetype = type
        # isComplete = flag
        specifiedAtom = spAtom
        # print residuetype, specifiedAtom
        for atom in self.atomDict.values():
            altLoc = atom.getAltLoc()
            resiname = atom.getResidueName().upper()
            atomname = atom.getPDBAtomName()
            chainID = atom.getChainID()
            resiNumber = atom.getResiSeqNum()
            # print 'ANALYSIS', resiname, atomname, chainID, resiNumber, altLoc
            # flag = atom.isComplete()[0]
            # if resiname == residuetype:
            #     print resiname, residuetype, atomname
            # if atomname == specifiedAtom:
            #     print atomname, specifiedAtom, residuetype
            if resiname == residuetype and atomname == specifiedAtom and not altLoc:
                # print 'specified atom', resiname, atomname, chainID, resiNumber
                self.specifiedAtomList.append(atom)

    def getSpecifiedAtoms(self):
        return self.specifiedAtomList

    def writeSpecifiedAtomFile(self):
        """
        All Information nescessary to identify the correct Residue in the file is written to a text file.
        :return:
        """
        specifiedAtoms = []
        for atom in self.specifiedAtomList:
            resiNumber = atom.getResiSeqNum()
            chainID = atom.getChainID()
            resiName = atom.getResidueName()
            atomName = 'CG'
            specifiedAtoms.append('{} {} {} {}'.format(atomName, resiName, chainID, resiNumber))

        f = open("SpecifiedResidues.txt", 'w')
        f.write('\n'.join(specifiedAtoms))
        f.close()

    def findNonAANearSpecAtom(self, flag):
        """
        Takes the specifiedAtomList and searches for non aa atoms near this atom. Perhaps an additional flag
        (now if it is water 'wat' or all non aa residues 'non')
        :return:
        """
        # print 'find all non aa or water near specified atom'
        if flag == 'non':
            secondList = self.getNonNaturalAAAtoms()
        else:  # also covers 'wat' as flag
            secondList = self.getWaterAtomList()

        for atom1 in secondList:
            for atom2 in self.getSpecifiedAtoms():
                if self.getAtomDistance(atom1.getAtomCoord(), atom2.getAtomCoord()) < 5:
                    self.removeAtomsList.append(atom1)
                    break

        # for i, atom1 in enumerate(self.getSpecifiedAtoms()):
        #     # altLoc1 = atom1.getAltLoc()
        #     for j in xrange(len(secondList)-i-1):
        #         j += i+1
        #         atom2 = secondList[j]
        #         # altLoc2 = atom2.getAltLoc()
        #         # print atom1.getAtomCoord(), atom2.getAtomCoord()
        #         if self.getAtomDistance(atom1.getAtomCoord(), atom2.getAtomCoord()) < 5:
        #             # self.ssBonds = True
        #             # chain1 = atom1.getChainID()
        #             # chain2 = atom2.getChainID()
        #             # num1 = atom1.getResiSeqNum()
        #             # num2 = atom2.getResiSeqNum()
        #             self.removeAtomsList.append(atom2)

    def getRemoveAtomList(self):
        return self.removeAtomsList

    def makeOMITrecords(self):
        """
        All atoms are writen to a list with the string need to omit the atom in shelxl. used here to remove atoms near
        a specified atom, to get a atom free electron density map near the specified atom.
        :return:
        """
        # print 'making omit records', self.removeAtomsList
        for atom in self.removeAtomsList:
            atomname = atom.getPDBAtomName()
            chainID = atom.getChainID()
            residueNumber = atom.getResiSeqNum()
            omitAtomString = 'OMIT ' + str(atomname) + '_' + str(chainID) + ':' + str(residueNumber).strip() + '\n'
            # print omitAtomString
            self.omitList.append(omitAtomString)

    def getOmitAtoms(self):
        return self.omitList

    # the former part is used for solvent model analysis.

    def findCysSAtoms(self):
        """
        Searches Atom Dict for sulfur atoms in Cys residues. All those atoms are appended to a List.
        :return: List of sulfur atoms in cysteine residues.
        """
        for atom in self.atomDict.values():
            resi = atom.getResidueName().upper()
            atomname = atom.getAtomElement()
            if resi == 'CYS' and atomname == 'S':
                self.cysSAtomsList.append(atom)

    def getCysSAtoms(self):
        return self.cysSAtomsList

    def getAtomDistance(self, atom1, atom2):
        """
        Takes the cartesian coordinates of two atoms and calculates their distance.
        :param atom1: np.array, three cartesian coordinates
        :param atom2: np.array, three cartesian coordinates
        :return: distance of two atoms
        """
        x1 = atom1[0]
        y1 = atom1[1]
        z1 = atom1[2]
        x2 = atom2[0]
        y2 = atom2[1]
        z2 = atom2[2]
        distance = ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
        return distance

    def findSSBonds(self):
        """
        calls upon findCysSAtoms function and creates a unique combination of all atoms in cysteine sulfur atom
        list. Of every sulfur atom pair the distance is calculated by getAtomDistance.
        If the atom distance is shorter than 2.5 Angstroms, a S-S bond is assumed and a DFIX instruction is written.
        :return:
        """
        alreadyInBondList = []
        self.findCysSAtoms()
        for i, atom1 in enumerate(self.getCysSAtoms()):
            altLoc1 = atom1.getAltLoc()
            for j in xrange(len(self.getCysSAtoms())-i-1):
                j += i+1
                atom2 = self.getCysSAtoms()[j]
                altLoc2 = atom2.getAltLoc()
                # print atom1.getAtomCoord(), atom2.getAtomCoord()
                if self.getAtomDistance(atom1.getAtomCoord(), atom2.getAtomCoord()) < 2.5:
                    # self.ssBonds = True
                    chain1 = atom1.getChainID()
                    chain2 = atom2.getChainID()
                    num1 = atom1.getResiSeqNum()
                    num2 = atom2.getResiSeqNum()
                    if chain1 == chain2 and num1 == num2:
                        pass
                    else:
                        if altLoc1 or altLoc2:
                            alreadyInBondList = self.makeSSBonds2(chain1, chain2, atom1, atom2, num1, num2, altLoc1,
                                                                  altLoc2, alreadyInBondList)
                        else:
                            alreadyInBondList = self.makeSSBonds(chain1, chain2, atom1, atom2, num1, num2,
                                                                 alreadyInBondList)
                    alreadyInBondList.append(atom1)
                    alreadyInBondList.append(atom2)

    def getSSBonds(self):
        return self.ssBondList

    def makeSSBonds2(self, chain1, chain2, atom1, atom2, num1, num2, altLoc1, altLoc2, alreadyInBondList):
        """
        All sulfur atoms with an alternate location code are handled here. These atoms will need Part instructions in
        the .ins file before and after teh atom information. Also the restraints for disulfide bridges will have to be
        addressed differently. After the chain:residuenumber a caret followed by a small letter assigning to the correct
        part is needed. e.g. 'PART 1' will need '^a' as appendix to give 'A:123^a'.
        The AltLoc given in the .pdb file can be a number or letter, so the symbol has to be assigned a fixed meaning
        for the .ins file. Numbers are kept to give the same number for the PART instruction referenced as given by
        altLocToRestraintDict for SSBONDS, letters are un-capitalized and taken for the SSBOND restraints and
        referenced for PART as before.
        It is not checked if the atom 'CB' has an altloc code and this can lead to errors while running shelxl ('No
        match for atom ... in ..').
        :param chain1:
        :param chain2:
        :param atom1:
        :param atom2:
        :param num1:
        :param num2:
        :param altLoc1:
        :param altLoc2:
        :return:
        """
        # altLocToRestraintDict = {'1': 'a', '2': 'b', '3': 'c', '4': 'd', '5': 'e', '6': 'f', '7': 'g', '8': 'h',
        #                          '9': 'i', '10': 'j', '11': 'k', '12': 'l', '13': 'm', '14': 'n', '15': 'o', '16': 'p',
        #                          '17': 'q', '18': 'r', '19': 's', '20': 't', '21': 'u', '22': 'v', '23': 'w', '24': 'x',
        #                          '25': 'y', '26': 'z'}
        altLocToRestraintDict = {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f', 7: 'g', 8: 'h',
                                 9: 'i', 10: 'j', 11: 'k', 12: 'l', 13: 'm', 14: 'n', 15: 'o', 16: 'p',
                                 17: 'q', 18: 'r', 19: 's', 20: 't', 21: 'u', 22: 'v', 23: 'w', 24: 'x',
                                 25: 'y', 26: 'z'}
        # altLoc1 = str(altLoc1)
        # altLoc2 = str(altLoc2)
        try:
            if altLoc1:
                altAltLoc1 = self.altLocDict[altLoc1]
            else:
                altAltLoc1 = None
        except KeyError:
            try:
                altLoc1 = int(altLoc1)
                altAltLoc1 = altLocToRestraintDict[altLoc1]  # the number gets the corresponding letter assigned
                self.altLocDict[altLoc1] = altAltLoc1  # the match is saved to the altLocDict
            except ValueError:
                altLoc1 = str(altLoc1)
                if altLoc1.lower() != altLoc1:
                    altAltLoc1 = altLoc1.lower()
                    self.altLocDict[altLoc1] = altAltLoc1  # the upper case letter has reference to the needed letter.
                # altLoc1 = ord(altLoc1.lower())-96  # the letter gets the appropriate number assigned
                else:
                    altAltLoc1 = altLoc1
                    self.altLocDict[altLoc1] = altAltLoc1  # here the altAltLoc and altLoc are the same, but no KeyError
        try:
            if altLoc2:
                altAltLoc2 = self.altLocDict[altLoc2]
            else:
                altAltLoc2 = None
        except KeyError:
            try:
                altLoc2 = int(altLoc2)
                altAltLoc2 = altLocToRestraintDict[altLoc2]  # the number gets the corresponding letter assigned
                self.altLocDict[altLoc2] = altAltLoc2  # the match is saved to the altLocDict
            except ValueError:
                altLoc2 = str(altLoc2)
                if altLoc2.lower() != altLoc1:
                    altAltLoc1 = altLoc2.lower()
                    self.altLocDict[altLoc2] = altAltLoc1  # the upper case letter has reference to the needed letter.
                # altLoc1 = ord(altLoc1.lower())-96  # the letter gets the appropriate number assigned
                else:
                    altAltLoc1 = altLoc1
                    self.altLocDict[
                        altLoc1] = altAltLoc1  # here the altAltLoc and altLoc are the same, but no KeyError
            # if altLoc1.isdigit():
            #     altAltLoc1 = altLocToRestraintDict[altLoc1]
            #     self.altLocDict[altLoc1] = altAltLoc1
            # elif altLoc1.isalpha():
            #     for key, value in altLocToRestraintDict:
            #         if value == altLoc1.lower():
            #             newAltLoc = key  # here should the PART number be assigned.
            #             if altLoc1.lower() != altLoc1:
            #                 altAltLoc1 = altLoc1.lower()  # altAltLoc should be lower case
            #                 self.altLocDict[altLoc1] = altAltLoc1  # now the lower case can be directly assigned in try.
            #                 self.altLocDict[altAltLoc1] = newAltLoc
            #             else:
            #                 self.altLocDict[altLoc1] = newAltLoc
            #                 altAltLoc1 = altLoc1
            # else:
            #     self.makeSSBonds(chain1, chain2, atom1, atom2, num1, num2, alreadyInBondList)
        # try:
        #     altAltLoc2 = self.altLocDict[altLoc2]
        # except KeyError:
        #     if altLoc2.isdigit():
        #         altAltLoc2 = altLocToRestraintDict[altLoc2]
        #         self.altLocDict[altLoc2] = altAltLoc2
        #     elif altLoc2.isalpha():
        #         for key, value in altLocToRestraintDict:
        #             if value == altLoc2.lower():
        #                 newAltLoc = key  # here should the PART number be assigned.
        #                 if altLoc2.lower() != altLoc1:
        #                     altAltLoc2 = altLoc2.lower()  # altAltLoc should be lower case
        #                     self.altLocDict[altLoc2] = altAltLoc2  # now the lower case can be directly assigned in try.
        #                     self.altLocDict[altAltLoc2] = newAltLoc
        #                 else:
        #                     self.altLocDict[altLoc2] = newAltLoc
        #                     altAltLoc2 =altLoc2
        #     else:
        #         print 'INFO: Problem handling alternate location code during creation of disulfide bond restraints.'
        #         self.makeSSBonds(chain1, chain2, atom1, atom2, num1, num2, alreadyInBondList)
        if altAltLoc1 and altAltLoc2:
            if atom1 not in alreadyInBondList and atom2 not in alreadyInBondList:
                self.ssBondList.append('\nDFIX 2.031 SG_{}:{}^{} SG_{}:{}^{}'
                                       .format(chain1, num1.strip(), altAltLoc1, chain2, num2.strip(), altAltLoc2))
                self.ssBondList.append('\nDANG 3.035 SG_{}:{}^{} CB_{}:{}^{}'
                                       .format(chain1, num1.strip(), altAltLoc1, chain2, num2.strip(), altAltLoc2))
                self.ssBondList.append('\nDANG 3.035 SG_{}:{}^{} CB_{}:{}^{}'
                                       .format(chain2, num2.strip(), altAltLoc2, chain1, num1.strip(), altAltLoc1))
                alreadyInBondList.append(atom1)
                alreadyInBondList.append(atom2)
                self.ssBonds = True
        elif altAltLoc1 and not altAltLoc2:
            if atom1 not in alreadyInBondList and atom2 not in alreadyInBondList:
                self.ssBondList.append('\nDFIX 2.031 SG_{}:{}^{} SG_{}:{}'
                                       .format(chain1, num1.strip(), altAltLoc1, chain2, num2.strip()))
                self.ssBondList.append('\nDANG 3.035 SG_{}:{}^{} CB_{}:{}'
                                       .format(chain1, num1.strip(), altAltLoc1, chain2, num2.strip()))
                self.ssBondList.append('\nDANG 3.035 SG_{}:{} CB_{}:{}^{}'
                                       .format(chain2, num2.strip(), chain1, num1.strip(), altAltLoc1))
                alreadyInBondList.append(atom1)
                alreadyInBondList.append(atom2)
                self.ssBonds = True
        elif not altAltLoc1 and altAltLoc2:
            if atom1 not in alreadyInBondList and atom2 not in alreadyInBondList:
                self.ssBondList.append('\nDFIX 2.031 SG_{}:{} SG_{}:{}^{}'
                                       .format(chain1, num1.strip(), chain2, num2.strip(), altAltLoc2))
                self.ssBondList.append('\nDANG 3.035 SG_{}:{} CB_{}:{}^{}'
                                       .format(chain1, num1.strip(), chain2, num2.strip(), altAltLoc2))
                self.ssBondList.append('\nDANG 3.035 SG_{}:{}^{} CB_{}:{}'
                                       .format(chain2, num2.strip(), altAltLoc2, chain1, num1.strip()))
                alreadyInBondList.append(atom1)
                alreadyInBondList.append(atom2)
                self.ssBonds = True
        return alreadyInBondList

    def makeSSBonds(self, chain1, chain2, atom1, atom2, num1, num2, alreadyInBondList):
        if atom1 not in alreadyInBondList and atom2 not in alreadyInBondList:
            self.ssBondList.append('\nDFIX 2.031 SG_{0}:{1:{4}f} SG_{2}:{3:{4}f}'
                                   .format(chain1, float(num1), chain2, float(num2), padding))
            self.ssBondList.append('\nDANG 3.035 SG_{0}:{1:{4}f} CB_{2}:{3:{4}f}'
                                   .format(chain1, float(num1), chain2, float(num2), padding))
            self.ssBondList.append('\nDANG 3.035 SG_{0}:{1:{4}f} CB_{2}:{3:{4}f}'
                                   .format(chain2, float(num2), chain1, float(num1), padding))
            alreadyInBondList.append(atom1)
            alreadyInBondList.append(atom2)
            self.ssBonds = True
        return alreadyInBondList

    def asShelxString(self, cell):
        """
        Takes all objects in Atomdict and sorts them according to chain ID and Residue number and priority.
        Since neg residue numbers were introduced, the sorting has to be done trice.
        But priority is only applied if no H or D atoms are given in natural aa residues. otherwise priority is not used
        to sort the atoms in one residue according to a given template.
        AtomObjectList is searched for the start of a new Residue and a new line RESI is introduced before.
        The RESI line gives a :RESI: statement followed by the Chain ID number for shelx : chain ID character is
         transformed into a number (first number of four) followed by the residue number
         (i.e. residue chain A residue number 23 yields 1023).
        The last element of the RESI line is the residue name; three letters representing the residues amino acid.
        chainID 'XXX' = offset for chains with more than 999 residues. (old)
        chainID 'ZZZ' = water residues (old)
        The function getLastAtomInResi(Residue) is called to find incomplete residues. those get a HFIX 0 instruction
        for the last present atom before the missing atom(s). this instruction is writen into the incompleteResiString.
        :return: All atoms sorted by chains and residues to give a string.
        """
        # print type(atom.getChainID()), type(atom.getResiSeqNum())
        # print self.neut
        if self.neut or self.keepHAtoms:
            atomObjectList = sorted(self.atomDict.values(), key=lambda atom: (atom.getChainID()+atom.getResiSeqNum()))
        else:
            atomObjectList = sorted(sorted(sorted([aO for aO in self.atomDict.values() if not aO.getAtomElement() == 'H'], key=lambda atom: (atom.getPriority()))
                                           , key=lambda atom: (atom.getResiSeqNumAsInt())), key=lambda atom: (atom.getChainID()))
        # maybe change function above to ignoring not only H atoms, but also D atoms?
        atomStringList = []
        makeHFIXfor = None
        residueBefore = None
        beware = False
        chainIDbefore = None
        waterResiCounter = 0
        # cTermOxList = ['OT', 'OXT', 'OX', 'OX1', 'OX2']
        # for chain in self.chains.values():
        #     for residue in chain.values():
        #         for atom in residue:
        #             if atom.getPDBAtomName() in cTermOxList:
        #                 residue.renameCtermOxygen()
        for atom in atomObjectList:
            # print atom.getChainID(), atom.getResidueName(), atom.getResiSeqNum()
            # if atom.getChainID() == 'HOH' and resetOccupancy:
            #     atom.resetOccupancy(1.00000)
            residueNew = atom.getResiSeqNum()
            # chainID = atom.getChainID()
            longchain = False
            water = False
            chainIDSet = self.chainIDSet
            chainID = atom.getChainID()
            # try:
            #     chainIDSet.remove('XXX')
            #     self.ligandChain = True
            # except KeyError:
            #     pass
            try:
                chainIDSet.remove('ZZZ')
                water = True
            except KeyError:
                pass
            # if len(chainID) == int(1):
            #     if ord(chainID.lower())-96 >= int(10):
            #         beware = True
            # if chainID == 'XXX' or chainID == 'ZZZ':
            #     if beware:
            #         if chainID == 'XXX':
            #             self.ligandChain = True
            #             atom.setChainID(chr(len(self.chainIDSet)+66))
            #         if chainID == 'ZZZ' and self.ligandChain:
            #             water = True
            #             atom.setChainID(chr(len(self.chainIDSet)+67))
            #         if chainID == 'ZZZ' and not self.ligandChain:
            #             water = True
            #             atom.setChainID(chr(len(self.chainIDSet)+66))
            #     else:
            #         if chainID == 'XXX':
            #             self.ligandChain = True
            #             atom.setChainID(chr(len(self.chainIDSet)+66))
            #         if chainID == 'ZZZ' and self.ligandChain:
            #             water = True
            #             atom.setChainID(chr(len(self.chainIDSet)+67))

            # This part only needed for HOH in separate chain, not needed if HOH not treated in special way.
            # if chainID == 'ZZZ':  # and not self.ligandChain:
            #     water = True
            #     atom.setChainID(chr(len(self.chainIDSet)+66))  # can create problems if file has more than 26 chains?

            # This part is older!
            # else:
            #     if ord(chainID.lower())-96 >= int(10) and len(chainIDSet) == int(1):
            #         atom.setChainID(chr(len(chainIDSet)+64))
            #     if chainID in self.overlongChainList and chainID == 'A':
            #         while self.printWarning:
            #             print 'The overlong chain A will be without offset.'
            #             self.printWarning = False
            #         longchain = True

            # This part only needed for HOH in separate chain, not needed if HOH not treated in special way.
            # if not residueNew == residueBefore or not chainID == chainIDbefore:
            #     if not atom.getResidueName() == 'HOH':
            #         # if longchain:
            #         #     atomStringList.append("\nRESI {} {}\n".format('{:0>3.0f}'.format(float(residueNew)),
            #         #                                                   atom.getResidueName()))
            #         if atom.getChainID() in self.chainIDSet:
            #             try:
            #                 makeHFIXfor = self.chains[atom.getChainID()][atom.getResiSeqNum()].getLastAtomInResi()
            #             except NoResidueError:
            #                 pass
            #         # else:
            #         atomStringList.append("\nRESI {} {}\n".format('{0}:{1:{2}f}'.format(atom.getChainID(),
            #                                                                             float(residueNew), padding),
            #                                                       atom.getResidueName()))
            #     else:
            #         # print 'water Resi found with new chain ID and string written: ', waterResiCounter, residueNew
            #         waterResiCounter += 1
            #         waterResiNum = waterResiCounter
            #         self.waterResiNumbers.append(waterResiNum)
            #         atomStringList.append("\nRESI {} {}\n".format('{0}:{1:{2}f}'.format(atom.getChainID(),
            #                                                                             float(waterResiNum), padding),
            #                                                       atom.getResidueName()))
            #         # print "\nRESI {} {}\n".format('{0}:{1:{2}f}'.format(atom.getChainID(), float(waterResiNum),
            #         #                                                     padding), atom.getResidueName())

            # This part is needed if HOH is not treated in a special way in an extra chain:
            if not residueNew == residueBefore or not chainID == chainIDbefore:
                # if longchain:
                #     atomStringList.append("\nRESI {} {}\n".format('{:0>3.0f}'.format(float(residueNew)),
                #                                                   atom.getResidueName()))
                if atom.getChainID() in self.chainIDSet and not atom.getResidueName() == 'HOH':
                    try:
                        makeHFIXfor = self.chains[atom.getChainID()][atom.getResiSeqNum()].getLastAtomInResi()
                    except NoResidueError:
                        pass
                # else:
                atomStringList.append("\nRESI {} {}\n".format('{0}:{1:{2}f}'.format(atom.getChainID(),
                                                                                    float(residueNew), padding),
                                                              atom.getResidueName()))

            if makeHFIXfor:
                if atom.getPDBAtomName() in makeHFIXfor:
                    self.incompleteResiString.append('REM HFIX 0 {}_{}\n'.format(atom.getPDBAtomName(),
                                                     '{0}:{1:{2}f}'.format(atom.getChainID(), float(residueNew),
                                                                           padding)))

            # This part is older!
            # if chainID == chainIDbefore and atom.getResidueName() == 'HOH':
            #     print 'water resi found with same chain ID and string written: ', waterResiCounter, residueNew
            #     waterResiCounter += 1
            #     waterResiNum = waterResiCounter
            #     self.waterResiNumbers.append(waterResiNum)
            #     atomStringList.append("\nRESI {} {}\n".format('{0}:{1:{2}f}'
            # .format(atom.getChainID(), float(waterResiNum), padding),
            #                                                   atom.getResidueName()))
            #     print "\nRESI {} {}\n".format('{0}:{1:{2}f}'.format(atom.getChainID(), float(waterResiNum), padding),
            #                                   atom.getResidueName())

            atomStringList.append(atom.asShelxString(self.elementList, cell))
            chainIDbefore = chainID
            residueBefore = str(residueNew)
        # atomStringList = [atom.asShelxString(self.elementList, cell) for atom in atomObjectList]
        return ''.join(atomStringList)


class Atom(object):
    forcedOccupancies = {}

    def __init__(self, line, *args, **kwargs):
        self.atomName = None
        self.isAnisou = False
        self.atomAnisou = None
        self.line = line
        self.PDBAtomName = None
        self.atomCoord = None
        self.atomCoordFrac = None
        self.atomOccupancy = None
        self.atomElement = None
        self.residueName = None
        self.atomSFAC = None
        self.chainID = None
        self.atomUIso = None
        self.altLoc = None
        self.insertionCode = None
        self.atomTempFactor = None
        self.createAtomName()
        self.residueList = None
        self.resiSeqNum = None
        self.atomString = None
        self.hasWater = False
        self.resiSeqOffsetDict = {}
        if line.startswith("ANISOU"):
            self.extractAnisou()
        else:
            self.extractAtomCoord()
            self.extractAtomElement()
            self.extractOccupancy()
            self.extractTempFactor()
            self.extractResidueName()
            self.extractChainID()
            self.extractResiSeqNum(self.getChainID())
            self.extractAltLoc()
            self.extractInsertionCode()

    def createAtomName(self):
        """
        Creates an unique name for each atom in the pdb file:
        takes the PDB atom name, PDB atom serial number, PDB residue name and PDB chain ID to create an atom name.
        :return: unique Atom name created from PDB data.
        """
        self.PDBAtomName = self.line[12:16].replace(' ', '')  # atom name in pdb
        while self.residueName == 'ILE':
            if self.PDBAtomName == 'CD':
                self.PDBAtomName = 'CD1'
            if self.PDBAtomName == 'OXT':
                self.PDBAtomName = 'OT2'
        atomSerial = self.line[6:11]  # atom serial number in pdb
        self.atomName = "{}_{}_{}_{}".format(self.PDBAtomName, atomSerial, self.getResidueName(), self.getChainID())

    def getAtomName(self):
        return self.atomName

    def getPDBAtomName(self):
        return self.PDBAtomName

    def getPriority(self):
        """
        returns the priority of an atom based on its residue and PDB atomname.
        The residueAtomDict in Residue class is used to assign priorities.
        Also the C terminal Oxygen atom name has to be handled.
        :return:
        """
        try:
            return '{:0>2}'.format(Residue.residueAtomDict2[self.residueName].index(self.getPDBAtomName()))
        except KeyError:
            return '00'
        except ValueError:
            cTermOxList = ['OT', 'OXT', 'OX', 'OX1', 'OX2']
            atomname = self.getPDBAtomName()
            if self.getPDBAtomName() in cTermOxList:
                if atomname == 'O' or atomname == 'OX' or atomname == 'OX1':
                    self.rename('OT1')
                elif atomname == 'OXT' or atomname == 'OT' or atomname == 'OX2':
                    self.rename('OT2')
            if self.getPDBAtomName() == 'OT2' or self.getPDBAtomName() == 'OT1':
                return '13'
            else:
                print ' *** ERROR: Illegal atom name for atom {} in ' \
                      'residue {}:{} {} ***'.format(self.getPDBAtomName(), self.getChainID(),
                                                    self.getResiSeqNum().lstrip(), self.residueName)
                # print self.getAtomElement(), self.getPDBAtomName(), '{}:{}'.format(self.getChainID(),
                #                                                                    self.getResiSeqNum().lstrip()), \
                #     self.residueName
                # raise ValueError
                print '*** PDB2INS is terminated without writing an .ins file. ***'
                exit()

    def rename(self, newName):
        """
        this function overrides the PDBAtomName with a new atom name.
        :param newName:
        :return:
        """
        self.PDBAtomName = newName

    def extractChainID(self):
        """
        extracts the chain identifier for the atom from the line of the pdb file.
        Water molecules are exracted separately.
        :return: chain ID (a character in pdb file)
        """
        self.chainID = self.line[21]
        # This part is only needed if water residues are treated in special way (in different chain).
        # if not self.residueName == "HOH":
        #     self.chainID = self.line[21]
        #     # self.chainID = self.line[21].upper()
        # else:
        #     self.hasWater = True
        #     self.chainID = 'ZZZ'

    def getChainID(self):
        return self.chainID

    def setChainID(self, ID):
        self.chainID = ID

    def extractAtomCoord(self):
        """
        extracts the atoms coordinates from the pdb file line as float.
        The coordinates are given in cartesian coordinates from pdb and transformed into fractional coordinates.
        :return: fractional atom coordinates as array.
        """
        x_coord = float(self.line[30:38])
        y_coord = float(self.line[38:46])
        z_coord = float(self.line[46:54])
        self.atomCoord = np.array((x_coord, y_coord, z_coord))

    def getAtomCoord(self):
        return self.atomCoord

    def extractAtomElement(self):
        """
        The atom's element is extracted from the line, all spaces are stripped. pdb gives the element as string.
        :return:atom element.
        """
        self.atomElement = self.line[76:78].replace(' ', '')
        # print self.atomElement
        if not self.atomElement:
            self.atomElement = self.line[12:16].replace(' ', '')[0]
        # print 'new', self.atomElement

    def getAtomElement(self):
        return self.atomElement

    def extractOccupancy(self):
        """
        The atom's occupancy is extracted from the line, all spaces are stripped and the occupancy is returned as float.
        the pdb file gives the atom's occupancy as real 6.2.
        :return: atom occupancy as float.
        """
        self.atomOccupancy = float(self.line[54:60].replace(' ', ''))

    def getOccupancy(self):
        """
        this function can override the atoms occupancies from the pdb file if the another occupancy is given in
        resetOccupancy.
        :return: float (occupancy)
        """
        try:
            occ = Atom.forcedOccupancies[self.residueName]
        except KeyError:
            occ = self.atomOccupancy
        return occ

    @staticmethod
    def resetOccupancy(resiname, newOccupancy):
        Atom.forcedOccupancies[resiname] = newOccupancy

    def extractTempFactor(self):
        """
        The atoms's temperature factor is extracted from the line and returned as float.
        The temperature factor is given as real 6.3 in the pdb file.
        :return: temperature factor as float.
        """
        self.atomTempFactor = float(self.line[60:66])

    def getTempFactor(self):
        return self.atomTempFactor

    def getUIso(self):
        """
        Calculates Uiso from the temperature factor in the pdb line.
        :return: float, uiso
        """
        self.atomUIso = (self.atomTempFactor / (8 * np.pi * np.pi))
        return float(self.atomUIso)

    def getUanis(self, cell):
        """
        The u anisou are transforms to fractional coordinates using transformation.py
        :param cell: Needed to transform the ADP.
        :return: fractional ADP in the correct SHELX order.
        """
        adp = transformations.cart2frac_ADP(self.atomTempFactor, cell)
        return [adp[0], adp[1], adp[2], adp[5], adp[4], adp[3]]

    def getTempFactorAsString(self, cell):
        """
        A string for U iso or anisou is produced.
        :param cell: Needed for transformation to fractional coordinates.
        :return: string of Uiso or Uanisou
        """
        if self.isAnisou:
            string = '{:6.5f} {:6.5f}=\n    {:6.5f} {:6.5f} {:6.5f} {:6.5f}'.format(*self.getUanis(cell))
        else:
            string = '{:6.5f}'.format(self.getUIso())
        return string

    def extractAnisou(self):
        """
        extracts the atoms anisotropic temperature factors for each atom from the ANISOU line of the pdb file.
        The Atom anisou data has to be corrected by 10-4 to give the correct magnitude.
        :return: array of the atoms's anisou temp factors.
        """
        atomAnisou11 = float(self.line[28:35])
        atomAnisou22 = float(self.line[35:42])
        atomAnisou33 = float(self.line[42:49])
        atomAnisou12 = float(self.line[49:56])
        atomAnisou13 = float(self.line[56:63])
        atomAnisou23 = float(self.line[63:70])
        self.atomAnisou = np.array((atomAnisou11, atomAnisou22, atomAnisou33, atomAnisou12, atomAnisou13, atomAnisou23))
        self.atomAnisou = self.atomAnisou * 0.0001

    def getRawADP(self):
        return self.atomAnisou

    def extractInsertionCode(self):
        try:
            self.insertionCode = self.line[26]
            if self.insertionCode == ' ':
                self.insertionCode = None
        except KeyError:
            pass

    def getInsertionCode(self):
        return self.insertionCode

    def extractAltLoc(self):
        """
        For every atom line the alternate location indicator is extracted. If none is given at the first position (17),
        the second position it was found during tests is used (27, code for insertion of residues).
        The alternate altloc is changed to a number representing the character if necessary.
        :return:
        """
        # translatedict = {'A': 1, 'B': 2, 'C': 3}
        try:
            self.altLoc = self.line[16]
            if self.altLoc == ' ':
                self.altLoc = None
                # try:
                #     self.altLoc = self.line[26]
                #     if self.altLoc == ' ':
                #
                # except KeyError:
                #     pass
        except KeyError:
            pass
        if self.altLoc:
            try:
                self.altLoc = int(self.altLoc)
            except ValueError:
                self.altLoc = ord(self.altLoc.lower())-96  # transforms letter to integer
                # except AttributeError:
                #     self.altLoc = None

    def getAltLoc(self):
        return self.altLoc

    def extractResidueName(self):
        """
        The residue name is extracted from the atom's pdb line.
        All water molecules are named HOH.
        :return: residue name as string
        """
        waterList = ['H2O', 'OH2', 'WAT']
        self.residueName = self.line[17:20]
        if self.residueName in waterList:
            self.residueName = 'HOH'

    def getResidueName(self):
        return self.residueName

    def extractResiSeqNum(self, chainID):
        """
        Extracts the residue sequence number from the atom's pdb line.
        :return: residue sequence number as string.
        """
        self.resiSeqNum = self.line[22:26]

    def setResiSeqNum(self, newNum):
        self.resiSeqNum = newNum

    def getResiSeqNum(self):
        return '{:>4}'.format(self.resiSeqNum)

    def getResiSeqNumAsInt(self):
        return int(self.resiSeqNum)

    def getAtomSFAC(self, elementList):
        """
        takes the elements collected in element list and
        assigns them a position number for the 'SFAC' line of the .ins file.
        :param elementList: List of all elements in the pdb file
        :return: the elements sfac number, starting with 1
        """
        self.atomSFAC = elementList.index(self.atomElement)+1

    def overrideADP(self, newADP):
        self.atomTempFactor = newADP
        self.isAnisou = True

    def asShelxString(self, elementList, cell):
        """
        Takes the cell coordinates of the molecule and transforms the atom coordinates
        from cartesian to fractional coordinates.
        If the atom has an altloc (alternate location indicator), the residues are divided into parts and all atoms with
        altloc get a part number coresponding to the initial altloc extracted fromt he pdb file.
        Number from 1 to 99 is legal.
        Also the atom string for the .ins file is created, containing the atoms PDB atom name, th   e atom's SFAC,
        the atom's fractional coordinates and the atom's occupancy with a preceding number one.
        :param elementList:  List of all elements in the pdb file
        :param cell: the molecules cell coordinates in fractional coordinates
        :return: string for one line in .ins file for each atom.
        """
        self.getAtomSFAC(elementList)
        self.atomCoordFrac = transformations.cart2frac(self.atomCoord, cell)
        tempFactor = self.getTempFactorAsString(cell)
        usedParts = []
        alternativeParts = {}
        if self.getAltLoc():
            altLoc = self.getAltLoc()
            if 1 <= altLoc <= 99:
                usedParts.append(altLoc)
            # if altLoc == 'A':
            #     self.atomString = "{:<6} {}  {}  1{:6.5f}  {}\n".format(self.PDBAtomName, self.atomSFAC,
            #                                                             '{:7.6f}  {:7.6f}  {:7.6f}'
            #                                                             .format(self.atomCoordFrac[0],
            #                                                                     self.atomCoordFrac[1],
            #                                                                     self.atomCoordFrac[2]),
            #                                                             self.getOccupancy(), tempFactor)
            # else:
                partBegin = altLoc
                # print partBegin
            else:
                try:
                    partBegin = alternativeParts[altLoc]
                except KeyError:
                    sortedUsedParts = sorted(usedParts)
                    partBegin = sortedUsedParts[-1] + 1
                    usedParts.append(partBegin)
                    alternativeParts[altLoc] = usedParts
            partEnd = 0
            self.atomString = "PART {}\n".format(partBegin)
            self.atomString += "{:<6} {}  {}  1{:6.5f}  {}\n".format(self.PDBAtomName, self.atomSFAC,
                                                                     '{:7.6f}  {:7.6f}  {:7.6f}'
                                                                     .format(self.atomCoordFrac[0],
                                                                             self.atomCoordFrac[1],
                                                                             self.atomCoordFrac[2]),
                                                                     self.getOccupancy(), tempFactor)
            self.atomString += "PART {}\n".format(partEnd)
        else:
            self.atomString = "{:<6} {}  {}  1{:6.5f}  {}\n".format(self.PDBAtomName, self.atomSFAC,
                                                                    '{:7.6f}  {:7.6f}  {:7.6f}'
                                                                    .format(self.atomCoordFrac[0],
                                                                            self.atomCoordFrac[1],
                                                                            self.atomCoordFrac[2]), self.getOccupancy(),
                                                                    tempFactor)
        return self.atomString

    # def registerResidue(self, residue):
    #     self.residue=residue
    #
    # def getResidueNumber(self):
    #     return self.residue.getMyNumber()


class Chain(dict):

    def __init__(self, *args, **kwargs):
        super(Chain, self).__init__(*args, **kwargs)
        self.resiChain = None
        self.residueList = []

    def finalize(self):
        """
        the dictionary resichain is sorted by keys and returned. lambda function to ensure the residues are sorted after
        numbers as integers and not string. String sorting would give false c terminal residues in getCTernResi()!
        :return: dict
        """
        self.resiChain = sorted(self.keys(), key=lambda x: int(x), reverse=False)

    def iterate(self):
        """
        the dictionary resichain is iterated by item and each item returned.
        :return:
        """
        # for key in self.resiChain:
        #     print key, self.resiChain[key]
        return [(item, self[item]) for item in self.resiChain]

    def __getitem__(self, item):
        try:
            return super(Chain, self).__getitem__(item)
        except KeyError:
            raise NoResidueError

    # def getChainID(self):
    #     return self.chainID

    def getNTermResi(self):
        """
        checks if the residues are natural amino acids. returns the amino acid only if true. the residues are iterated
        in normal order, so that the first residue from the N terminal side for which this criteria is true is found.
        :return:
        """
        # self.finalize()
        for resiNumber, residue in self.iterate():
            if residue.isAminoAcid():
                # print 'This is the resinumber', resiNumber, 'and this is the residue', residue
                return residue

    def getResiList(self):
        """
        Should return a list containing the residue sequence for the chain.
        :return:
        """
        for resiNumber, residue in self.iterate():
            self.residueList.append(residue.getResiName())
        return self.residueList

    def getCTermResi(self):
        """
        checks if the residues are natural amino acids. Returns the residue only if true. The residues are iterated in
        reverse, so that the first residue from the c terminal side for which this criteria is true is found.
        :return: string (residue)
        """
        # self.finalize()
        for resiNumber, residue in reversed(self.iterate()):
            if residue.isAminoAcid() and not residue.hasInsertionCode():
                return residue

    def getChainLength(self):
        return len(self)

    # def getMyNumber(self, residueID):
    #     self.container


class Residue(list):

    naturalAA = [x.upper() for x in ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys',
                                     'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']]
    naturalAAAtomNumbers = {'ALA': 6, 'ARG': 12, 'ASN': 9, 'ASP': 9, 'CYS': 7, 'GLN': 10, 'GLU': 10, 'GLY': 5, 'HIS': 10,
                            'ILE': 9, 'LEU': 9, 'LYS': 9, 'MET': 9, 'PHE': 12, 'PRO': 8, 'SER': 7, 'THR': 8, 'TRP': 15,
                            'TYR': 13, 'VAL': 8}
    residueAtomDict = {'ALA': ['C', 'CA', 'N', 'O', 'CB'],
                       'ARG': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                       'ASN': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'OD1', 'ND2'],
                       'ASP': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'OD1', 'OD2'],
                       'CYS': ['C', 'CA', 'N', 'O', 'CB', 'SG'],
                       'GLN': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
                       'GLU': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
                       'GLY': ['C', 'CA', 'N', 'O'],
                       'HIS': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                       'ILE': ['C', 'CA', 'N', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
                       'LEU': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD1', 'CD2'],
                       'LYS': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
                       'MET': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'SD', 'CE'],
                       'PHE': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                       'PRO': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD'],
                       'SER': ['C', 'CA', 'N', 'O', 'CB', 'OG'],
                       'THR': ['C', 'CA', 'N', 'O', 'CB', 'OG1', 'CG2'],
                       'TRP': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3'],
                       'TYR': ['C', 'CA', 'N', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
                       'VAL': ['C', 'CA', 'N', 'O', 'CB', 'CG1', 'CG2']
                       }
    residueAtomDict2 = {'ALA': ['N', 'CA', 'CB', 'C', 'O'],
                        'ARG': ['N', 'CA', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'C', 'O'],
                        'ASN': ['N', 'CA', 'CB', 'CG', 'OD1', 'ND2', 'C', 'O'],
                        'ASP': ['N', 'CA', 'CB', 'CG', 'OD1', 'OD2', 'C', 'O'],
                        'CYS': ['N', 'CA', 'CB', 'SG', 'C', 'O'],
                        'GLN': ['N', 'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'C', 'O'],
                        'GLU': ['N', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'C', 'O'],
                        'GLY': ['N', 'CA', 'C', 'O'],
                        'HIS': ['N', 'CA', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', 'C', 'O'],
                        'ILE': ['N', 'CA', 'CB', 'CG1', 'CG2', 'CD1', 'C', 'O'],
                        'LEU': ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'C', 'O'],
                        'LYS': ['N', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'C', 'O'],
                        'MET': ['N', 'CA', 'CB', 'CG', 'SD', 'CE', 'C', 'O'],
                        'PHE': ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'C', 'O'],
                        'PRO': ['N', 'CA', 'CB', 'CG', 'CD', 'C', 'O'],
                        'SER': ['N', 'CA', 'CB', 'OG', 'C', 'O'],
                        'THR': ['N', 'CA', 'CB', 'OG1', 'CG2', 'C', 'O'],
                        'TRP': ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3', 'C', 'O'],
                        'TYR': ['N', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'C', 'O'],
                        'VAL': ['N', 'CA', 'CB', 'CG1', 'CG2', 'C', 'O']
                        }
    resiAtomTupleDict = {'ALA': [('CA', 'CB')],
                         'ARG': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'NE'), ('NE', 'CZ'), ('CZ', 'NH1'),
                                 ('CZ', 'NH2')],
                         'ASN': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'OD1'), ('CG', 'NG2')],
                         'ASP': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'OD1'), ('CG', 'OD2')],
                         'CYS': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'SD')],
                         'GLN': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'OE1'), ('CD', 'NE2')],
                         'GLU': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'OE1'), ('CD', 'OE2')],
                         'GLY': [],
                         'HIS': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'ND1'), ('CG', 'CD2'), ('ND1', 'CE1'),
                                 ('CD2', 'NE2'), ('NE2', 'CE1'), ('CE1', 'NE2')],
                         'ILE': [('CA', 'CB'), ('CB', 'CG1'), ('CG1', 'CD1')],
                         'LEU': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2')],
                         'LYS': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'CE'), ('CE', 'NZ')],
                         'MET': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'SD'), ('SD', 'CE')],
                         'PHE': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'CE1'),
                                 ('CD2', 'CE2'), ('CE1', 'CZ'), ('CE2', 'CZ')],
                         'PRO': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD')],
                         'SER': [('CA', 'CB'), ('CB', 'OG')],
                         'THR': [('CA', 'CB'), ('CB', 'OG1'), ('CB', 'CG2')],
                         'TRP': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2'), ('CG1', 'NE1'),
                                 ('NE1', 'CE2'), ('CD2', 'CE2'), ('CD2', 'CE3'), ('CE2', 'CZ2'), ('CZ2', 'CH2'),
                                 ('CH2', 'CZ3'), ('CZ3', 'CH2')],
                         'TYR': [('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2'), ('CD1', 'CE1'),
                                 ('CD2', 'CE2'), ('CE1', 'CZ'), ('CE2', 'CZ'), ('CZ', 'OH')],
                         'VAL': [('CA', 'CB'), ('CB', 'CG1'), ('CB', 'CG2')]
                         }

    def getAtomNames(self):
        atomnames = []
        for atom in self:
            atomnames.append(atom.getPDBAtomName())
        return atomnames

    def renameCtermOxygen(self):
        """
        Now all O atoms called 'OXT' are automatically renamed to 'OT1' or 'OT2' when they fit the naming scheme for
        c terminal oxygen atoms.
        It is assumed that all oxygen atoms called OXT, OT, OX1, OX2, OX are c terminal oxygen atoms.
        ***This function should be superfluous after changes in PRIORITY ASSIGNMENT.***
        :return:
        """
        for atom in self:
            atomname = atom.getPDBAtomName()
            if atomname == 'O' or atomname == 'OX' or atomname == 'OX1':
                atom.rename('OT1')
            elif atomname == 'OXT' or atomname == 'OT' or atomname == 'OX2':
                atom.rename('OT2')
            else:
                pass

    def getMissingAtoms(self):
        """
        Checks the residues for missing atoms. The missing atom names are found by subtracting the present atoms for
        a specific residue from the expected atoms given in the dictionary resiAtomDict (specific entries for each
        general amino acid).
        :return:Boolean (residueIncomplete), list (missingAtomNames)
        """
        residueIncomplete = False
        residueName = self[0].getResidueName()
        try:
            missingAtomNames = set(self.residueAtomDict[residueName]) - set(self.getAtomNames())
        except KeyError:
            return residueIncomplete, []
        else:
            if missingAtomNames:
                residueIncomplete = True
            return residueIncomplete, missingAtomNames

    def getLastAtomInResi(self):
        """
        mAN = missingAtomNames from the function getMissingAtoms.
        Finds the last not missing atom in an incomplete residue by checking against the resiAtomTupleDict (specific
        entries for each general amino acid).
        :return: pdbAtomName
        """
        residueIncomplete, mAN = self.getMissingAtoms()
        residueName = self[0].getResidueName()
        if residueIncomplete:
            makeHFixfor = [pair[0] for pair in self.resiAtomTupleDict[residueName] if pair[1] in mAN and not pair[0] in mAN]
            residueSeqNum = self[0].getResiSeqNum()
            residueChainID = self[0].getChainID()
            return makeHFixfor
        else:
            return []

    def getResiName(self):
        return self[0].getResidueName()

    def getResiSeqNum(self):
        return self[0].getResiSeqNum()

    def hasInsertionCode(self):
        """
        Checks if the residue has an insertion code an therefore could have the highest residue number because of the
        applied offset.
        :return:
        """
        insertionCode = self[0].getInsertionCode()
        if not insertionCode:
            return False
        else:
            return True

    def isAminoAcid(self):
        """
        checks whether the residue name is a natural amino acid by checking the residue name against the list naturalAA.
        :return: boolean
        """
        residueName = self[0].getResidueName()
        isAA = False
        if residueName in Residue.naturalAA:
            isAA = True
        return isAA

    def isRightAtom(self):
        pass

    def isComplete(self):
        """
        should check whether the residue is complete (has the right number of atoms for the amino acid. With the dict
        'naturalAAAtomNumbers'.
        :return:
        """
        # return len(self) == Residue.naturalAAAtomNumbers[self[0].getResidueName()]
        residueName = self[0].getResidueName()
        missingResis = Residue.residueAtomDict[residueName]

    # def append(self, atom):
    #     super(Residue, self).append(atom)
    #     atom.registerResidue(self)
    #
    # def getMyNumber(self):
    #     return self.chain.getMyNumber(self.id)

    def append(self, atom, force=False):
        if not atom.getPDBAtomName() in self.getAtomNames():
            super(Residue, self).append(atom)
            return True
        if force:
            super(Residue, self).append(atom)
            return True
        return False


class NoResidueError(Exception):
    pass


def main():
    global start_time
    print 'INFO: +++ starting PDB2INS. +++ '
    start_time = time.time()
    main.__doc__ = head

    print main.__doc__
    #if version[:3] < '2.7':
    #    print 'Please consider upgrading your python. This program requires python 2.7'
    #    exit()
    parser = CommandlineParser()
    global options
    if not options:
        options = parser()
        # print options
    Data()


if __name__ == '__main__':
    main()


