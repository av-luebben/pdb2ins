__author__ = 'anna'
"""
second project pdb2hkl
by Anna Vera Luebben
start July 2015

read structure factor pdb file and generate .hkl file for SHELXL.
"""

import os
from sys import exit

buildin_raw_Input = raw_input


def raw_input(*args, **kwargs):
    inputString = buildin_raw_Input(*args, **kwargs)
    if inputString.lower() == 'q' or inputString.lower() == 'exit':
        print '*** pdb2hkl has been terminated ***'
        exit()
    return inputString


class CommandLineParser(object):

    def __init__(self):
        self.options = {'filename': None,
                        'o': None,
                        'i': False}
        self.validates = {'filename': self._validateDummy,
                          'o': self._validateDummy,
                          'i': self._validateDummy}

    def __call__(self, *args, **kwargs):
        from sys import argv
        try:
            self.options['filename'] = argv[1]
        except IndexError:
            self.options['filename'] = None
            return self.options
        activeOption = False
        for arg in argv[2:]:
            # print 'iterating over argv: ', arg
            if not activeOption:
                # print 'i am in not active options'
                activeOption = self._slash(arg)
                # print 'this is the active option', activeOption
                if not activeOption:
                    # print 'im in not active options.'
                    self.options[arg.lstrip('-')] = True
                else:
                    pass
            else:
                # print self.validates[activeOption](arg)
                if self.validates[activeOption](arg):
                    self.options[activeOption] = arg
                    activeOption = False
                else:
                    print 'ERROR: Wrong type of argument for option {}'.format(activeOption)
                    exit(1)
        return self.options

    def _slash(self, arg):
        arg = arg.lstrip('-')
        if arg not in self.options:
            print 'ERROR: Unknown cmd line option: {}'.format(arg)
            exit(2)
        else:
            if not self.options[arg]:
                return False
            else:
                return arg

    def _validateDummy(self, arg):
        #print 'validating dummy for ', arg
        return True


# class Pdb2insInterface(object):
#
#     # def __init__(self):
#     #     self.filename = None
#
# def pdb2insInterface(filename):


class IO(object):

    def __init__(self):
        self.workfile = None
        self.filename = None
        self.outputFilename = None
        self.datafile = None

    def askFilename(self):
        """

        :return:
        """
        self.filename = options['filename']
        if str(self.filename).endswith('.pdb'):
            sfFileNameTry = ''.join(str(self.filename).split('.')[:-1]) + '-sf.cif'
            if not os.path.isfile(sfFileNameTry):
                self.filename = None
            else:
                self.filename = sfFileNameTry
                print 'INFO: Using file {}'.format(sfFileNameTry)
        try:
            if not os.path.isfile(self.filename) and '@' not in self.filename:
                newstring = str(self.filename[:4].upper())+(self.filename[-4:])
                if not os.path.isfile(newstring) and not os.path.isfile(self.filename.lower()):
                    print 'INFO: Given filename not valid.'
                    self.filename = None
                if os.path.isfile(newstring):
                    self.filename = newstring
                try:
                    os.path.isfile(self.filename.lower())
                    self.filename = self.filename.lower()
                except:
                    print 'ERROR in pdb2hkl with: ', self.filename
                    pass
                    # print self.filename
                # if os.path.isfile(self.filename.lower()):
                #     self.filename = self.filename.lower()
        except TypeError:
            self.filename = None
        if not self.filename:
            # the next to lines are only there for Automated call of pdb2ins with a filename, otherwise not used.
            # print '** ERROR Filename not correct. **'
            # exit()
            while True:
                self.filename = raw_input("\nEnter name of a structure factor pdb file to read. "
                                          "To download a structure factor pdb file enter \'@<PDBCODE>\': ")
                if not os.path.isfile(self.filename) and not '@' in self.filename:
                    newstring = str(self.filename[:4].upper())+(self.filename[-4:])
                    if not os.path.isfile(newstring) and not os.path.isfile(self.filename.lower()):
                        print 'INFO: File \'{}\' not found.'.format(self.filename)
                    if os.path.isfile(newstring):
                        self.filename = newstring
                        break
                    if os.path.isfile(self.filename.lower()):
                        self.filename = self.filename.lower()
                        break
                else:
                    break
        else:
            self.filename = self.filename
        if self.filename.startswith('@'):
            if not len(self.filename[1:]) == 4:
                print 'ERROR: Given PDB code not valid.'
                exit(2)
            else:
                self.filename = self.fetchPDB(self.filename[1:])
        else:
            pass

    def fetchPDB(self, pdbCode, force=False):
        """
        http://www.rcsb.org/pdb/files/r3FBXsf.ent.gz
        Given a PDB code this function retrieves the PDB file from the PDB database directly.
        :param pdbCode: A valid pdb code
        :param force: If force is true an existing file with the 'pdbfile' name is overwritten by the new one.
        :return:
        """
        import urllib
        import gzip
        import os
        import string

        if options['o']:
            pdbFile = ''.join(str(options['o']).split('.')[:-1]) + '-sf.cif'
        else:
            pdbFile = pdbCode.lower() + '-sf.cif'
        remoteCode = string.upper(pdbCode)
        # if not os.path.exists(pdb_dir):
        #     os.mkdir(pdb_dir)
        if not os.path.exists(pdbFile) or force:
            try:
                filename = urllib.urlretrieve(
                    'http://www.rcsb.org/pdb/files/r' +
                    remoteCode + 'sf.ent.gz')[0]
            except:
                print "WARNING: {} not found.\n".format(pdbCode)
            else:
                if os.path.getsize(filename) > 0:  # If 0, then pdb code was invalid
                    try:
                        open(pdbFile, 'w').write(gzip.open(filename).read())
                        print "INFO: Fetched structure factor file for PDB code: {}".format(pdbCode)
                    except IOError:
                        print 'IO ERROR. No file found. \nNo structure factor file available for this PDB code.'
                        os.remove(pdbFile)
                else:
                    print "WARNING: {} not valid.\n".format(pdbCode)
                os.remove(filename)
        return pdbFile

    def read(self):
        """

        :return:
        """
        self.askFilename()
        # if not options['filename']:
        #     self.askFilename()
        # else:
        #     self.filename = options['filename']
        try:
            self.workfile = open(self.filename, 'r')
        except IOError:
            self.askFilename()
        else:
            self.datafile = self.workfile.readlines()
            # for line in self.datafile:
            #     print 'x', line

    def askOutputFilename(self):
        """

        :return:
        """
        self.filename = options['filename']
        if self.filename.startswith('@'):
            self.filename = self.filename[1:]
        defaultName = os.path.splitext(self.filename)[0] + '.hkl'
        if not options['i']:
            self.outputFilename = raw_input("\nEnter name of .hkl file to be created [{}]: ".format(defaultName))
            if not self.outputFilename:
                self.outputFilename = defaultName
            elif not self.outputFilename.endswith('.hkl'):
                self.outputFilename += '.hkl'
        else:
            self.outputFilename = options['o']
            if not self.outputFilename:
                self.outputFilename = defaultName
            elif not self.outputFilename.endswith('.hkl'):
                self.outputFilename += '.hkl'

    def writeFile(self, data):
        """
        writes the joined strings from DATA in an output file either named by the user in askOutputFilename or calls it
        newfile.ins.
        :param data: Gives joined string of all data needed for the .ins file
        :return: A new file in .hkl
        """
        self.askOutputFilename()
        if not self.outputFilename:
            print "INFO: FILE newfile.hkl was written."
            f = open("newfile.hkl", 'w')
            f.write(''.join(data))
        else:
            print "INFO: File {} was written.".format(self.outputFilename)
            f = open("{}".format(self.outputFilename), 'w')
            f.write(''.join(data))
        f.close()


class Data(object):

    def __init__(self):
        self.io = IO()
        self.io.read()
        self.loop = False
        self.underscore = False
        self.hklf = None
        self.parameterDict = {}
        self.dataDict = {}
        self.loopcounterBefore = None
        self.loopDataDict = {}
        self.foundUnmerged = False
        self.parameter = []
        self.parameterCount = {}
        self.dataString = []
        self.outfile = []
        self.biggestPlus = 0
        self.biggestMinus = 0
        self.biggestMeas = 0
        self.rescale = False
        self.factor = None
        self.errormessage = False
        self.readContent()
        self.checkRescaleData()
        self.writeString()
        self.io.writeFile(self.dataString)
        # self.makeHKLFglobal()

    # def makeHKLFglobal(self):
    #     global hklf
    #     hklf = self.hklf

    def getHKLF(self):
        return self.hklf

    def readContent(self):
        """
        Takes all lines from the original file. If the line is empty, nothing is done with it. A starting loop is
        detected when the line starts with 'loop_' which has to be followed by the list of parameters this loop
        contains.
        The parameters must start with '_' and the first line without a parameter at the beginning will be considered
        as indicator the at the loop can end (from here on only the values for the aforementioned parameters are given).
        Loops are counted. All parameters of loops are saved into parameterDict with the loop number it occurred in.
         All values/data lines from a loop (not parameters) are given to readLoopDataLine function.
        :return:
        """
        parametercounter = 0
        loopcounter = 0
        if not self.io.datafile:
            exit()
        for i, line in enumerate(self.io.datafile):
            if len(line.lstrip()) <= 1:
                continue
            elif line[0] == "#":
                loopEnds = True
                parametercounter = 0
                continue
            elif line[:5] == 'loop_':
                # print 'here the loops starts.', i, line
                loopcounter += 1
                self.loop = True
                loopEnds = False
            elif line.lstrip()[0] == '_':
                # print '_line found', i, line
                if self.loop and loopEnds:
                    # print 'here the loop ends.'
                    self.loop = False
                    self.underscore = True
                if self.loop:
                    # print 'loop', i, line, \
                    #       'loopcounter', loopcounter
                    self.parameterDict[line.split(' ')[0].rstrip('\n')] = (loopcounter, parametercounter)
                    # print self.parameterDict
                    parametercounter += 1
                    # print 'parametercounter', parametercounter
                    self.underscore = True
                else:
                    # print 'else line', i, line
                    pass
            else:
                # if i <= 80:
                    # print 'this arrives in else', i, line

                if self.loop:
                    # if not line[:5] == 'loop_':
                    #     print 'this is the loopcounter: {} and this is loop {},  ' \
                    #           'loopEnds is {} and parametercounter {}'.format(loopcounter, self.loop, loopEnds,
                    #                                                           parametercounter)

                    # if i <= 80:
                    #     print 'this arrives in else, within a loop', i, line
                    if not loopEnds:
                        self.parameterCount[loopcounter] = parametercounter
                        # print self.parameterCount
                        self.underscore = False
                        self.readLoopDataLine(line, loopcounter)
                # else:
                #     loopEnds = True
        # for k, v in self.parameterDict.items():
        #     print 'parameterDict', k, v

    def readLoopDataLine(self, line, loopcounter):
        """
        All lines with data within a loop arrive here. The lines are saved to loopDataDict with the loop number they
        come from as key.
        :param line:
        :param loopcounter:
        :return:
        """
        # print type(line), line
        # if '?' not in line:
        # if loopcounter == 3:
        #     print 'this is loop data line', line, 'this is the loopcounter', loopcounter
        #     exit()
        try:
            self.loopDataDict[loopcounter] += line[:-1].split()
        except KeyError:
            self.loopDataDict[loopcounter] = line[:-1].split()
        # self.loopDataDict[self.loopcounterBefore] = looplist
        # self.loopcounterBefore = loopcounter

    def findHKL(self):
        if self.isH()[0] and self.isK()[0] and self.isL()[0]:
            hLoop, hNumber = self.isH()[1:]
            kLoop, kNumber = self.isK()[1:]
            lLoop, lNumber = self.isL()[1:]
            # print hLoop, hNumber, kLoop, kNumber, lLoop, lNumber
            return hLoop, hNumber, kLoop, kNumber, lLoop, lNumber
        else:
            print 'ERROR: hkl data not found. Please check file.'
            exit()

    def isH(self):
        # print self.parameterDict.keys()
        # print 'isH', self.parameterDict['_refln.index_h'], self.parameterDict['_refln.index_h'][0]
        try:
            hLoop = self.parameterDict['_refln.index_h'][0]
            hParam = self.parameterDict['_refln.index_h'][1]
            # print 'this is were h is in parameterDict', hLoop, hParam
            return True, hLoop, hParam
        except KeyError:
            return False, None, None

    def isK(self):
        try:
            kLoop = self.parameterDict['_refln.index_k'][0]
            kParam = self.parameterDict['_refln.index_k'][1]
            # print 'this is were k is in parameterDict', kLoop, kParam, self.parameterDict
            return True, kLoop, kParam
        except KeyError:
            return False, None

    def isL(self):
        try:
            lLoop = self.parameterDict['_refln.index_l'][0]
            lParam = self.parameterDict['_refln.index_l'][1]
            # print 'this is were l is in parameterDict', lLoop, lParam
            return True, lLoop, lParam
        except KeyError:
            return False, None

    def findRightParameter(self):
        if self.isIminus()[0] and self.isIplus()[0]:
            print 'INFO: Found unmerged I data. Please use HKLF 4 in SHELXL.'
            self.hklf = 4
            IminusLoop, IminusParam = self.isIminus()[1]
            IminusSigmaLoop, IminusSigmaParam = self.isIminus()[2]
            IplusLoop, IplusParam = self.isIplus()[1]
            IplusSigmaLoop, IplusSigmaParam = self.isIplus()[2]
            if IminusLoop == IplusLoop:
                self.foundUnmerged = True
                return IminusLoop, IminusParam, IminusSigmaLoop, IminusSigmaParam, IplusLoop, IplusParam, IplusSigmaLoop, IplusSigmaParam
            else:
                print'ERROR: Found unmerged data in different loops. Please handle.'
                self.hklf = None
                exit(1)
        elif self.isImeas()[0]:
            print 'INFO: Found I data. Please use HKLF 4 in SHELXL.'
            self.hklf = 4
            ImeasLoop, ImeasParam = self.isImeas()[1]
            ImeasSigmaLoop, ImeasSigmaParam = self.isImeas()[2]
            return ImeasLoop, ImeasParam, ImeasSigmaLoop, ImeasSigmaParam
        elif self.isFminus()[0] and self.isFplus()[0]:
            print 'INFO: Found unmerged F data. Please use HKLF 3 in SHELXL.'
            self.hklf = 3
            FminusLoop = self.isFminus()[1][0]
            FminusParam = self.isFminus()[1][1]
            FminusSigmaLoop = self.isFminus()[2][0]
            FminusSigmaParam = self.isFminus()[2][1]
            FplusLoop = self.isFplus()[1][0]
            FplusParam = self.isFplus()[1][1]
            FplusSigmaLoop = self.isFplus()[2][0]
            FplusSigmaParam = self.isFplus()[2][1]
            if FminusLoop == FplusLoop:
                self.foundUnmerged = True
                return FminusLoop, FminusParam, FminusSigmaLoop, FminusSigmaParam, FplusLoop, FplusParam, FplusSigmaLoop, FplusSigmaParam
            else:
                print'ERROR: Found unmerged data in different loops. Please handle.'
                self.hklf = None
                exit(2)
        elif self.isFmeas()[0]:
            print 'INFO: Found F data. Please use HKLF 3 in SHELXL.'
            self.hklf = 3
            Floop = self.isFmeas()[1][0]
            FmeasParam = self.isFmeas()[1][1]
            FSigmaLoop = self.isFmeas()[2][0]
            FSigmaParam = self.isFmeas()[2][1]
            return Floop, FmeasParam, FSigmaLoop, FSigmaParam
        else:
            print 'ERROR: No complete set of structure factor data found.'
            self.hklf = None
            exit(3)

    def isIplus(self):
        try:
            iplusLoop = self.parameterDict['_refln.pdbx_I_plus']
            iplusSigmaLoop = self.parameterDict['_refln.pdbx_I_plus_sigma']
            return True, iplusLoop, iplusSigmaLoop
        except KeyError:
            # if not self.errormessage:
            #     print 'ERROR: Data is incomplete. Sigma for I plus might be missing.'
            #     self.errormessage = True
            return False, None

    def isIminus(self):
        try:
            iminusLoop = self.parameterDict['_refln.pdbx_I_minus']
            iminusSigmaLoop = self.parameterDict['_refln.pdbx_I_minus_sigma']
            return True, iminusLoop, iminusSigmaLoop
        except KeyError:
            # if not self.errormessage:
            #     print 'ERROR: Data is incomplete. Sigma for I minus might be missing.'
            #     self.errormessage = True
            return False, None

    def isImeas(self):
        # keywords = ['_refln.I_meas_au', '_refln.I_meas', '_refln.F_squared_meas', '_refln.intensity_meas']
        # keywords2 = ['_refln.I_meas_sigma_au', '_refln.I_meas_sigma', '_refln.I_meas_sigma', '_refln.intensity_sigma']
        # for i in keywords:
        #     try:
        #         imeasLoop = self.parameterDict[i]
        #     except KeyError:
        #         continue
        try:
            imeasLoop = self.parameterDict['_refln.I_meas_au']
            imeasSigmaLoop = self.parameterDict['_refln.I_meas_sigma_au']
            return True, imeasLoop, imeasSigmaLoop
        except KeyError:
            try:
                imeasLoop = self.parameterDict['_refln.I_meas']
                imeasSigmaLoop = self.parameterDict['_refln.I_meas_sigma']
                return True, imeasLoop, imeasSigmaLoop
            except KeyError:
                try:
                    imeasLoop = self.parameterDict['_refln.F_squared_meas']
                    imeasSigmaLoop = self.parameterDict['_refln.F_squared_sigma']
                    return True, imeasLoop, imeasSigmaLoop
                except KeyError:
                    try:
                        imeasLoop = self.parameterDict['_refln.intensity_meas']
                        imeasSigmaLoop = self.parameterDict['_refln.intensity_sigma']
                        return True, imeasLoop, imeasSigmaLoop
                    except KeyError:
                        # if not self.errormessage:
                        #     print 'ERROR: Data is incomplete. Sigma for I meas might be missing.'
                        #     self.errormessage = True
                        return False, None

    def isFplus(self):
        try:
            fplusLoop = self.parameterDict['_refln.pdbx_F_plus']
            fplusSigmaLoop = self.parameterDict['_refln.pdbx_F_plus_sigma']
            return True, fplusLoop, fplusSigmaLoop
        except KeyError:
            # if not self.errormessage:
            #     print 'ERROR: Data is incomplete. Sigma for F plus might be missing.'
            #     self.errormessage = True
            return False, None

    def isFminus(self):
        try:
            fminusLoop = self.parameterDict['_refln.pdbx_F_minus']
            fminusSigmaLoop = self.parameterDict['_refln.pdbx_F_minus_sigma']
            return True, fminusLoop, fminusSigmaLoop
        except KeyError:
            # if not self.errormessage:
            #     print 'ERROR: Data is incomplete. Sigma for F minus might be missing.'
            #     self.errormessage = True
            return False, None

    def isFmeas(self):
        try:
            fmeasLoop = self.parameterDict['_refln.F_meas_au']
            fmeasSigmaLoop = self.parameterDict['_refln.F_meas_sigma_au']
            # print 'this is were Fmeas is in the parameterDict', fmeasLoop
            return True, fmeasLoop, fmeasSigmaLoop
        except KeyError:
            try:
                fmeasLoop = self.parameterDict['_refln.F_meas']
                fmeasSigmaLoop = self.parameterDict['_refln.F_meas_sigma']
                return True, fmeasLoop, fmeasSigmaLoop
            except KeyError:
                # if not self.errormessage:
                #     print 'ERROR: Data is incomplete. Sigma for F meas might be missing.'
                #     self.errormessage = True
                return False, None

    def findFlack(self):
        try:
            flackLoop = self.parameterDict['_refln.status']
            # print flackLoop[0], flackLoop[1]
            return flackLoop[0], flackLoop[1]
        except KeyError:
            print 'INFO: No free R flack found.'
            return None

    def checkRescaleData(self):
        self.parameter = self.findRightParameter()
        # hasFlack = self.findFlack()
        # if hasFlack:
        #     self.flackLoop = hasFlack[0]
        #     self.flackPosition = hasFlack[1]

        # print ' rescale', self.rescale
        if self.foundUnmerged:
            minusLoop = self.parameter[0]
            # minusNum = self.parameter[1]
            minusSigmaLoop = self.parameter[2]
            # minusSigmaNum = self.parameter[3]
            plusLoop = self.parameter[4]
            # plusNum = self.parameter[5]
            plusSigmaLoop = self.parameter[6]
            # plusSigmaNum = self.parameter[7]
            enumerator = self.parameterCount[minusLoop]
            dataline = []
            for word in self.loopDataDict[minusLoop]:
                dataline.append(word)
                if len(dataline) == enumerator:
                    self.control(dataline)
                    dataline = []
        if not self.foundUnmerged:
            measLoop = self.parameter[0]
            # measNum = self.parameter[1]
            measSigmaLoop = self.parameter[2]
            # measSigmaNum = self.parameter[3]
            # flackLoop, flackNum = self.findFlack()
            enumerator = self.parameterCount[measLoop]
            dataline = []
            for word in self.loopDataDict[measLoop]:
                dataline.append(word)
                if len(dataline) == enumerator:
                    self.control(dataline)
                    dataline = []
        # print 'rescale reloaded', self.rescale
        if self.rescale:
            self.rescaleData()

    def rescaleData(self):
        if self.biggestMeas > 9999999:
            # print 'big meas'
            self.factor = 9999999/self.biggestMeas
        elif self.biggestPlus > 9999999 and self.biggestPlus >= self.biggestMinus:
            # print 'big plus'
            self.factor = 9999999/self.biggestPlus
        elif self.biggestMinus > 9999999 and self.biggestMinus > self.biggestPlus:
            # print 'big minus'
            self.factor = 9999999/self.biggestMinus
        # print self.factor

    def getRescaleFactor(self):
        return float(self.factor)

    def control(self, dataline):
        """
        Takes the dataline and extracts h, k, l values. depending on data, also merged or unmerged data is extracted.
        If a '?' is found, the line is not used for the .hkl file. Is minus ans plus present in the line and only one
        of those or their errror has a '?', the other one is used.
        If one of the values is greater than 9 999 999, an
        error message is printed. Here the data should be rescaled in the future.
        """
        notMinus = False
        notPlus = False
        if self.foundUnmerged:
            h = dataline[self.findHKL()[1]]
            k = dataline[self.findHKL()[3]]
            l = dataline[self.findHKL()[5]]
            minus = dataline[self.parameter[1]]
            minussigma = dataline[self.parameter[3]]
            plus = dataline[self.parameter[5]]
            plussigma = dataline[self.parameter[7]]
            if '?' in h or '?' in k or '?' in l:
                return
            if '?' in minus or '?' in minussigma:
                notMinus = True
            if '?' in plus and '?' in plussigma:
                notPlus = True
            if notPlus and notMinus:
                return
            if not notMinus:
                minus = float(minus)
                # if minus > 9999999.:
                    # print 'minus', minus
                if minus > 9999999. and self.biggestMinus <= minus:
                    self.biggestMinus = minus
                    # print self.biggestMinus
                    self.rescale = True
            if not notPlus:
                plus = float(plus)
                # if plus > 9999999.:
                    # print 'plus', plus
                if plus > 9999999. and self.biggestPlus <= plus:
                    self.biggestPlus = plus
                    # print self.biggestPlus
                    self.rescale = True
        else:
            noFlack = False
            h = dataline[self.findHKL()[1]]
            k = dataline[self.findHKL()[3]]
            l = dataline[self.findHKL()[5]]
            # print 'this is position ', self.findHKL()[3], 'and this is data at that position in dataline', k
            # print 'this is the dataline', dataline
            # exit()
            meas = dataline[self.parameter[1]]
            meassigma = dataline[self.parameter[3]]
            if '?' in h or '?' in k or '?' in l or '?' in meas or '?' in meassigma:
                return
            # print h, k, l, meas, meassigma
            try:
                meas = float(meas)
            except ValueError:
                print '\nERROR: Syntax error, the file contains numerical parameters mixed with letters: ', meas
                print 'Attention: The program will not write an .hkl file from insufficient data.\n'
                exit(1)
            # if meas > 9999999.:
            #     print 'meas', meas
            if meas > 9999999. and self.biggestMeas <= meas:
                self.biggestMeas = meas
                # print self.biggestMeas
                self.rescale = True
        # print 'rescale 2.0', self.rescale
            # except ValueError:
            #     print 'File format could not be read.'
            #     exit()

    def writeString(self):
        """
        the function findhkl is called to check whether hkl data is available. If not, the program is terminated.
        Next the 'find rigth parameter' function is called to find the appropriate data for the .hkl file. This function
        distinguishes between merged and unmerged data. The loop numbers for the data is extracted and saved to
        variables.
        The loopdataDict is called with the correct loop number to extract the data for the chosen parameter.
        The loopdata is reassembled from chunks(words) into lines which are in the length of the number of parameters
        given in the loop (enumerator). Each complete line of data is given to the flush function.
        """
        # print 'started write string.'
        hLoop, hnum, kLoop, kNum, lLoop, lNum = self.findHKL()
        # print 'i am in write string and have hkl:', self.findHKL()
        if not hLoop == kLoop == lLoop:
            print 'ERROR: Not all data necessary in one loop.'
            exit(2)
        # self.parameter = self.findRightParameter()
        if self.foundUnmerged:
            # print 'found unmerged'
            minusLoop = self.parameter[0]
            # minusNum = self.parameter[1]
            minusSigmaLoop = self.parameter[2]
            # minusSigmaNum = self.parameter[3]
            plusLoop = self.parameter[4]
            # plusNum = self.parameter[5]
            plusSigmaLoop = self.parameter[6]
            # plusSigmaNum = self.parameter[7]
            enumerator = self.parameterCount[minusLoop]
            dataline = []
            if not minusLoop == minusSigmaLoop == plusLoop == plusSigmaLoop:
                print 'ERROR: Not all data necessary in one loop.'
                exit(2)
            for word in self.loopDataDict[minusLoop]:
                dataline.append(word)
                if len(dataline) == enumerator:
                    self.flush(dataline)
                    dataline = []
        else:
            # print 'found merged'
            measLoop = self.parameter[0]
            # measNum = self.parameter[1]
            measSigmaLoop = self.parameter[2]
            # measSigmaNum = self.parameter[3]
            # flackLoop, flackNum = self.findFlack()
            enumerator = self.parameterCount[measLoop]
            # print 'measloop', measLoop, measSigmaLoop, enumerator
            # print self.parameterCount
            # for key, value in self.parameterCount.items():
            #     print 'parameterCount', key, value
            if not measLoop == measSigmaLoop:
                print 'ERROR: Not all data necessary found in one loop.'
                exit(2)
            dataline = []
            for word in self.loopDataDict[measLoop]:
                dataline.append(word)
                # print dataline[-1], len(dataline), enumerator
                if len(dataline) == enumerator:
                    # print 'flush dataline', enumerator
                    self.flush(dataline)
                    dataline = []
            if self.foundUnmerged:
                self.dataString.append('   0   0   0       0       0')
            if not self.foundUnmerged:
                self.dataString.append('   0   0   0       0       0       0       0')

    def flush(self, dataline):
        """
        Takes the dataline and extracts h, k, l values. depending on data, also merged or unmerged data is extracted.
        If a '?' is found, the line is not used for the .hkl file. If one of the values is greater than 9 999 999, an
        error message is printed. Here the data should be rescaled in the future.
        """
        notMinus = False
        notPlus = False
        plusString = None
        minusString = None
        noFlack = False
        # print 'I am in flush', dataline
        if self.foundUnmerged:
            # print dataline
            h = dataline[self.findHKL()[1]]
            k = dataline[self.findHKL()[3]]
            l = dataline[self.findHKL()[5]]
            minus = dataline[self.parameter[1]]
            minussigma = dataline[self.parameter[3]]
            plus = dataline[self.parameter[5]]
            plussigma = dataline[self.parameter[7]]
            try:
                f = dataline[self.findFlack()[1]]
                if f == 'f':
                    flag = -1
                else:
                    flag = 1
            except IndexError:
                noFlack = True
            except TypeError:
                noFlack = True
            # First, all values with an '?' where a number should be have to be rejected.
            if '?' in h or '?' in k or '?' in l:
                return
            try:
                h = int(h)
                k = int(k)
                l = int(l)
            except ValueError:
                print '\nERROR: Syntax error in structure factor file. Please handle.\n' \
                      'Attention: The program will not write an .hkl file from insufficient data.\n'
                exit()
            # Second, all minus and plus of one line have to be checked for '?' also. Part of the line can be kept,
            # if only one of them (plus or minus) has a '?'. If both have, the hole line is rejected.
            if '?' in minus or '?' in minussigma:
                notMinus = True
            if '?' in plus and '?' in plussigma:
                notPlus = True
            if notPlus and notMinus:
                return
            # Only those without a '?' can be converted to float, all others raise a ValueError.
            try:
                if not notMinus:
                    minus = float(minus)
                    minussigma = float(minussigma)
                if not notPlus:
                    plus = float(plus)
                    plussigma = float(plussigma)
            except ValueError:
                print '\nERROR: Syntax error in structure factor file. Please handle.\n' \
                      'Attention: The program will not write an .hkl file from insufficient data.\n'
                exit()
            # Third, if munis or plus has the value of zero, the part should be rejected, too.
            if plus == 0:
                notPlus = True
            if minus == 0:
                notMinus = True
            if notPlus and notMinus:
                return
            if self.rescale:
                # minus = float(minus)
                # plus = float(plus)
                factor = self.getRescaleFactor()
                if not notMinus:
                    minus = factor * minus
                    minussigma = factor * minussigma
                if not notPlus:
                    plus = factor * plus
                    plussigma = factor * plussigma
            # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {} {:>7.6n} {:>3} \n'.format(h, k, l, myformat(minus),
            #                                                                                 minussigma, int(-1)))
            # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {} {:>7.6n} {:>3} \n'.format(h, k, l, myformat(plus),
            #                                                                                 plussigma, int(-1)))
            # if len(minus) >= 7:
            # print type(h), type(k), type(l), type(minus), type(minussigma), type(int(-1))
            # exit()
            # if minus > 999999 and plus > 999999:
            #     self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.7n} {:>7.6n} {:>3} \n'.format(h, k, l, minus, minussigma,
            #                                                                              int(-1)))
            #     self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.7n} {:>7.6n} {:>3} \n'.format(h, k, l, plus, plussigma,
            #                                                                              int(1)))
            if not notMinus and not notPlus:
                if minus > 999999 and plus > 999999:
                    minusString = '{: >7.7n}'.format(minus)
                    plusString = '{: >7.7n}'.format(plus)
                elif minus > 999999:
                    minusString = '{: >7.7n}'.format(minus)
                    plusString = '{: >7.6n}'.format(plus)
                elif plus > 999999:
                    plusString = '{: >7.7n}'.format(plus)
                    minusString = '{: >7.6n}'.format(minus)
                elif int(minus) < 0 and int(plus) < 0:
                    minusString = '{: >-6.4n}'.format(minus)
                    plusString = '{: >-6.4n}'.format(plus)
                elif int(minus) < 0:
                    minusString = '{: >-6.4n}'.format(minus)
                    plusString = '{: >7.6n}'.format(plus)
                elif int(plus) < 0:
                    plusString = '{: >-6.4n}'.format(plus)
                    minusString = '{: >7.6n}'.format(minus)
                else:
                    plusString = '{: >7.5n}'.format(plus)
                    minusString = '{: >7.5n}'.format(minus)
            if not notMinus and notPlus:
                plusString = None
                if minus > 999999:
                    minusString = '{: >7.7n}'.format(minus)
                elif int(minus) < 0:
                    minusString = '{: >-6.4n}'.format(minus)
                else:
                    minusString = '{: >7.5n}'.format(minus)
            if notMinus and not notPlus:
                minusString = None
                if plus > 999999:
                    plusString = '{: >7.7n}'.format(plus)
                elif int(plus) < 0:
                    plusString = '{: >-6.4n}'.format(plus)
                else:
                    plusString = '{: >7.5n}'.format(plus)
            if len(str(plussigma).strip()) > 8:
                plussigma = float(str(plussigma)[:8])
            if len(str(minussigma).strip()) > 8:
                minussigma = float(str(minussigma)[:8])
                # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.7n} {:>7.6n} {:>3} \n'.format(h, k, l, minus, minussigma,
                #                                                                          int(-1)))
                # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.6n} {:>7.6n} {:>3} \n'.format(h, k, l, plus, plussigma,
                #                                                                          int(1)))
            # elif plus > 999999:
                # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.6n} {:>7.6n} {:>3} \n'.format(h, k, l, minus, minussigma,
                #                                                                          int(-1)))
                # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.7n} {:>7.6n} {:>3} \n'.format(h, k, l, plus, plussigma,
                #                                                                          int(1)))
            # elif int(minus) < 0 and int(plus) < 0:
            #     self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >-6.5n} {:>7.6n} {:>3} \n'.format(h, k, l, minus, minussigma,
            #                                                                              int(-1)))
            #     self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >-6.5n} {:>7.6n} {:>3} \n'.format(h, k, l, plus, plussigma,
            #                                                                              int(1)))

                # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >-6.5n} {:>7.6n} {:>3} \n'.format(h, k, l, minus, minussigma,
                #                                                                          int(-1)))
                # self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.6n} {:>7.6n} {:>3} \n'.format(h, k, l, plus, plussigma,
                #                                                                          int(1)))
            # elif int(plus) < 0:
            #     self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >7.6n} {:>7.6n} {:>3} \n'.format(h, k, l, minus, minussigma,
            #                                                                              int(-1)))
            #     self.dataString.append('{:>4.0f} {:>3.0f} {:>3.0f} {: >-6.5n} {:>7.6n} {:>3} \n'.format(h, k, l, plus, plussigma,
            #                                                                              int(1)))

            # print plusString, minusString
            # if len(plusString) > 7:
            #     plusString = plusString[:7]
            # if len(minusString) > 7:
            #     minusString = minusString[:7]
            hminus = int(h) * (-1)
            kminus = int(k) * (-1)
            lminus = int(l) * (-1)
            # print dataline
            if noFlack:
                flag = 1
            # print flack
            if plusString:
                if len(plusString) > 8:
                    plusString = plusString[:8]
                self.dataString.append('{:>4.0f}{:>4.0f}{:>4.0f}{: >8s}{:> 8n}{:>4} \n'.format(h, k, l, plusString,
                                                                                                   plussigma, flag))
            if minusString:
                if len(minusString) > 8:
                    minusString = minusString[:8]
                self.dataString.append('{:>4.0f}{:>4.0f}{:>4.0f}{: >8s}{:> 8n}{:>4} \n'.format(hminus, kminus,
                                                                                                   lminus, minusString,
                                                                                                   minussigma, flag))
            # print '{: >7s}'.format(plusString), '{: >7s}'.format(minusString)
        else:
            noFlack = False
            h = dataline[self.findHKL()[1]]
            k = dataline[self.findHKL()[3]]
            l = dataline[self.findHKL()[5]]
            # print 'this is position ', self.findHKL()[3], 'and this is data at that position in dataline', k
            # print 'this is the dataline', dataline
            # exit()
            meas = dataline[self.parameter[1]]
            meassigma = dataline[self.parameter[3]]
            if '?' in h or '?' in k or '?' in l or '?' in meas or '?' in meassigma:
                return
            try:
                h = int(h)
                k = int(k)
                l = int(l)
                meas = float(meas)
                meassigma = float(meassigma)
            except ValueError:
                print '\nERROR: Syntax error in structure factor file. Please handle.\n' \
                      'Attention: The program will not write an .hkl file from insufficient data.\n'
                quit()
            try:
                f = dataline[self.findFlack()[1]]
                if f == 'f':
                    flag = -1
                else:
                    flag = 1
            except IndexError:
                noFlack = True
            except TypeError:
                noFlack = True
            if meas == 0:
                return
            if self.rescale:
                meassigma = float(meassigma)
                factor = self.getRescaleFactor()
                meas = factor * meas
                meassigma = factor * meassigma
            # try:
            #     ttt = float(meas)
            #     if ttt > 9999999:
            #         print 'This files data has to be scaled. The format of the hkl file will be wrong.'
            # except ValueError:
            #     print 'File format could not be read.'
            #     exit()

            # if meas > 999999:
            #     measString = '{: >7.7n}'.format(meas)
            # elif meas > 0:
            #     measString = '{: >-6.4n}'.format(meas)
            # else:
            #     measString = '{: >7.6n}'.format(meas)
            # if len(measString) > 7:
            #     measString = measString[:7]
            if 'e' not in str(meas):
                if len(str(meas).strip()) >= 8:
                    if str(meas).startswith('0'):
                        meas = float(str(meas)[:7])
                    else:
                        meas = float(str(meas)[:8])
                if len(str(meassigma).strip()) >= 8:
                    if str(meassigma).startswith('0'):
                        meassigma = float(str(meassigma)[:7])
                    else:
                        meassigma = float(str(meassigma)[:8])
            else:  # due to scientific notation, small numbers get an exponential number starting with e-5.
                x = '{:f}'.format(meas)
                meas = float(x[:7])
            # print meas, type(meas), meassigma, type(meassigma)

            # print h, k, l, meas
            # print '{:> 8.5n} {:> 8.5g} {:> 8n}'.format(meas, meas, meas)
            # try:
            #     flag = dataline[self.findFlack()[2]]
            # except IndexError:
            #     noFlack = True
            # except TypeError:
            #     noFlack = True
            if noFlack:
                self.dataString.append('{:>4}{:>4}{:>4}{:> 8n}{:> 8n} \n'.format(h, k, l, meas, meassigma))
                # if meas > 999999:
                # elif int(meas) < 0:
                #     self.dataString.append('{:>4} {:>3} {:>3} {:>-6.5n} {:>7.6n} \n'.format(h, k, l, meas, meassigma))
                # else:
                #     self.dataString.append('{:>4} {:>3} {:>3} {:>7.6n} {:>7.6n} \n'.format(h, k, l, meas, meassigma))
                # # self.dataString.append('{:>4} {:>3} {:>3} {:>7f} {:>7f} \n'.format(h, k, l, meas, meassigma))
            if not noFlack:
                # if meas > 999999:
                # print meas, type(meas), meassigma, type(meassigma)
                self.dataString.append('{:>4}{:>4}{:>4}{:> 8n}{:> 8n}{:> 4} \n'.format(h, k, l, meas, meassigma,
                                                                                           flag))
                # elif int(meas) < 0:
                #     self.dataString.append('{:>4} {:>3} {:>3} {:>-6.5n} {:>7.6n} {:>3} \n'.format(h, k, l, meas,
                #                                                                                   meassigma, flack))
                # else:
                #     self.dataString.append('{:>4} {:>3} {:>3} {:>7.6n} {:>7.6n} {:>3}\n'.format(h, k, l, meas,
                #                                                                                 meassigma, flack))
                # self.dataString.append('{:>4} {:>3} {:>3} {:>7f} {:>7f} {:>3} \n'.format(h, k, l, meas, meassigma,
                #                                                                             flack))


        # self.outFile.write('{} {} {} {} {}'.format())

# def myformat(x, l=8):
#     """
#     takes a number and formats it into fixed format standard for hkl files, so the whole number, including sign (only
#     given when number negative) takes up 7 digits.
#     :param x: float
#     :param l: jength of total string
#     :return:
#     """
#     t = '{:>-8.7f}'.format(x)[:l]
#     i1 = -2
#     i2 = -3
#     e = ''
#     try:
#         l = int(t[-2])
#     except ValueError:
#         l = int(t[-1])
#         i1 = -3
#         i2 = -4
#         e = '.'
#     p = 1 * l >= 5
#     try:
#         return t[:i1] + str(int(t[i1]) + p) + e
#     except ValueError:
#         return t[:i2] + str(int(t[i2]) + p)+'.'


def run(forceOptions=None):

    """
    ########################################################################
    #                                 PDB2HKL                              #
    #                     by Anna V. Luebben (Version 2017/1)              #
    ########################################################################

    Reads a pdb (.cif) file and generates a .hkl file for SHELXL.

    Usage:
        pdb2hkl <filename|@pdbcode> [options]

    <filename|@pdbcode> Exact name of a file in PDB format or a legal pdb code
    prefixed by '@'.
    """
    global options
    if not forceOptions:
        print run.__doc__
        parser = CommandLineParser()
        options = parser()
    else:
        options = forceOptions
    hklf = Data().getHKLF()
    if forceOptions:
        print 'INFO: Closing pdb2hkl.'
        return hklf


if __name__ == '__main__':
    run()