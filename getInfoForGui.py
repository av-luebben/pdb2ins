__author__ = 'anna'
"""
This is an addition to pdb2ins.py
This .py will open a .pdb file and extract the parameters that can be added or changed while running pdb2ins.
"""
# import sys
# from os import listdir
# print(listdir(sys._MEIPASS+'/include'))
import os

from cmd import CommandlineParser
from spagsydata import getSymmCards, testSpaceGroup
from pdb2ins import IO
from spagsydata import testSpaceGroup


def test():
    # main.__doc__ = head
    #
    # print main.__doc__
    # if version[:3] < '2.7':
    #    print 'Please consider upgrading your python. This program requires python 2.7'
    #    exit()
    parser = CommandlineParser()
    global options
    options = parser()
    Info(options)


class Info(object):

    def __init__(self, options):
        self.options = options
        self.filename = None
        self.datafile = None
        self.anis = False
        self.neutron = False
        self.incompleteFile = False
        self.crystLine = None
        self.cell = None
        self.spaceGroup = None
        self.shortenedSpaceGroup = None
        self.zValue = None
        self.remark200Lines = []
        self.everythingOkay = True
        self.wavelengthLine = None
        self.wavelength = None
        self.checkFile()
        self.readFile()
        self.extractWavelength()

    def checkFile(self):
        """
        calls the function in pdb2ins that opens and reads the workfile after the function askfilename confirmed that
        the file exists or has downloaded it via fromDali.
        datafile = self.workfile.readlines()
        :return:
        """
        io = IO(self.options)
        try:
            io.read()
        except AttributeError:
            print 'ERROR: File could not be opened in Get Info For GUI.'
            self.everythingOkay = False
            return
        self.filename = io.workfile
        self.datafile = io.dataf

    def readFile(self):
        try:
            for line in self.datafile:
                if line[0] == "#":
                    continue
                if line[:6] == 'EXPDTA':
                    self.getCrystData(line)
                if line[:6] == "ATOM  " or line[:6] == "HETATM":
                    continue
                if line[:6] == "HET   ":
                    continue
                if line[:6] == "ANISOU":
                    self.anis = True
                else:
                    if line[:6] == "CRYST1":
                        self.crystLine = line
                        self.extractCell()
                        self.extractSpaceGroup()
                        self.extractZvalue()
                    if line[:10] == "REMARK 200":
                        self.remark200Lines.append(line)
                    if line[:6] == "SEQRES":
                        pass
                        # self.sequenceLines.append(line)
                    if line[:5] == "SCALE":
                        pass
                    else:
                        pass
        except TypeError:
            print 'ERROR: No file found in get Info for GUI.'
            self.everythingOkay = False
            return

    def getCrystData(self, line):
        if 'NEUTRON' in line:
            self.neutron = True
        if 'X-RAY DIFFRACTION' not in line and 'NEUTRON' not in line:
            self.incompleteFile = True

    # def getInfoFromLines(self, line):
    #     if line[:6] == "CRYST1":
    #         self.crystLine = line
    #         self.extractCell()
    #         self.extractSpaceGroup()
    #         self.extractZvalue()
    #     if line[:10] == "REMARK 200":
    #         self.remark200Lines.append(line)
    #     if line[:6] == "SEQRES":
    #         pass
    #         # self.sequenceLines.append(line)
    #     if line[:5] == "SCALE":
    #         pass
    #         # self.scaleLine.append(line)
    #     # header = Header()
    #     # header.interpretLine(line)

    def extractCell(self):
        """
        If the entry describes a structure determined by a technique other than X-ray crystallography,
         CRYST1 contains a = b = c = 1.0, alpha = beta = gamma = 90 degrees, space group = P 1, and Z = 1.
         """
        try:
            cell_a = float(self.crystLine[6:15])
            cell_b = float(self.crystLine[15:24])
            cell_c = float(self.crystLine[24:33])
            alpha = float(self.crystLine[33:40])
            beta = float(self.crystLine[40:47])
            gamma = float(self.crystLine[47:54])
            self.cell = [cell_a, cell_b, cell_c, alpha, beta, gamma]
        except ValueError:
            self.cell = None
        if self.cell:
            if (cell_a + cell_b + cell_c) < 1:
                self.cell = None
            if cell_a <= 20.00 or not 20 < alpha < 160:
                print "ERROR: Cell may not be correct! Please check."
                self.cell = None

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
        try:
            self.spaceGroup = self.crystLine[55:66].strip('\n')
            if self.spaceGroup.strip() not in doNotReplace1:
                self.spaceGroup = self.spaceGroup.replace(' 1 ', '').replace(' 1', '').lstrip()
            if self.spaceGroup[1] == "R":
                x = self.cell[6] - self.cell[5]
                if x >= 20:
                    self.spaceGroup[1] = "H"
            self.shortenedSpaceGroup = self.spaceGroup.replace(" ", "")
            if not testSpaceGroup(self.shortenedSpaceGroup):
                self.spaceGroup = None
        except ValueError:
            pass

    def extractZvalue(self):
        """The Z value is the number of polymeric chains in a unit cell. In the case of heteropolymers,
        Z is the number of occurrences of the most populous chain."""
        # validZ = ['0.5', '1', '2', '3', '4', '5', '6', '8', '9', '10', '11', '12', '16', '24', '36']
        # print 'this is zvalue',  self.crystLine[66:70].strip()
        try:
            self.zValue = int(self.crystLine[66:70])
        except KeyError:
            self.zValue = None
        except ValueError:
            self.zValue = None

    def extractWavelength(self):
        """
        searches the remark200 lines of the pdb file for the wavelength. If the wavelength is found,
        the user is asked if the wavelength is correct and has a chance to correct it.
        If a wavelength could not ge found, the user is asked to enter a wavelength.
        :return: wavelength in Angstrom.
        """
        # if options['w']:
        #     try:
        #         self.wavelength = float(options['w'])
        #     except TypeError:
        #         self.wavelength = None
        #     except ValueError:
        #         self.wavelength = None
        if not self.remark200Lines:
            print 'File might not contain X-ray data.'
            # self.incompleteFile = True
            self.wavelength = None
            # self.everythingOkay = False
            return
        for line in self.remark200Lines:
            if not self.wavelength and "WAVELENGTH" in line:
                try:
                    wavelengthExtracted = float(line.split(':')[-1].lstrip(' ').split(' ')[0].partition(';')[0])
                except ValueError:
                    try:
                        self.wavelengthLine += line
                    except TypeError:
                        self.wavelengthLine = line
                    continue
                if self.validateWavelength(wavelengthExtracted):
                    self.wavelength = wavelengthExtracted
                    return

    @staticmethod
    def validateWavelength(wavelength):
        if 0.2 <= float(wavelength) <= 5.0:
            return True
        else:
            return False


if __name__ == '__main__':
    test()
