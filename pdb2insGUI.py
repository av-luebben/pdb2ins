__author__ = 'Anna'
"""
Created in April and May 2016 by Anna Luebben as GUI for the program PDB2INS.
This program creates a GUI where a file can be loaded, information displayed and added, and a .ins and .hkl created
from a pdb file or PDB code using the program PDB2INS and its subroutines.
"""

from Tkinter import *
import os
# import pexpect
from tkFileDialog import askopenfilename
from subprocess import Popen, PIPE
from time import sleep
from getInfoForGui import Info
from ScrolledText import ScrolledText
from threading import Thread
from spagsydata import testSpaceGroup

import pdb2ins
import time
import sys
from cmd import CommandlineParser

stdout_bkp = sys.stdout

ALIVE = True


def make_executable(path):
    """
    makes file 'path' executable.
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)


def selectAll(event):
    event.widget.selection_range(0, END)


def openPDBFile():
    newFileName = askopenfilename()
    fileName.set(newFileName)


def openSFFile():
    newFileName = askopenfilename()
    sfFilename.set(newFileName)


def SlaveCallback(*args):
    """
    Callback function to be given to the pdb2ins module to receive output from 'print' statements instead of sys.stdout.
    :param args: list of strings.
    :return: None
    """
    for arg in args:
        arg = arg[0]
        if 'INFO:' in arg or 'ERROR:' in arg:
            t.insert(END, arg+'\n')


class Runner(Thread):
    def __init__(self):
        root.after(500, update)
        super(Runner, self).__init__()
        self.start()

    def run(self):
        """
        This function is called by pressing the 'Create INS file" button in the GUI.
        First the functions '...Missing' is called to make sure a value for all important fields are given.
        (pdb2ins must have: HKLF, if not using pdb2hkl subroutine, and the following if not in file: wavelength, cell,
        space group and zValue.)
        Afterwards the cmd string is completed
        """
        if not self.validateEverything():
            return
        tc = takeCell()
        cmdString = 'None {file} -w {waveLength} {cell} -i {zValue}' \
                    ' {HKLF} {anis} {redo} {pdb2hkl}' \
                    ' {sfFile} {output}'.format(file=fileName.get(),
                                                waveLength=waveLength.get(),
                                                cell='-c ' + cell.get().replace(' ', '') if tc else '',
                                                zValue='-z ' + zValue.get() if not zValue.get() else '',
                                                HKLF='-h ' + hKLF.get() if not pdb2hkl.get() else '',
                                                anis='-a' if anis.get() else '',
                                                redo='-r' if redo.get() else '',
                                                pdb2hkl='-b' if pdb2hkl.get() else '',

                                                sfFile='-d ' + sfFilename.get() if sfFilename.get() else '',
                                                output='-o ' + output.get() if output.get() else '')
        opt = CommandlineParser()
        opt = opt(cmdString.split())
        pdb2ins.setSlaveMode(opt, SlaveCallback)
        pdb2ins.main()

    # def run2(self):
    #     """
    #     This function is called by pressing the 'Create INS file" button in the GUI.
    #     First the functions '...Missing' is called to make sure a value for all important fields are given.
    #     (pdb2ins must have: HKLF, if not using pdb2hkl subroutine, and the following if not in file: wavelength, cell,
    #     space group and zValue.)
    #     Afterwards the cmd string is completed
    #     With this cmd string the program pdb2ins is called.
    #     While pdb2ins (or the subroutines) are running, the output is monitored for keywords (index).
    #     The keywords trigger either a command line response, break or an output to the textfield 't' in the GUI.
    #     """
    #     # print self.validateEverything()
    #     # exit()
    #     if not self.validateEverything():
    #         return
    #     # if not pdb2hkl.get() == 1:
    #     #     variable = isHklfMissing()
    #     #     if not variable:
    #     #         return
    #     # setSFFilename()
    #     # if pdb2hkl.get() == 1 and sfFilename.get() == '':
    #     #     if not fileName.get().startswith('@'):
    #     #         sfFilenameMissing()
    #     #         return
    #     # if zValue.get() == '':
    #     #     zValueMissing()
    #     #     return
    #     # w = waveLength.get()
    #     # if w == 'None' or w == '':
    #     #     wavelengthMissing()
    #     #     return
    #     # if not validateSpaceGroup(spaceGroup.get(), cell.get().replace(' ', '')):
    #     #     spaceGroup.set('')
    #     #     spaceGroupMissing()
    #     #     return
    #     # if not cell.get().count(',') == 5:
    #     #     validateCell()
    #     #     return
    #     tc = takeCell()
    #     try:
    #         cmdString = '-w {waveLength} {cell} -i {zValue}' \
    #                     ' {HKLF} {anis} {redo} {pdb2hkl}' \
    #                     ' {sfFile} {output}'.format(waveLength=waveLength.get(),
    #                                             cell='-c ' + cell.get().replace(' ', '') if tc else '',
    #                                             zValue='-z ' + zValue.get() if not zValue.get() else '',
    #                                             HKLF='-h ' + hKLF.get() if not pdb2hkl.get() else '',
    #                                             anis='-a' if anis.get() else '',
    #                                             redo='-r' if redo.get() else '',
    #                                             pdb2hkl='-b' if pdb2hkl.get() else '',
    #
    #                                             sfFile='-d ' + sfFilename.get() if sfFilename.get() else '',
    #                                             output='-o ' + output.get() if output.get() else '')
    #     except ValueError:
    #         return
    #     # print(cmdString)
    #     #cmdString = ['C:\Python27\python.exe' 'F:\Daten\pdb2ins\PDB2INS\pdb2ins.py ']
    #     # p = Popen(cmdString, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    #
    #     exeDir = os.path.dirname(sys.executable)
    #     # cmdString = ['{}/pdb2ins'.format(exeDir), fileName.get()] + [_ for _ in cmdString.split(' ') if _]
    #     try:
    #         make_executable(os.path.join(sys._MEIPASS, 'pdb2ins'))
    #     except AttributeError:
    #         cmdString = ['{}/pdb2ins'.format(exeDir), fileName.get()] + [_ for _ in cmdString.split(' ') if _]
    #     else:
    #         cmdString = [os.path.join(sys._MEIPASS, 'pdb2ins'),
    #                      fileName.get()] + [_ for _ in cmdString.split(' ') if _]
    #     c = ' '.join(cmdString)
    #     try:
    #         child = pexpect.spawn(c)
    #     except pexpect.ExceptionPexpect:
    #         global t
    #         t.delete(1.0, END)
    #         t.insert(END, '##############################\n'
    #                       'ERROR: PDB2INS is not installed\n'
    #                       '##############################\n'
    #                       'Please make sure that the pdb2ins executable is in the same directory as the PDB2INSGUI'
    #                       'executable.\n')
    #         return
    #         # cmdString = ['python2', '/home/anna/pdbtools/PDB2INS/pdb2ins.py',
    #                        fileName.get()] + [_ for _ in cmdString.split(' ') if _]
    #         # c = ' '.join(cmdString)
    #         # child = pexpect.spawn(c)
    #     while True:
    #         index = child.expect(['HKLF', 'water', ' filename ', 'created.', 'Hydrogen', 'X-RAY', 'anisotropic',
    #                               pexpect.EOF, 'INFO:', 'ERROR:'])
    #         # print child
    #         # cmdOutput = StringVar
    #         if index == 0:
    #             child.sendline('\n')
    #         elif index == 1:
    #             child.sendline('\n')
    #         elif index == 2:
    #             child.sendline('\n')
    #         elif index == 3:
    #             pass
    #         elif index == 4:
    #             child.sendline('\n')
    #         elif index == 5:
    #             child.sendline('\n')
    #             break
    #         elif index == 6:
    #             child.sendline('\n')
    #         elif index == 7:
    #             global ALIVE
    #             root.after(500, update)
    #             ALIVE = False
    #             # cmdOutput.append('PDB2INS finished running.')
    #             break
    #         elif index == 8:
    #             global cmdOutput
    #             # cmdOutput.append('INFO' + child.readline())
    #             cmdOutput.append(child.readline())
    #             # global t
    #             # t.insert(END, cmdOutput)
    #             # print '####', cmdOutput
    #             # raw_input()
    #             # for i in child.read():
    #             #     if 'INFO' in str(i):
    #             #         cmdOutput += str(i)
    #
    #             # print sys.stdout
    #             # print child.before
    #             # cmdOutput += str(sys.stdout)
    #         elif index == 9:
    #             # print 'ERROR found.'
    #             # global t
    #             t.insert(END, 'ERROR' + child.readline())
    #             root.after(50, update)
    #             child.sendline('q\n')

    def validateEverything(self):
        """
        Is called at the beginning of the run function and calls all validate or Missing function if variables are not
        as expected.
        return: Boolean
        """
        resetEntryColors()
        global t
        t.delete(1.0, END)
        fine = True
        if not pdb2hkl.get() == 1:
            fine = not all((isHklfMissing(), fine))
        setSFFilename()
        if pdb2hkl.get() == 1 and sfFilename.get() == '':
            if not fileName.get().startswith('@'):
                sfFilenameMissing()
                fine = False
        if zValue.get() == '':
            zValueMissing()
            fine = False
        w = waveLength.get()
        if w == 'None' or w == '':
            wavelengthMissing()
            fine = False
        c = cell.get().replace(' ', '')
        if not c.count(',') == 5:
            validateCell()
            fine = False
        s = spaceGroup.get()
        if s == '':
            spaceGroupMissing()
            fine = False
        if not fine:
            root.after(100, self.validateEverything)
        else:
            if not validateSpaceGroup(s, c):
                # spaceGroup.set('')
                spaceGroupWrong()
                root.after(500, self.validateEverything)
                fine = False
        return fine

        # else:
        #     return True


def resetEntryColors():
    """
    Resets all entry colors that can be validated to white when called.
    :return:
    """
    sfEntry['bg'] = 'white'
    zEntry['bg'] = 'white'
    hKLFEntry['bg'] = 'white'
    waveLengthEntry['bg'] = 'white'
    sGEntry['bg'] = 'white'
    cellEntry['bg'] = 'white'


def sfFilenameMissing():
    """
    called if the sf file entry is missing. Changes backhround colour of structure file entry field.
    Also a text message is printed into the text box.
    :return:
    """
    # sfEntry['bg'] = 'white'
    # global t
    # t.delete(1.0, END)
    if not sfFilename.get():
        t.insert(END, 'ATTENTION! Please enter a structure factor file!')
        sfEntry['bg'] = 'red2'
        # root.after(500, sfFilenameMissing)
        # if not sfFilename.get() == '':
    # t.delete(1.0, END)


def zValueMissing():
    """
    If the zValue is missing, this function will be called through validateEverything(). The
    :return:
    """
    # zEntry['bg'] = 'white'
    # global t
    # t.delete(1.0, END)
    if not zValue.get():
        zEntry['bg'] = 'red2'
        t.insert(END, 'ATTENTION! Z Value is missing. Please enter a value for z.\n')
        # root.after(500, zValueMissing)


def isHklfMissing():
    """
    Called at the beginning of run() through validateEverything()
    While the Entry for HKLF is missing the function is in a loop and the Entry field will be red.
    Is the HKLF entry given the loop is broken and run() will be callable.
    The loop is not entered when a hkl file will be created by pdb2hkl, because the subroutine can set hklf variable!
    """
    # defaultbg = root.cget('bg') bg is 'light grey' under linux, which is not the color the entry field had originally!
    # hKLFEntry['bg'] = 'white'
    if not pdb2hkl.get() == 1:
        if not hKLF.get() == '':
            return False
        if hKLF.get() == '':
            # global t
            # t.delete(1.0, END)
            t.insert(END, 'ATTENTION! Value for HKLF missing. Please enter a value for HKLF.\n')
            # root.after(500, update)
            hKLFEntry['bg'] = 'red2'
            # root.after(100, isHklfMissing)
            return True
        # if not hKLFEntry.get():
        #     # root.after(500, update)
        #     t.insert(END, 'ATTENTION! Value for HKLF missing. Please enter a value for HKLF.\n')
        #     hKLFEntry['bg'] = 'red2'
            # root.after(100, isHklfMissing)
            # if not hKLFEntry.get() == '':
            #     hklfboolean.set('True')
            #     hKLFEntry['bg'] = 'white'
            #     t.delete(1.0, END)
            #     return True
            # return True
        return False
    else:
        return False
        # while True:
    #     sleep(1)
    #     quesiton = p.communicate()
    #     for _ in quesiton:
    #         print _
    #     if 'HKLF code' in quesiton[0]:
    #         p.communicate(input='4\n')


def wavelengthMissing():
    """
    Function is called if the wavelength is missing. Changes background of wave length entry field and writes a text
    message to the text box.
    :return:
    """
    # waveLength.set('')
    # waveLengthEntry['bg'] = 'white'
    # global t
    # t.delete(1.0, END)
    if not waveLength.get():
        waveLengthEntry['bg'] = 'red2'
        t.insert(END, 'ATTENTION! Wavelength is missing. Please enter a wavelength in Angstrom.\n')
        # root.after(500, wavelengthMissing)
        # if not waveLength.get() == '':
        #     hklfboolean.set('True')
        #     waveLengthEntry['bg'] = 'white'
    # t.delete(1.0, END)


def spaceGroupMissing():
    """
    Changes the background of the space group entry field when the space group is missing.
    Writes a text message into the text box.
    :return:
    """
    # sGEntry['bg'] = 'white'
    # global t
    # t.delete(1.0, END)
    if not spaceGroup.get():
        sGEntry['bg'] = 'red2'
        t.insert(END, 'ATTENTION! Space group is missing. Please enter a space group in Hermann-Mauguin '
                      'nomenclature.\n')
        # root.after(500, spaceGroupMissing)
        # if not spaceGroup.get() == '':
        #     hklfboolean.set('True')
    # t.delete(1.0, END)


def spaceGroupWrong():
    """
    Changes the background of the space group entry field to red, if the space group is missing. Also the space group is
    validated by the function validateSpaceGroup. The validate function needs the cell and the space group.
    :return:
    """
    # sGEntry['bg'] = 'white'
    # global t
    # t.delete(1.0, END)
    if not validateSpaceGroup(spaceGroup.get(), cell.get()):
        # root.after(500, update)
        sGEntry['bg'] = 'red2'
        t.insert(END, 'ATTENTION! Spacegroup is not correct. Please enter a correct space group in Hermann-Mauguin '
                      'nomenclature.\n')
        # root.after(100, spaceGroupWrong)
        # if validateSpaceGroup(spaceGroup.get(), cell.get()):
        #     hklfboolean.set('True')
        #     sGEntry['bg'] = 'white'
    # t.delete(1.0, END)


def validateSpaceGroup(sg, ce):
    """
    The space group extracted from the PDB file or entered by the user is tested if it is a legal space group.
    For space group test the cell is necessary. space group is abbreviated as in pdb2ins and given to the
    testSpaceGroup function from the pdb2ins program subroutine spagsydata. spagsydata return a boolean describing if
    the space group is a valid space group or not.
    :param sg (spacegroup), ce (cell)
    :return: Boolean
    """
    if not sg:
        return False
    try:
        arg = ce.split(',')
        arg = [float(i) for i in arg if i]
        cell_a = arg[0]
        cell_b = arg[1]
        cell_c = arg[2]
        alpha = arg[3]
        beta = arg[4]
        gamma = arg[5]
        cellNewFormat = [cell_a, cell_b, cell_c, alpha, beta, gamma]

        doNotReplace1 = ['P 1', 'A 1', 'B 1', 'C 1', 'I 1', 'F 1']
        if sg.strip() not in doNotReplace1:
            sg = sg.replace(' 1 ', '').lstrip()
        if sg[1] == "R":
            x = cellNewFormat[6] - cellNewFormat[5]
            if x >= 20:
                sg[1] = "H"
        shortenedSpaceGroup = sg.replace(' ', '')
        return testSpaceGroup(shortenedSpaceGroup)
    except IndexError:
        return False
    # if testSpaceGroup(shortenedSpaceGroup):
    #     return True
    # else:
    #     return False


def validateCell():
    """
    Called from the validate everything function.
    Function changes cell entry field to red while not six parameters are given for cell, separator ','.
    :return:
    """
    # cellEntry['bg'] = 'white'
    # global t
    # t.delete(1.0, END)

    if not cell.get().count(',') == 5:
        cellEntry['bg'] = 'red2'
        t.insert(END, 'ERROR: cell must have six arguments.\n '
                      'Please enter cell again in the format "a, b, c, alpha, beta, gamma".\n')
        # root.after(500, validateCell)
        # if cell.get().count(',') == 5:
        #     cellEntry['bg'] = 'white'
        #     root.after(500, update)
    # root.after(500, update)
    # t.delete(1.0, END)
    """
    Called from the run function.
    This functions makes sure the variable 'cell' has 6 entries (a,b,c,alpha,beta,gamma). If not ValueError is raised.
    :param cell: StringVar
    :return:cell or ValueError
    """
    # # defaultbg = root.cget('fg')
    # if not cell.count(',') == 5:
    #     cellEntry['bg'] = 'red2'
    #     CustomValueError()
    # else:
    #     cellEntry['bg'] = 'white'
    #     return cell


def takeCell():
    """
    This function checks whether the cell has ben changed. If the cell is the same as extracted from the pdb file
    return is False, else return True
    :return: Boolean
    """
    cellnew = cell.get().replace(' ', '')
    cellold = originalcell.get().replace(' ', '')
    if cellnew == cellold:
        return False
    else:
        return True


def CustomValueError():
    t.insert(END, 'ERROR: cell must have six arguments.\n '
                  'Please enter cell again in the format: a, b, c, alpha, beta, gamma\n')


# def buildcmdString():
#     options = {'filename': fileName.get()}
#     info = Info(options)
#     w = wavelength.get()
#     cell = validateCell(cell.get().replace(' ', ''))
#     if not w == info.wavelength:
#         cmdstring =+ '-w {waveLength} '.format(waveLength=w)
#     if not cell == info.cell:
#         cmdstring += '-c {cell} '.format(cell=cell)


def resetGui():
    """
    This function resets all GUI fields to the default when started. This way no artifacts from a previous run prevail.
    :return:
    """
    global t
    t.delete(1.0, END)
    redoCheck.deselect()
    cell.set('')
    spaceGroup.set('')
    zValue.set('')
    output.set('')
    hKLF.set('')
    anisCheck.deselect()
    sfFilename.set('@PDBCode / FileName')
    waveLength.set('')
    createHKLCheck.deselect()
    controlSFfield()


def isFilename():
    """
    this function checks if a filename or PDB code was entered into the filename entry.
    :return:
    """
    f = fileName.get()
    if f == '@PDBCode / FileName':
        t.delete(1.0, END)
        t.insert(END, 'ATTENTION! Filename or PDB code missing! '
                      'Please enter a filename or "@" followed by a four character PDB code.\n')
        return False
    if f.startswith('@'):
        if not len(f) == 5:
            t.delete(1.0, END)
            t.insert(END, 'ATTENTION! The PDB Code given was not valid. Please enter "@" followed by a four '
                          'character PDB code\n')
            return False
    return True


def loadFile():
    """
    This functions calls the python file getInfoForGui (Info) and gives 'options'.
    Options must contain a filename and, if the filename starts with an '@', also 'r'.
    Otherwise pdb2ins doesn't know which website to use for the pdb code.
    getInfoForGui will now return incomplete File if the file does not contain the remark200Line. could be problematic!
    """
    resetGui()
    if not isFilename():
        return
    options = {'filename': fileName.get(), 'r': checkredo(), 'o': None, 'i': True, 'd': None}
    info = Info(options)
    if not info.everythingOkay:
        t.delete(1.0, END)
        t.insert(END, '####### ERROR: File/PDB Code not found! #######')
        root.after(500, update)
        return
    if info.incompleteFile:
        t.insert(END, '####### ATTENTION: File might not contain X-ray diffraction data! #######')
        root.after(500, update)
        # return
    if not info.wavelength:
        if not info.wavelengthLine == 'None' and info.wavelengthLine:
            t.insert(END, 'ERROR: Could not find valid wavelength in pdb file. Please enter wavelength manually.\n'
                          'Information found in pdb file:\n' + info.wavelengthLine + '\n')
    else:
        waveLength.set(info.wavelength)
    cell.set(str(info.cell)[1:-1])
    originalcell.set(str(info.cell)[1:-1])
    spaceGroup.set(info.spaceGroup)
    zValue.set(info.zValue)
    output.set(setOutputFilename())
    # print ''
    root.after(500, update)


def update():
    """
    This functions automatically updates the Textfield for the command line output every 500 milliseconds, but only
    while pdb2ins is running.
    :return:
    """
    global t, cmdOutput
    for cmd in cmdOutput:
        t.insert(END, cmd)
    cmdOutput = []
    if ALIVE:
        root.after(500, update)


def checkredo():
    """
    Called by loadfile function.
    Simple True or False return if the redo command should be included in the cmd string.
    :return:Boolean
    """
    if redo.get() != 0:
        return True
    else:
        return False


def setOutputFilename():
    """
    Called by loadFile().
    The output filename is generated from the input filename. If the input filename is a pdb Code, the '@' is cut out
    and '.ins' inserted at the end. If the filename ends with a '.pdb' this is replaced by '.ins'.
    """
    f = fileName.get()
    if f[-4:] == '_a.pdb':
        return f.replace('_a.pdb', '.ins')
    if f[0] == '@':
        new = f.replace('@', '') + '.ins'
        return new
    else:
        new2 = "".join(f.split('.')[:-1]) + '.ins'
        return new2


def controlSFfield():
    """
    The PDB2HKL checkbox is connected to the sf input file Entry field. By default the Entry box is disabled.
    Should the pdb2hkl checkbox be activated, the entry box is only activated if the filename is not a pdb code,
    but a .pdb file.
    """
    # print 'this is pdb2hkl status', pdb2hkl.get()
    if pdb2hkl.get() == 0:
        sfEntry.configure(state='disable')
    else:
        global fileName, t
        if '@' in fileName.get():
            sfEntry.configure(state='disable')
            if '@PDBCode' not in fileName.get():
                t.insert(END, 'Structure factor file will be downloaded automatically via given pdb code.\n')
        else:
            sfEntry.configure(state='normal')


def setSFFilename():
    f = fileName.get()
    if sfFilename.get() == '@PDBCode / FileName':
        if not f.startswith('@'):
            sfFilename.set('')
            global t
            # t.delete(0, END)
            if pdb2hkl.get() == 1:
                t.insert(END, 'Attention! Please enter a structure factor file to continue. ')
    if f.startswith('@'):
        sfFilename.set('')
        # print 'sf filename should be reset to f: ', sfFilename.get()


class printInfoText(object):
    variableDict = {'Load': 'Either a file in PDB format or "@" followed by a legal 4 character PDB code can be '
                            'entered.\n'
                            'When a PDB code is entered the .pdb file will be downloaded automatically from the RCSB '
                            'Protein Data Bank\n'
                            'and saved in the form "PDB code" + "_a.pdb" to the same folder as the .ins file.\n'
                            'After selecting a file or entering a PDB code, press the "LOAD" button.\n'
                            'Important information will be displayed and missing information can be added manually.\n'
                            'If a PDB code is given, the checkbox "PDB redo" can be selected to use the PDB redo server'
                            ' files instead of the RCSB Protein Data Bank.\n'
                            'The pdb file given should conform to the Protein Data Bank notes "Atomic Coordinate and '
                            'Bibliographic Entry Format Description 3.30".\n',
                    'PDBredo': 'When the PDB redo box is checked the PDB file will be downloaded from the PDB_REDO '
                               'server.\n PDB2INS will download the BEST TLS version of the selected PDB code.\n'
                               'For more information visit the PDB_REDO website.\n',
                    'wavelength': 'The wavelength of the X-rays used to collect the data in Angstroms.\n',
                    'zValue': 'The Z value is in general the number of units per cell.\n'
                              'Here Z value can be understood as the number of polymeric chains in a unit cell.\n'
                              'In the case of heteropolymers, Z is the number of occurrences of the most populous '
                              'chain.\n',
                    'cell': 'Enter the unit cell dimensions in the form a, b, c, alpha, beta, gamma  in Angstroms and '
                            'Degrees.\n'
                            'Six parameters are expected.\n',
                    'hklf': 'THe HKLF value defines the format of the reflection data in the .hkl file.\n'
                            'Most common is 3 or 4 corresponding to a reflection file containing Amplitudes or '
                            'Intensities, respectively, but a number from 1 to 6 is legal.\n'
                            'Please consult the SHELX websites for more information on this SHELXL instruction.\n'
                            'Should the "Create HKL file" option be used, the HKLF value will be derived from the '
                            'structure factor file directly.\n',
                    'spacegroup': 'The space group is displayed in the format it was given in the .pdb file.\n'
                                  'If missing, enter the space group in Hermann-Mauguin notation, e.g. P 1 21 1.\n',
                    'outputfile': 'You can choose your output filename and enter it here. \n'
                                  'Should no filename be given, the input filename will be used with the '
                                  'ending ".ins".\n'
                                  'When a PDB code was given as source, the default output will be the PDB code + '
                                  '".ins."\n'
                                  'The output filename will also be given to all automatically downloaded files with '
                                  'the appropriate suffix.\n'
                                  'For example a PDB code "@1a4m" is given and the output filename "hydrolase1a4m.ins"'
                                  'is chosen,\n'
                                  'PDB2INS writes the files "hydrolase1a4m.ins" and "hydrolase1a4m_a.pdb".\n'
                                  'This will also apply to the structure factor file and .hkl file, if the option '
                                  '"Create HKL file" is active.\n',
                    'anis': 'When the PDB file contains anisotropic data, selecting the anisotropic check box will '
                            'include the anisotropic atom data. \n'
                            'The default is converting anisotropic data to isotropic. \n'
                            'Anisotropic refinement is possible afterwards in SHELXL with the "ANIS" instruction.\n',
                    'pdb2hkl': 'PDB2INS can create a .hkl file from a structure factor file as given in the RCSB '
                               'Protein Data Base.\n'
                               'Should a PDB code be given as input file, the structure factor file will be downloaded '
                               'automatically together with the .pdb file.\n'
                               'The output file name will be derived from the filename specified in the "Output file" '
                               'field.\n'
                               'A structure factor file downloaded from the RCSB Protein Data Bank will be saved with '
                               'the suffix "-sf.cif".\n'}

    def __init__(self, variable):
        self.text = printInfoText.variableDict[variable]

    def __call__(self, *args, **kwargs):
        global t
        t.delete(1.0, END)
        t.insert(END, '******************************\n')
        t.insert(END, self.text)
        t.insert(END, '******************************\n')

root = Tk()
cmdOutput = []
root.wm_title('PDB2INS')

"""
Setting the needed variables.
"""
fileName = StringVar()
fileName.set('@PDBCode / FileName')
sfFilename = StringVar()
sfFilename.set('@PDBCode / FileName')
# fileName.set('/home/anna/pdbtools/PDB2INS/testfiles/1a4m.pdb')
hklfboolean = BooleanVar()
hklfboolean.set('False')
child = StringVar()
waveLength = StringVar()
# waveLength.set()
zValue = StringVar()
cell = StringVar()
originalcell = StringVar()
# cell.set()
hKLF = StringVar()
# hKLF.set('4')
spaceGroup = StringVar()
output = StringVar()
anis = IntVar()
redo = IntVar()
pdb2hkl = IntVar()

"""
Here the GUI front is defined.
"""
Label(root, text='PDB File:', justify=LEFT).grid(row=0, column=0, pady=20, sticky=S)
fileEntry = Entry(root, textvariable=fileName, width=55)
fileEntry.bind("<FocusIn>", selectAll)
fileEntry.grid(row=0, column=1, columnspan=4, pady=20, sticky=W+S)
Button(root, text='Browse', justify=LEFT, command=openPDBFile).grid(row=0, column=4, pady=20, sticky=S)
Frame(root, width=20).grid(row=0, column=5, pady=20, sticky=S)
Button(root, text='Load', justify=LEFT, command=loadFile).grid(row=0, column=6, pady=20, sticky=E+W+S)
Frame(root, width=20).grid(row=0, column=7, pady=20, sticky=S)
Button(root, text='?', command=printInfoText('Load')).grid(row=0, column=8, pady=20, sticky=S)

Label(root, text='PDB redo').grid(row=1, column=0)
redoCheck = Checkbutton(root, var=redo)
redoCheck.grid(row=1, column=1, sticky=W)
Frame(root, width=20).grid(row=1, column=2)
Button(root, text='?', command=printInfoText('PDBredo')).grid(row=1, column=3)

Frame(root, height=15).grid(row=2)

Label(root, text='Wavelength').grid(row=3, column=0)
waveLengthEntry = Entry(root, width=50, textvariable=waveLength)
waveLengthEntry.grid(row=3, column=1)
Frame(root, width=20).grid(row=3, column=2)
Button(root, text='?', command=printInfoText('wavelength')).grid(row=3, column=3)
Frame(root, width=20).grid(row=3, column=4)
Label(root, text='Z Value').grid(row=3, column=5)
zEntry = Entry(root, width=50, textvariable=zValue)
zEntry.grid(row=3, column=6)
Frame(root, width=20).grid(row=3, column=7)
Button(root, text='?', command=printInfoText('zValue')).grid(row=3, column=8)


Label(root, text='Cell').grid(row=4, column=0)
cellEntry = Entry(root, width=50, textvariable=cell)
cellEntry.grid(row=4, column=1)
Frame(root, width=20).grid(row=4, column=2)
Button(root, text='?', command=printInfoText('cell')).grid(row=4, column=3)
Frame(root, width=20).grid(row=4, column=4)
Label(root, text='HKLF').grid(row=4, column=5)
hKLFEntry = Entry(root, width=50, textvariable=hKLF)
hKLFEntry.grid(row=4, column=6)
Frame(root, width=20).grid(row=4, column=7)
Button(root, text='?', command=printInfoText('hklf')).grid(row=4, column=8)


Label(root, text='Space group').grid(row=5, column=0)
sGEntry = Entry(root, width=50, textvariable=spaceGroup)
sGEntry.grid(row=5, column=1)
Frame(root, width=20).grid(row=5, column=2)
Button(root, text='?', command=printInfoText('spacegroup')).grid(row=5, column=3)
Frame(root, width=20).grid(row=5, column=4)
Label(root, text='Output File').grid(row=5, column=5)
outputEntry = Entry(root, width=50, textvariable=output)
outputEntry.grid(row=5, column=6)
Frame(root, width=20).grid(row=5, column=7)
Button(root, text='?', command=printInfoText('outputfile')).grid(row=5, column=8)


Label(root, text='Anisotropic').grid(row=6, column=0,)
anisCheck = Checkbutton(root, var=anis)
anisCheck.grid(row=6, column=1, sticky=W)
Frame(root, width=20).grid(row=6, column=2)
Button(root, text='?', command=printInfoText('anis')).grid(row=6, column=3)
Frame(root, width=20).grid(row=6, column=4)
# Label(root, text='PDB redo').grid(row=5, column=5)
# redoCheck = Checkbutton(root, var=redo)
# redoCheck.grid(row=5, column=6, sticky=W)
# Frame(root, width=20).grid(row=5, column=7)
# Button(root, text='?').grid(row=5, column=8)

Frame(root, height=10).grid(row=7)


sfEntry = Entry(root, textvariable=sfFilename, width=50)
sfEntry.bind("<FocusIn>", selectAll)
sfEntry.grid(row=8, column=6)
Label(root, text='Create HKL').grid(row=8, column=0,)
createHKLCheck = Checkbutton(root, command=controlSFfield, var=pdb2hkl)
createHKLCheck.grid(row=8, column=1, sticky=W)
Frame(root, width=20).grid(row=8, column=2)
Button(root, text='?', command=printInfoText('pdb2hkl')).grid(row=8, column=3)
Frame(root, width=20).grid(row=8, column=4)
Label(root, text='SF File').grid(row=8, column=5)

# sfFileEntry = Entry(root, textvariable=sfFilename, width=20)
# sfEntry.bind("<FocusIn>", selectAll)
# sfFileEntry.grid(row=8, column=6)  # , columnspan=4, pady=20, sticky=W+S
Frame(root, width=20).grid(row=8, column=7)
Button(root, text='Browse', command=openSFFile).grid(row=8, column=8)

Frame(root, height=20).grid(row=9)

Button(root, text='Create INS File', command=Runner).grid(row=10, columnspan=9, sticky=N+S+E+W)


t = ScrolledText(root)
# t = Text(root)
# t.config(state=DISABLED)
t.grid(row=11, column=0, rowspan=8, columnspan=9, sticky=N+S+E+W)

# t.insert(END, cmdOutput.get())
controlSFfield()
root.mainloop()
