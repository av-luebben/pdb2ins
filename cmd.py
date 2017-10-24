__author__ = 'anna'

"""
This is an addition to pdb2ins.py, the pdb2ins program.
"""
from sys import exit


class CommandlineParser(object):

    def __init__(self):
        """
        ######################################################################################
        #                                    PDB2INS                                         #
        #                        by Anna V. Luebben (Version 2016/3)                         #
        ######################################################################################

        Usage:
            pdb2ins <filename|@pdbcode> [options]

        <filename|@pdbcode> Exact name of a file in PDB format or legal pdb code prefixed by '@'.

        Options:
            The following options are available with '-' as a prefix:
            'w' followed by the wavelength in Angstroms to enter a wavelength manually.
            'h' followed by a value to enter hkl format in SHELX syntax.
            'c' followed by the six coordinates (a,b,c,alpha,beta,gamma) in Angstroms and Degree
                to enter the cell. (No spaces allowed.)
            's' followed by the space group. (No spaces allowed.)
            'i' skip user input interrupts.
            'a' use anisotropic displacement data if available.
            'b' create a .hkl file from a PDB structure factor file (.cif) or pdb code.
            'd' followed by a filename to specify structure factor file in cif format.
            'e' all hydrogen atoms are transferred from the PDB file to the INS file.
            'o' followed by a filename to specify the output filename.
            'r' if a pdb code is given with '@', this option fetches the PDB file from REDO server.
            'z' followed by a number specifying the z value (number of molecules per cell).

        Example:
            pdb2ins myPDBFile.pdb -w 1.54178 -h 4 -c 100,20.5,30,90,90,90 -i -a -o myFile.ins
        """
        self.options = {'filename': None,
                        'w': None,  # enter wavelength
                        'h': None,  # give hklf value
                        's': None,  # enter spacegroup (no space)
                        'r': False,  # use pdb redo
                        'c': None,  # cell, 6 parameters (comma separated, no space)
                        'i': False,  # interactive modus (questions)
                        'a': False,  # use anis data?
                        'z': None,  # z-value, number of units per cell
                        'b': False,  # create hkl file
                        'd': None,  # give structure factor file name, if not pdb code given as filename
                        'e': False,  # keep H atoms
                        'o': None,  # specify output filename
                        'GUI': False}
        self.validates = {'filename': self._validate_dummy,
                          'w': self._validate_w,
                          'h': self._validate_hklf,
                          'c': self._validate_c,
                          's': self._validate_s,
                          'i': self._validate_dummy,
                          'a': self._validate_dummy,
                          'o': self._validate_dummy,
                          'z': self._validate_dummy,
                          'b': self._validate_dummy,
                          'd': self._validate_d,
                          'e': self._validate_dummy,
                          'r': self._validate_dummy,
                          'GUI': self._validate_dummy}

    def __call__(self, override=None):
        """
        Takes the arguments given after the program is called and interprets them.
        :param override: list of strings to be used instead of sys.argv.
        :return:
        """
        if not override:
            from sys import argv
        else:
            argv = override
            # self.options['GUI'] = True
        if '--help' in argv:
            # print __init__.__doc__
            print CommandlineParser.__init__.__doc__
            exit()
        try:
            self.options['filename'] = argv[1]
        except IndexError:
            self.options['filename'] = None
            return self.options
        activeOption = False
        for arg in argv[2:]:
            # print 'reading arg', arg
            if not activeOption:
                activeOption = self._slash(arg)
                if not activeOption:
                    self.options[arg.lstrip('-')] = True
                else:
                    pass
                    # print 'found option', arg, ', waiting for next argument to assign value.'
            else:
                # print 'assigning value', arg, 'to option', activeOption
                if self.validates[activeOption](arg):
                    self.options[activeOption] = arg
                    activeOption = False
                else:
                    print 'ERROR: Wrong type of argument for option: {}'.format(activeOption)
                    exit(1)
        return self.options

    def _slash(self, arg):
        arg = arg.lstrip('-')
        if arg not in self.options:
            print 'ERROR: Unknown cmd line option: {}'.format(arg)
            exit(2)
        else:
            if self.options[arg] == False:
                return False
            else:
                return arg

    def _validate_w(self, arg):
        """
        Controls the input given after 'w'. The wavelength entered must be between 0.2 and 5.0 Angstroms. Otherwise the
        wavelength will be dismissed as not sensible.
        :param arg: string, wavelength
        :return: boolean
        """
        # print arg, type(arg), float(arg), type(float(arg))
        if 0.2 <= float(arg) <= 5.0:
            return True
        else:
            return False

    def _validate_c(self, arg):
        """
        Controls the parameters given after 'c' (cell).
        There must be six arguments given for the command to be accepted. Also the cell should have an a greater than
        20 Angstroms and alpha should be between 20 and 160 degrees.
        :param arg:
        :return: boolean
        """
        # '20,20,20,90,90,90'
        try:
            arg = arg.split(',')
            arg = [float(i) for i in arg if i]
            cell_a = arg[0]
            alpha = arg[3]
            if cell_a <= 2.00 or not 20 < alpha < 160:
                return False
        except:
            return False
        if len(arg) == 6:
            return True

    def _validate_d(self, arg):
        """
        Contrals parameter given in cmd after 'd' (structure factor file)
        :param arg:
        :return:
        """
        if self.options['filename']:
            filename = str(self.options['filename'])
            if filename.startswith('@'):
                arg = filename
                return False
            else:
                return True

    def _validate_s(self, arg):
        """
        Controls parameter given in cmd after 's' (spacegroup).
        :param arg:
        :return: boolean
        """
        from spagsydata import testSpaceGroup
        doNotReplace1 = ['P 1', 'A 1', 'B 1', 'C 1', 'I 1', 'F 1']
        try:
            if arg.strip() not in doNotReplace1:
                arg = arg.replace(' 1 ', '').replace(' 1', '').lstrip()
            # if arg[1] == "R":
            #     x = self.cell[6] - self.cell[5]
            #     if x >= 20:
            #         self.spaceGroup[1] = "H"
            arg = arg.replace(' ', '')
        except:
            return False
        if testSpaceGroup(arg):
            return True
        else:
            return False

    def _validate_hklf(self, arg):
        """
        Controls if the option entered under h (hklf code) cmd is valid.
        There must be given 3 for F or 4 for F-squared data.
        :param arg:
        :return:
        """
        if arg in ('3', '4'):
            return True
        else:
            return False

    def _validate_dummy(self, arg):
        """
        Changes the Boolian from Flase to True for the option in question.
        :param arg:
        :return:
        """
        return True


if __name__ == '__main__':
    cP = CommandlineParser()
    options = cP()
    cell = options['c']

