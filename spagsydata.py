__author__ = 'anna'
"""
Written by Anna Vera Luebben in March 2015 for pdb2ins.py.
This is an addition to pdb2ins.py, the pdb2ins program.
It takes the space group of a data set an generates all symmetry operations from a
dictionary of symmetry generators.
The function getSymmCards takes the Space group in abbreviated shelx format and returns a string as SYMM line.
Version July 2016
"""

from numpy import matrix, array


def testSpaceGroup(spaceGroup):
    """
    Validates if the given space group (abbreviated) is correct by trying to find it in spagsy.
    :param spaceGroup: string, abbreviated space group
    :return: Boolean
    """
    try:
        generators = spagsy[spaceGroup]
        return True
    except KeyError:
        return False
    except ValueError:
        return False


def getSymmCards(spaceGroup):
    """
    Takes the shortened Space Group, extracts the Space groups' generators from spagsy and generates all symmetry
    operations possible for the space group out of the generators in the form of matrices and vectors.
    :param spaceGroup: shortened space group symbol
    :return: string, symmetry operations as matrices and vectors
    """
    generators = spagsy[spaceGroup]
    unity = generators[0]
    symmElements = [unity]
    for i, generator in enumerate(generators[1:]):
        symmElements.append(generator)
        symmopsOfGenerator = symmElements[1:]
        for element in symmopsOfGenerator:
            while True:
                newElement = element * generator  # Implement Generator multiplication
                if newElement == unity or newElement == generator or newElement == element or newElement in symmopsOfGenerator:
                    break
                # print newElement
                symmopsOfGenerator.append(newElement)
                element = newElement
        symmElements += symmopsOfGenerator

    for symmElement in symmElements:
        symmElement.normalize()
    # print [symmElement.asString() for symmElement in symmElements]
    symmElements = list(set(symmElements))
    return [symmElement.asString() for symmElement in symmElements]

# def getSymmCards(spaceGroup):
#     """
#     Takes the shortened Space Group, extracts the Space groups' generators from spagsy and generates all symmetry
#     operations possible for the space group out of the generators in the form of matrices and vectors.
#     :param spaceGroup: shortened space group symbol
#     :return: string, symmetry operations as matrices and vectors
#     """
#     generators = spagsy[spaceGroup]
#     unity = generators[0]
#     symmElements = [unity]
#     for i, generator in enumerate(generators[1:]):
#         i += 1
#         currentSecondGenerator = unity
#         j = 1
#         symmElements.append(generator)
#         # unityI = len(symmElements) - 1
#         do = True
#         while do:
#             newElement = symmElements[-1] * generator  # Implement Generator multiplication
#             if newElement == unity or newElement == generator or newElement == currentSecondGenerator:
#                 if j == i:
#                     j += 1
#                 try:
#                     symmElements.append(generators[j])
#                 except IndexError:
#                     break
#                 currentSecondGenerator = generators[j]
#                 j += 1
#             else:
#                 symmElements.append(newElement)
#     # print [symmElement.asString() for symmElement in symmElements]
#     symmElements = list(set(symmElements))
#
#     return [symmElement.asString() for symmElement in symmElements]


class Generator(object):
    def __init__(self, nmatrix, nvector):
        self.matrix = nmatrix
        self.vector = nvector

    def __str__(self):
        return '\n'.join((str(self.matrix), str(self.vector)))

    def __mul__(self, other):
        """
        Defines the multiplication of matrices as present in the symmetry operators, which are listed in spagsydata.
        Also the new vector is calculated here.
        :param other: one of the matrices or vectors needed for the calculation
        :return: new matrix and new vector created from the symmetry generators
        """
        newmatrix = self.matrix * other.matrix
        a11 = self.matrix[0, 0]
        a12 = self.matrix[0, 1]
        a13 = self.matrix[0, 2]
        a21 = self.matrix[1, 0]
        a22 = self.matrix[1, 1]
        a23 = self.matrix[1, 2]
        a31 = self.matrix[2, 0]
        a32 = self.matrix[2, 1]
        a33 = self.matrix[2, 2]
        t1 = other.vector[0]
        t2 = other.vector[1]
        t3 = other.vector[2]
        newvector = array([a11 * t1 + a12 * t2 + a13 * t3 + self.vector[0],
                           a21 * t1 + a22 * t2 + a23 * t3 + self.vector[1],
                           a31 * t1 + a32 * t2 + a33 * t3 + self.vector[2]])
        for i, x in enumerate(newvector):
            newvector[i] = (x + 99) % 1
            # if x < 0:
            #     newvector[i] = 1 + x
            # elif x == 1:
            #     newvector[i] = 0
            # elif x > 1:
            #     newvector[i] = x-1
            # else:
            #     pass
        # print newvector
        newGen = Generator(newmatrix, newvector)
        return newGen

    def __eq__(self, other):
        """
        tests the list of symmetry operation for duplicates: in combination with hash, all duplicates are removed.
        :param other: the  new matrix has to be compered with this matrix
        :return: all matrices that are unique.
        """
        # print 'Test for equality was preformed.'
        return (self.matrix == other.matrix).all()

    def __hash__(self):
        """
        checks all generated symmetry elements for doubles and returns a list with all unique symmetry operations
        :return:
        """
        strings = []
        for x in self.matrix:
            for y in x:
                strings.append(str(y))
        for x in self.vector:
            strings.append(str(x))
        return hash(''.join(strings))

    def normalize(self):
        for i, x in enumerate(self.vector):
            self.vector[i] = (x + 99) % 1

    def asString(self):
        """
        Transforms the Matrices of all possible symmetry operations into symmetry instruction for SYMM line in .ins.
        Each line of the 3x3 matrix is read and transformed to one port of the SYMM operation.
        Also the string for the .ins file is generated.
        :return: string, SYMM line
        """
        ref = ['x', 'y', 'z', '']
        allsymmlines = []
        # badSymm = False
        for k in xrange(3):
            vector = self.vector[k]
            # if vector < 0:
            #     badSymm = True
            # if vector >= 1:
            #     badSymm = True
            partOfSymmOp = [self.matrix[k, 0], self.matrix[k, 1], self.matrix[k, 2], vector]
            # print 'this is part of symmOp: ', partOfSymmOp
            symmline = ['{:+}{}'.format(j, ref[i]) if j else '' for i, j in enumerate(partOfSymmOp)]
            # print 'This is symmline: ', symmline
            # symmline = [block.replace('+1', '').replace('1', '') if i < 3 else block for i,
            #                                                                              block in enumerate(symmline)]
            symmline = [block.replace('1', '') if i < 3 else block for i, block in enumerate(symmline)]
            # print 'This is symm line again: ', symmline
            allsymmlines.append(''.join(symmline).lstrip('+'))
        if not all([el == ref[i] for i, el in enumerate(allsymmlines)]):  # and '1.0' not in allsymmlines[-1]:
            # if not badSymm:
            return 'SYMM {}, {}, {}\n'.format(*allsymmlines)
            # else:
            #     return ''
        else:
            return ''


spagsy = {'P1': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],  # sg1
          'A1': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],
          'B1': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],
          'C1': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],
          'I1': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],
          'F1': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],
          'P2': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (1) P 1 2 1 sg3
                 Generator(matrix([[-1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, -1]]), array([0, 0, 0]))],  # (2) 2 0,y,0
          'P21': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 1 21 1 sg4
                  Generator(matrix([[-1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, -1]]), array([0, 0.5, 0]))],  # (2) 2(0,0.5,0) 0,y,0
          'C2': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (1) C 1 2 1  is t(0.5,0.5,0) included? sg5
                 Generator(matrix([[-1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, -1]]), array([0, 0, 0]))],  # (2) 2 0,y,0
          'I2': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (1) C 1 2 1 mit t(0.5, 0.5, 0.5) sg5 ?
                 Generator(matrix([[-1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, -1]]), array([0, 0, 0]))],  # (2) 2 0,y,0
          'P222': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) P 1 2 1 sg16
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (3) 2 0,y,0
          'P2221': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # P 2 2 21 sg17
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0.5, 0, 0]))],  # (3) 2 0,y,0.25
          'P2122': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # P 21 2 2
                    Generator(matrix([[1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, -1]]), array([0.5, 0, 0])),  # (2) 2(0.5,0,0) x,0,0
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0.5, 0, 0]))],  # (3) 2 0.25,0,z
          'A2122': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # A 21 2 2
                    Generator(matrix([[1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, -1]]), array([0.5, 0, 0])),  # (2) 2(0.5,0,0) x,0,0
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0.5, 0, 0]))],  # (3) 2 0.25,0,z
          'P2212': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # P 2 21 2
                    Generator(matrix([[1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, -1]]), array([0, 0.5, 0])),  # (2) 2(0,0.5,0) x,0,0
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0.5, 0]))],  # (3) 2 0,y+0.25,0
          'B2212': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # B 2 21 2
                    Generator(matrix([[1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, -1]]), array([0, 0.5, 0])),  # (2) 2(0,0.5,0) x,0,0
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0.5, 0]))],  # (3) 2 0,y+0.25,0
          'P21212': [Generator(matrix([[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # P 21 21 2 sg18
                     Generator(matrix([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # (2) 2(0,0,0) 0,0,z
                     Generator(matrix([[-1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, -1]]), array([0.5, 0.5, 0]))],  # (3) 2(0,0.5,0) 0.25,y,0
          'P22121': [Generator(matrix([[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # P 2 21 21
                     Generator(matrix([[1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]), array([0, 0, 0])),  # (2) 2 x,0,0
                     Generator(matrix([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, 1]]), array([0, 0.5, 0.5]))],  # (3) 2(0,0,0.5) 0,0.25,z
          'A222': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # A 2 2 2
                   Generator(matrix([[1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, -1]]), array([0, 0, 0])),  # (2) 2 x,0,0
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0]))],  # (3) 2 0,0,z
          'P21221': [Generator(matrix([[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # P 21 2 21
                     Generator(matrix([[-1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, -1]]), array([0, 0, 0])),  # (2) 2 0,y,0
                     Generator(matrix([[1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, -1]]), array([0.5, 0, 0.5]))],  # (3) 2(0.5,0,0) x,0,0.25
          'B222': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # B 2 2 2
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0])),  # (2) 2 0,y,0
                   Generator(matrix([[1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (3) 2 x,0,0
          'P212121': [Generator(matrix([[1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, 1]]), array([0, 0, 0])),  # P 21 21 21 sg19
                      Generator(matrix([[-1, 0, 0],
                                        [0, -1, 0],
                                        [0, 0, 1]]), array([0.5, 0, 0.5])),  # (2) 2(0,0,0.5) 0.25,0,z
                      Generator(matrix([[-1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, -1]]), array([0, 0.5, 0.5]))],  # (3) 2(0,0.5,0) 0,y,0.25
          'C2221': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # C 2 2 21 sg20
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0, 0.5]))],  # (3) 2 0,y,0.25
          'C222': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # C 2 2 2 sg21
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (3) 2 0,y,0
          'F222': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # F 2 2 2 sg22
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (3) 2 0,y,0
          'I222': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # I 2 2 2 sg23
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (3) 2 0,y,0
          'I212121': [Generator(matrix([[1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, 1]]), array([0, 0, 0])),  # I 21 21 21 sg24
                      Generator(matrix([[-1, 0, 0],
                                        [0, -1, 0],
                                        [0, 0, 1]]), array([0.5, 0, 0.5])),  # (2) 2(0,0,0.5) 0.25,0,z
                      Generator(matrix([[-1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, -1]]), array([0, 0.5, 0.5]))],  # (3) 2(0,0.5,0) 0,y,0.25
          'P4': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # P 2 sg75 George has just one additional matrix?
                 Generator(matrix([[-1, 0, 0],
                                   [0, -1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                 Generator(matrix([[0, -1, 0],
                                   [1, 0, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],  # (3) 4+ 0,0,z
          'P41': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # P 41 sg76 George has just one additional matrix?
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                  Generator(matrix([[0, -1, 0],
                                    [1, 0, 0],
                                    [0, 0, 1]]), array([0, 0, 0.25]))],  # (3) 4+(0,0,0.25) 0,0,z
          'P42': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # P 43 sg77 George has just one additional matrix?
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (2) 2(0,0,0) 0,0,z
                  Generator(matrix([[0, -1, 0],
                                    [1, 0, 0],
                                    [0, 0, 1]]), array([0, 0, 0.5]))],  # (3) 4+(0,0,0.5) 0,0,z
          'P43': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # P 43 sg78 George has just one additional matrix?
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                  Generator(matrix([[0, -1, 0],
                                    [1, 0, 0],
                                    [0, 0, 1]]), array([0, 0, 0.75]))],  # (3) 4+(0,0,0.75) 0,0,z
          'I4': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # I 2 sg79 George has just one additional matrix?
                 Generator(matrix([[-1, 0, 0],
                                   [0, -1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                 Generator(matrix([[0, -1, 0],
                                   [1, 0, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],  # (3) 4+ 0,0,z
          'I41': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # I 41 sg80 George has just one additional matrix?
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0.5, 0.5, 0.5])),  # (2) 2(0,0,0.5) 0.25,0.25,z
                  Generator(matrix([[0, -1, 0],
                                    [1, 0, 0],
                                    [0, 0, 1]]), array([0, 0.5, 0.25]))],  # (3) 4+(0,0,0.25) -0.25,0.25,z
          'P422': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # P 422 sg89
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[0, -1, 0],
                                     [1, 0, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (3) 4+ 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (5) 2 0,y,0
          'P4212': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # P 4 21 2 sg90
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                    Generator(matrix([[0, -1, 0],
                                      [1, 0, 0],
                                      [0, 0, 1]]), array([0.5, 0.5, 0])),  # (3) 4+ 0,0.5,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0.5, 0.5, 0]))],  # (5) 2(0,0.5,0) 0.25,y,0
          'P4122': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # P 41 2 2 sg91
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                    Generator(matrix([[0, -1, 0],
                                      [1, 0, 0],
                                      [0, 0, 1]]), array([0, 0, 0.25])),  # (3) 4+(0,0,0.25) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (5) 2 0,y,0
          'P41212': [Generator(matrix([[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # P 41 21 2 sg92
                     Generator(matrix([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                     Generator(matrix([[0, -1, 0],
                                       [1, 0, 0],
                                       [0, 0, 1]]), array([0.5, 0.5, 0.25])),  # (3) 4+(0,0,0.25) 0,0.5,z
                     Generator(matrix([[-1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, -1]]), array([0.5, 0.5, 0.25]))],  # (5) 2(0,0.5,0) 0.25,y,0.125
          'P4222': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # P 42 2 2 sg93
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (2) 2(0,0,0) 0,0,z
                    Generator(matrix([[0, -1, 0],
                                      [1, 0, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (3) 4+(0,0,0.5) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (5) 2 0,y,0
          'P42212': [Generator(matrix([[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # P 42 21 2 sg94
                     Generator(matrix([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # (2) 2(0,0,0) 0,0,z
                     Generator(matrix([[0, -1, 0],
                                       [1, 0, 0],
                                       [0, 0, 1]]), array([0.5, 0.5, 0.5])),  # (3) 4+(0,0,0.5) 0,0.5,z
                     Generator(matrix([[-1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, -1]]), array([0.5, 0.5, 0.5]))],  # (5) 2(0,0.5,0) 0.25,y,0.25
          'P4322': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # P 43 2 2 sg95
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                    Generator(matrix([[0, -1, 0],
                                      [1, 0, 0],
                                      [0, 0, 1]]), array([0, 0, 0.75])),  # (3) 4+(0,0,0.75) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (5) 2 0,y,0
          'P43212': [Generator(matrix([[1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, 1]]), array([0, 0, 0])),  # P 43 21 2 sg96
                     Generator(matrix([[-1, 0, 0],
                                       [0, -1, 0],
                                       [0, 0, 1]]), array([0, 0, 0.5])),  # (2) 2(0,0,0.5) 0,0,z
                     Generator(matrix([[0, -1, 0],
                                       [1, 0, 0],
                                       [0, 0, 1]]), array([0.5, 0.5, 0.75])),  # (3) 4+(0,0,0.75) 0,0.5,z
                     Generator(matrix([[-1, 0, 0],
                                       [0, 1, 0],
                                       [0, 0, -1]]), array([0.5, 0.5, 0.75]))],  # (5) 2(0,0.5,0) 0.25,y,0.375
          'I422': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # I 4 2 2 sg97
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[0, -1, 0],
                                     [1, 0, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (3) 4+ 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (5) 2 0,y,0
          'I4122': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # I 41 2 2 sg98
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0.5, 0.5, 0.5])),  # (2) 2(0,0,0.5) 0.25,0.25,z
                    Generator(matrix([[0, -1, 0],
                                      [1, 0, 0],
                                      [0, 0, 1]]), array([0, 0.5, 0.25])),  # (3) 4+(0,0,0.25) -0.25,0.25,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0.5, 0, 0.75]))],  # (5) 2 0.25,y,0.375
          'P3': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (1) P 3 sg143
                 Generator(matrix([[0, -1, 0],
                                   [1, -1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],  # (2) 3+ 0,0,z
          'P31': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 31 sg144
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 1. / 3.]))],  # (2) 2(0,0,1/3) 0,0,z
          'P32': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 32 sg145
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 2. / 3.]))],  # (2) 2(0,0,2/3) 0,0,z
          'R3': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (1) R 3 sg146 rhombohedral axes
                 Generator(matrix([[0, 0, 1],
                                   [1, 0, 0],
                                   [0, 1, 0]]), array([0, 0, 0]))],  # (2) 3+ x,x,x
          'H3': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (1) R 3 sg146 hexagonal axes
                 Generator(matrix([[0, -1, 0],
                                   [1, -1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],  # (2) 3+ 0,0,z
          'P312': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) P 3 1 2 sg149
                   Generator(matrix([[0, -1, 0],
                                     [1, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 3+ 0,0,z
                   Generator(matrix([[0, -1, 0],
                                     [-1, 0, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (4) 2 x,-x,0
          'P321': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) P 3 2 1 sg150
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (2) 3+ 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (4) 2 x,x,0
          'P3112': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 31 1 2 sg151
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 1. / 3.])),  # (2) 3+(0,0,1/3) 0,0,z
                    Generator(matrix([[0, -1, 0],
                                      [-1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 2. / 3.]))],  # (4) 2(0,0,1/3) x,-x,0
          'P3121': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 31 2 1 sg152
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 1. / 3.])),  # (2) 3+(0,0,1/3) 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (4) 2 x,x,0
          'P3212': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 32 1 2 sg153
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 2. / 3.])),  # (2) 3+(0,0,2/3) 0,0,z
                    Generator(matrix([[0, -1, 0],
                                      [-1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 1. / 3.]))],  # (4) 2 x,-x,1/6
          'P322': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 32 2 1 sg154
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 2. / 3.])),  # (2) 3+(0,0,2/3) 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (4) 2 x,x,0
          'P3221': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 32 1 2 sg154
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 2. / 3.])),  # (2) 3+(0,0,2/3) 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (4) 2 x,x,0
          'R32': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) R 32 sg155 rhombohedral axes
                  Generator(matrix([[0, 0, 1],
                                    [1, 0, 0],
                                    [0, 1, 0]]), array([0, 0, 0])),  # (2) 3+ x,x,x
                  Generator(matrix([[0, 0, -1],
                                    [0, -1, 0],
                                    [-1, 0, 0]]), array([0, 0, 0]))],  # (4) 2 -x,0,x (anders in den alten Tables!!!!)
          'H32': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) H 32  hexagonal axes (sg155)
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (2) 3+ 0,0,z
                  Generator(matrix([[0, 1, 0],
                                    [1, 0, 0],
                                    [0, 0, -1]]), array([0, 0, 0]))],  # (4) 2 x,x,0
          'P6': [Generator(matrix([[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (1) P 6 sg168
                 Generator(matrix([[0, -1, 0],
                                   [1, -1, 0],
                                   [0, 0, 1]]), array([0, 0, 0])),  # (2) 3+ 0,0,z
                 Generator(matrix([[-1, 0, 0],
                                   [0, -1, 0],
                                   [0, 0, 1]]), array([0, 0, 0]))],  # (4) 2 0,0,z
          'P61': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 61 sg169
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 1. / 3.])),  # (2) 3+(0,0,1/3) 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0.5]))],  # (4) 2(0,0,0.5) 0,0,z
          'P65': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 65 sg170
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 2. / 3.])),  # (2) 3+(0,0,2/3) 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0.5]))],  # (4) 2(0,0,0.5) 0,0,z
          'P62': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 62 sg171
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 2. / 3.])),  # (2) 3+(0,0,2/3) 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0]))],  # (4) 2 0,0,z
          'P64': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 64 sg172
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 1. / 3.])),  # (2) 3+(0,0,1/3) 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0]))],  # (4) 2 0,0,z
          'P63': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P 63 sg173
                  Generator(matrix([[0, -1, 0],
                                    [1, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (2) 3+ 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0.5]))],  # (4) 2(0,0,0.5) 0,0,z
          'P622': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) P 6 2 2 sg177
                   Generator(matrix([[0, -1, 0],
                                     [1, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 3+ 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (4) 2 0,0,z
                   Generator(matrix([[0, 1, 0],
                                     [1, 0, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (7) 2 x,x,0
          'P6122': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 61 2 2 sg178
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 1. / 3.])),  # (2) 3+(0,0,1/3) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (4) 2(0,0,0.5) 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 1. / 3.]))],  # (7) 2 x,x,1/6
          'P6522': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 65 2 2 sg179
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 2. / 3.])),  # (2) 3+(0,0,2/3) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (4) 2(0,0,0.5) 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 2. / 3.]))],  # (7) 2 x,x,1/3
          'P6222': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 62 2 2 sg180
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 2. / 3.])),  # (2) 3+(0,0,2/3) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (4) 2 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 2. / 3.]))],  # (7) 2 x,x,1/3
          'P6422': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 64 2 2 sg181
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 1. / 3.])),  # (2) 3+(0,0,1/3) 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (4) 2 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 1. / 3.]))],  # (7) 2 x,x,1/6
          'P6322': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P 63 2 2 sg182
                    Generator(matrix([[0, -1, 0],
                                      [1, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (2) 3+ 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0.5])),  # (4) 2(0,0,0.5) 0,0,z
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0, 0, 0]))],  # (7) 2 x,x,0
          'P23': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) P23 sg195
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, -1]]), array([0, 0, 0])),  # (3) 2 0,y,0
                  Generator(matrix([[0, 0, 1],
                                    [1, 0, 0],
                                    [0, 1, 0]]), array([0, 0, 0]))],  # (5) 3+ x,x,x
          'F23': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) F23 sg196
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, -1]]), array([0, 0, 0])),  # (3) 2 0,y,0
                  Generator(matrix([[0, 0, 1],
                                    [1, 0, 0],
                                    [0, 1, 0]]), array([0, 0, 0]))],  # (5) 3+ x,x,x
          'I23': [Generator(matrix([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (1) I23 sg197
                  Generator(matrix([[-1, 0, 0],
                                    [0, -1, 0],
                                    [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                  Generator(matrix([[-1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, -1]]), array([0, 0, 0])),  # (3) 2 0,y,0
                  Generator(matrix([[0, 0, 1],
                                    [1, 0, 0],
                                    [0, 1, 0]]), array([0, 0, 0]))],  # (5) 3+ x,x,x
          'P213': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) P21 3 sg198
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0.5, 0, 0.5])),  # (2) 2(0,0,0.5) 0.25,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0.5, 0.5])),  # (3) 2(0,0.5,0) 0,y,0.25
                   Generator(matrix([[0, 0, 1],
                                     [1, 0, 0],
                                     [0, 1, 0]]), array([0, 0, 0]))],  # (5) 3+ x,x,x
          'I213': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) I21 3 sg199
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0.5, 0.5])),  # (3) 2(0,0.5,0) 0,y,0.25
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0.5, 0, 0.5])),  # (2) 2(0,0,0.5) 0.25,0,z
                   Generator(matrix([[0, 0, 1],
                                     [1, 0, 0],
                                     [0, 1, 0]]), array([0, 0, 0]))],  # (5) 3+ x,x,x
          'P432': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) P4 3 2 sg207
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0])),  # (3) 2 0,y,0
                   Generator(matrix([[0, 0, 1],
                                     [1, 0, 0],
                                     [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                   Generator(matrix([[0, 1, 0],
                                     [1, 0, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (13) 2 x,x,0
          'P4232': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P42 3 2 sg208
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0, 0])),  # (3) 2 0,y,0
                    Generator(matrix([[0, 0, 1],
                                      [1, 0, 0],
                                      [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0.5, 0.5, 0.5]))],  # (13) 2(0.5,0.5,0) x,x,0.25
          'F432': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) I4 3 2 sg209
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0])),  # (3) 2 0,y,0
                   Generator(matrix([[0, 0, 1],
                                     [1, 0, 0],
                                     [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                   Generator(matrix([[0, 1, 0],
                                     [1, 0, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (13) 2 x,x,0
          'F4132': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) F41 3 2 sg210
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0, 0.5, 0.5])),  # (2) 2(0,0,0.5) 0,0.25,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0.5, 0.5, 0])),  # (3) 2(0,0.5,0) 0.25,y,0
                    Generator(matrix([[0, 0, 1],
                                      [1, 0, 0],
                                      [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0.75, 0.25, 0.75]))],  # (13) 2(0.5,0.5,0) x,x-0.25,3/8
          'I432': [Generator(matrix([[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (1) I4 3 2 sg209
                   Generator(matrix([[-1, 0, 0],
                                     [0, -1, 0],
                                     [0, 0, 1]]), array([0, 0, 0])),  # (2) 2 0,0,z
                   Generator(matrix([[-1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, -1]]), array([0, 0, 0])),  # (3) 2 0,y,0
                   Generator(matrix([[0, 0, 1],
                                     [1, 0, 0],
                                     [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                   Generator(matrix([[0, 1, 0],
                                     [1, 0, 0],
                                     [0, 0, -1]]), array([0, 0, 0]))],  # (13) 2 x,x,0
          'P4332': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P43 3 2 sg212
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0.5, 0, 0.5])),  # (2) 2(0,0,0.5) 0.25,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0.5, 0.5])),  # (3) 2(0,0.5,0) 0,y,0.25
                    Generator(matrix([[0, 0, 1],
                                      [1, 0, 0],
                                      [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([0.25, 0.75, 0.75]))],  # (13) 2(0.5,0.5,0) x,x+0.25,3/8
          'P4132': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) P41 3 2 sg213
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0.5, 0, 0.5])),  # (2) 2(0,0,0.5) 0.25,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0.5, 0.5])),  # (3) 2(0,0.5,0) 0,y,0.25
                    Generator(matrix([[0, 0, 1],
                                      [1, 0, 0],
                                      [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([-0.25, 0.25, 0.25]))],  # (13) 2(0.5,0.5,0) x,x-0.25,1/8
          'I4132': [Generator(matrix([[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]]), array([0, 0, 0])),  # (1) I41 3 2 sg214
                    Generator(matrix([[-1, 0, 0],
                                      [0, -1, 0],
                                      [0, 0, 1]]), array([0.5, 0, 0.5])),  # (2) 2(0,0,0.5) 0.25,0,z
                    Generator(matrix([[-1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, -1]]), array([0, 0.5, 0.5])),  # (3) 2(0,0.5,0) 0,y,0.25
                    Generator(matrix([[0, 0, 1],
                                      [1, 0, 0],
                                      [0, 1, 0]]), array([0, 0, 0])),  # (5) 3+ x,x,x
                    Generator(matrix([[0, 1, 0],
                                      [1, 0, 0],
                                      [0, 0, -1]]), array([-0.25, 0.25, 0.25]))],  # (13) 2(0.5,0.5,0) x,x-0.25,1/8


}






