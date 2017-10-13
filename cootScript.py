__author__ = 'Anna'
"""
custom coot script.
-calling coot from console
-open shelx project
-display/generate electron density maps
-set correct map settings
-export difference density map
-exit coot
"""

rootdir = '/home/anna/pdbFiles/pdbFilesTill17Res'

filename = rootdir + folder + filename
#
allow_duplicate_sequence_numbers()

#
read_shelx_ins_file(filenameRes)

#
read_shelx_fcf_file(filenameFCF)

#
handle_shelx_fcf_file(filenameFCF)

#
display_maps(1)

#
map_is_difference_map(3)

#
set_diff_map_iso_level_increment(0.001)

#
set_last_map_contour_level_0.05

#
set_map_is_difference_map(3)

#
eport_map(3, FullExportPath)

#
coot_real_exit(0)

