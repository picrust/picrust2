#!/usr/bin/python

# ==================================================================================================
#     Flags
# ==================================================================================================
tgl_edge_num = True
tgl_likelihood = True
tgl_lwr = True
tgl_proximal = False
tgl_pendant = False

# ==================================================================================================
#     Script Init.
# ==================================================================================================

import sys
from genesis import *

utils.Logging.log_to_stdout()

if __name__ != "__main__":
    print "Not running as main program."
    exit()

if len(sys.argv) != 4:
    mode = "wrong"
else:
    mode = sys.argv[1]

if not mode in [ "-p", "-v", "-s" ]:
    print "Usage:"
    print "       " + sys.argv[0], "-[pvs] file1 file2\n"
    print "with options:"
    print "       -p: Print both jplace file contents as tables."
    print "       -v: Validate that they are equal. Output a comparison table."
    print "       -s: Validate that they are equal. No output (except warnings etc).\n"
    print "The return value of the script is:"
    print "        0: Files are equal."
    print "        1: Files are not equal."
    exit()

jplace_lhs = sys.argv[2]
jplace_rhs = sys.argv[3]

# ==================================================================================================
#     Prepare: Load and check files.
# ==================================================================================================

# Load files.

jpp = placement.JplaceReader()
jpp.report_invalid_numbers = True

pmap_lhs = placement.SampleSet()
pmap_rhs = placement.SampleSet()

jpp.from_file(jplace_lhs, pmap_lhs)
jpp.from_file(jplace_rhs, pmap_rhs)

# Do some checks.

if not pmap_lhs.validate() or not pmap_rhs.validate():
    print "PlacementMaps not valid."
    exit()

if not placement.has_correct_edge_nums(pmap_lhs) or not placement.has_correct_edge_nums(pmap_rhs):
    print "Edge nums not correct."
    exit()

if not placement.compatible_trees( pmap_lhs, pmap_rhs ):
    print "Trees not compatible."
    exit()

# Check sizes.

if pmap_lhs.pquery_size() != pmap_rhs.pquery_size():
    print "Different number of pqueries."
    exit()

# Check if both contain the same names.
# Also, for simplicity, we expect that each Pquery contains exaclty one name.

lhs_names = []
for i in range(0, pmap_lhs.pquery_size()):
    if pmap_lhs.pquery(i).name_size() != 1:
        print "Can only compare Pquries with one name each."
        exit()
    lhs_names.append( pmap_lhs.pquery(i).name_at(0).name )

rhs_names = []
for i in range(0, pmap_rhs.pquery_size()):
    if pmap_rhs.pquery(i).name_size() != 1:
        print "Can only compare Pquries with one name each."
        exit()
    rhs_names.append( pmap_rhs.pquery(i).name_at(0).name )

if set(lhs_names) != set(rhs_names):
    print "Pquery names differ."
    exit()

# ==================================================================================================
#     Compare Pqueries.
# ==================================================================================================

# Sort the Placements per Pquery by their like weight ratio.

placement.sort_placements_by_like_weight_ratio( pmap_lhs )
placement.sort_placements_by_like_weight_ratio( pmap_rhs )

# In print mode: Simply print both tables.

if mode == "-p":
    print pmap_lhs
    print pmap_rhs
    exit()

# Prepare the output table with columns.

if mode == "-v":
    tab = utils.text.Table()
    tab.add_column("Rank")
    tab.add_column("Name")

    if tgl_edge_num:
        tab.add_column("edge_num L")
        tab.add_column("edge_num R")

    if tgl_likelihood:
        tab.add_column("likelihood L")
        tab.add_column("likelihood R")

    if tgl_lwr:
        tab.add_column("like_weight_ratio L")
        tab.add_column("like_weight_ratio R")

    if tgl_proximal:
        tab.add_column("proximal_length L")
        tab.add_column("proximal_length R")

    if tgl_pendant:
        tab.add_column("pendant_length L")
        tab.add_column("pendant_length R")

    tab.add_column("Correct?")

# Iterate all pqueries and compare them!

status = 0
correct = 0
failed = 0
for i in range(0, pmap_lhs.pquery_size()):

    # Get the Pqueries

    pqry_lhs = pmap_lhs.pquery(i)
    pqry_rhs = placement.find_pquery( pmap_rhs, pqry_lhs.name_at(0).name )
    rank = 0

    for j in range(min( pqry_lhs.placement_size(), pqry_rhs.placement_size() )):
        place_lhs = pqry_lhs.placement_at(j)
        place_rhs = pqry_rhs.placement_at(j)

        rank += 1

        # Fill the table.

        if mode == "-v":
            tab.append(str(rank))
            tab.append(pqry_lhs.name_at(0).name)

            if tgl_edge_num:
                tab.append(str(place_lhs.edge_num))
                tab.append(str(place_rhs.edge_num))

            if tgl_likelihood:
                tab.append(str(place_lhs.likelihood))
                tab.append(str(place_rhs.likelihood))

            if tgl_lwr:
                tab.append(str(place_lhs.like_weight_ratio))
                tab.append(str(place_rhs.like_weight_ratio))

            if tgl_proximal:
                tab.append(str(place_lhs.proximal_length))
                tab.append(str(place_rhs.proximal_length))

            if tgl_pendant:
                tab.append(str(place_lhs.pendant_length))
                tab.append(str(place_rhs.pendant_length))

        # Do the comparison and act accordingly.

        if place_lhs.edge_num == place_rhs.edge_num:
            if mode == "-v":
                tab.append("\x1b[01;32mYES\x1b[00m")
            correct += 1
        else:
            if mode == "-v":
                tab.append("\x1b[01;31mNO\x1b[00m")
            status = 1
            failed += 1

# Final output: table and exit code.

total = correct + failed
percent_failed = (float(failed) / float(total)) * 100.0

if mode == "-v":
    print tab.to_string( utils.text.simple_layout() )

print "Compared {0:d} edge placements, {1:d} of which were incorrect ({2:.2f}%)".format(total, failed, percent_failed)

sys.exit(status)
