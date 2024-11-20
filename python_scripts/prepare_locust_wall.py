"""
This file reads in a CDB file and outputs a LOCUST wall file.
"""

import os
import time

tic = time.time()

repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cdb_fname = os.path.join(repository_path, "input_data", "Stepmesh.cdb")
locust_fname = os.path.join(repository_path, "input_data", "SPR-068.cdb.locust")

tet_counter = 0

cdb_file = open(cdb_fname, 'r', encoding="utf-8")
locust_file = open(locust_fname, 'w', encoding="utf-8")

# Process NBLOCK
for line in cdb_file:
    if 'NBLOCK' in line:
        print('\nNBLOCK found.')
        locust_file.write('nblock\n')
        break

print('CDB NBLOCK (Alleged) Format = ' + cdb_file.readline().strip('\n'))
fmt = '{:8d}' + 3 * '{:20.9e}' + '\n'
print('LOCUST NBLOCK Format = ' + fmt.strip('\n'))
locust_file.write('(1i8,3e20.9e3)\n')

for line in cdb_file:
    split_line = line.split()
    if '-1,' in split_line[-1]:
        print('End of NBLOCK reached.')
        locust_file.write('-1\n')
        break
    else:
        locust_file.write(fmt.format(int(split_line[0]),
                                     float(split_line[-3]),
                                     float(split_line[-2]),
                                     float(split_line[-1])))

# Process multiple EBLOCKs
cmpnt_no = -1
fmt = 19 * '{:9d}' + '\n'

while True:
    for line in cdb_file:
        if 'EBLOCK' in line:
            cmpnt_no += 1
            print(f'\nStart of EBLOCK {cmpnt_no} reached.')
            locust_file.write('/wb,elem,start\n')
            locust_file.write(f'/com,*********** Elements for Body {cmpnt_no} ***********\n')
            locust_file.write('eblock\n')
            break
    else:
        # Break the outer loop if no more EBLOCKs are found
        break

    print('CDB EBLOCK (Alleged) Format = ' + cdb_file.readline().strip('\n'))
    print('LOCUST EBLOCK Format = ' + fmt.strip('\n'))
    locust_file.write('(19i9)\n')

    for line in cdb_file:
        if '-1' in line:
            print(f'End of EBLOCK {cmpnt_no} reached.')
            locust_file.write('-1\n')
            break
        else:
            split_line = list(map(int, line.split()))
            tet_counter += 1
            duplicates = [split_line[-1]] * 5
            locust_file.write(fmt.format(cmpnt_no, cmpnt_no, 1, cmpnt_no, 0, 0, 0, 0, 8, 0,
                                         *split_line[10:14], *duplicates))

# End EBLOCK processing
locust_file.write('/wb,elem,end\n')

cdb_file.close()
locust_file.close()

toc = time.time()
print(f'Time taken to process: {toc - tic:.2f} seconds')
