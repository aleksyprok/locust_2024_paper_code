"""
This file reads in a CDB file and outputs a LOCUST wall file.
"""

import os
import time
import matplotlib.pyplot as plt

def prepare_locust_wall(cdb_fname, locust_fname):
    """
    This function reads in a CDB file and outputs a LOCUST wall file.
    """

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
                # locust_file.write(fmt.format(cmpnt_no, cmpnt_no, 1, cmpnt_no, 0, 0, 0, 0, 8, 0,
                #                              *split_line[10:14], *duplicates))
                locust_file.write(fmt.format(cmpnt_no, cmpnt_no, 1, cmpnt_no, 0, 0, 0, 0, 8, 0,
                                             tet_counter, *split_line[11:14], *duplicates))

    # End EBLOCK processing
    locust_file.write('/wb,elem,end\n')

    cdb_file.close()
    locust_file.close()


def prepare_continuous_wall(csv_fname, dat_fname):
    with open(csv_fname, 'r') as file:
        lines = file.readlines()
    
    # Initialize variables for the different sections
    main_wall_coords = []
    lower_dome_coords = []
    upper_dome_coords = []
    current_coords = None

    # Parse the file line by line
    for line in lines:
        line = line.strip()
        if not line:  # Skip empty lines
            continue

        if line == "Main wall":
            current_coords = main_wall_coords
        elif line == "Lower Dome":
            current_coords = lower_dome_coords
        elif line == "Upper Dome":
            current_coords = upper_dome_coords
        elif "," in line and "R,Z" not in line:  # Coordinate line
            r, z = map(float, line.split(","))
            current_coords.append((r / 1000, z / 1000))  # Convert to meters

    # Plot each section
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_section(ax, "Main Wall", main_wall_coords)
    plot_section(ax, "Lower Dome", lower_dome_coords)
    plot_section(ax, "Upper Dome", upper_dome_coords)

    # Set plot labels and legend
    ax.set_xlabel("R (mm)")
    ax.set_ylabel("Z (mm)")
    ax.legend()
    ax.set_title("Continuous Wall and Domes")
    ax.grid(True)
    ax.set_aspect("equal")

    continuous_wall_coords = []
    continuous_wall_coords.append((main_wall_coords[0][0], 0))
    continuous_wall_coords.extend(main_wall_coords[0:7])
    # continuous_wall_coords.append((2.5366, -8.1801))
    continuous_wall_coords.append((2.68688, -8.22998))
    continuous_wall_coords.append((2.73642, -8.07076))
    continuous_wall_coords.extend(lower_dome_coords[0:5])
    continuous_wall_coords.append((2.737563, -8.071124))
    continuous_wall_coords.append((2.687022, -8.230092))
    continuous_wall_coords.extend(main_wall_coords[7:41])
    continuous_wall_coords.append((2.751634, 8.251497))
    continuous_wall_coords.append((2.792549, 8.089293))
    continuous_wall_coords.extend(upper_dome_coords[0:5][::-1])
    continuous_wall_coords.append((2.79028, 8.08725))
    continuous_wall_coords.append((2.751378, 8.251434))
    continuous_wall_coords.extend((main_wall_coords[41:]))
    continuous_wall_coords.append(continuous_wall_coords[0])

    print(continuous_wall_coords)

    plot_section(ax, "Continuous Wall", continuous_wall_coords)

    with open(dat_fname, 'w') as dat_file:
        for r, z in continuous_wall_coords:
            dat_file.write(f"{r:8.6f} {z:8.6f}\n")

    plt.show()

def plot_section(ax, name, coords):
    """Helper function to plot a section with annotated indices."""
    r_values, z_values = zip(*coords)
    ax.plot(r_values, z_values, marker='o', label=name)
    for idx, (r, z) in enumerate(coords):
        ax.annotate(str(idx), (r, z), textcoords="offset points", xytext=(5, -5), fontsize=8)

if __name__ == "__main__":
    tic = time.time()
    REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    CDB_FNAME = os.path.join(REPOSITORY_PATH, "input_data", "Stepmesh.cdb")
    LOCUST_FNAME = os.path.join(REPOSITORY_PATH, "input_data", "SPR-068.cdb.locust")
    prepare_locust_wall(CDB_FNAME, LOCUST_FNAME)
    CSV_FNAME = os.path.join(REPOSITORY_PATH, "input_data", "SPR68_2D_Wall.csv")
    DAT_FNAME = os.path.join(REPOSITORY_PATH, "input_data", "SPR-068_wall.dat")
    prepare_continuous_wall(CSV_FNAME, DAT_FNAME)
    toc = time.time()
    print(f'Time taken to process: {toc - tic:.2f} seconds')