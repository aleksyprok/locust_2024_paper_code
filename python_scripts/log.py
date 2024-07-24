"""
This module contains routines for updating a Run object from the
runs.py module with information from a LOCUST LOG*.out file.
"""
from typing import Optional

class AttributeParser:
    """
    This class is used to parse a line from a LOCUST .log file and extract
    an attribute from it.
    """
    def __init__(self, attribute_name, parse_action):
        self.attribute_name = attribute_name
        self.parse_action = parse_action

class Log:
    """
    This class contains routines for updating a Run object from the
    runs.py module with information from a LOCUST .log file.
    """
    def __init__(self, log_path):
        self.analytic_ripple: Optional[bool] = None
        self.bplasma_file: Optional[str] = None
        self.bplasma: Optional[bool] = None
        self.axisymmetric: Optional[bool] = None
        self.ncoil: Optional[int] = None
        self.rcoil: Optional[float] = None
        self.rcoil_inner: Optional[float] = None
        self.bplasma_n: Optional[int] = None
        self.coil_set: Optional[str] = None
        self.rmp_current: Optional[float] = None
        self.rmp_phase: Optional[float] = None
        self.rmp_response: Optional[bool] = None
        self.rwm_gain: Optional[float] = None
        self.rwm_bscale: Optional[float] = None
        # total_stopped_power_locust, total_stopped_power_locust_error, \
        # max_energy_flux_locust, Pinj
        self.total_stopped_power_string: Optional[str] = None
        self.total_stopped_power: Optional[float] = None
        self.total_stopped_power_error: Optional[float] = None
        self.pinj: Optional[float] = None # Total injected power [MW]
        self.simulation_time = None
        self.eqdsk_fname: Optional[str] = None
        self.i3dr: Optional[int] = None
        parsing_map = {
            ':locust_info : analytical TF ripple field':
                AttributeParser('analytic_ripple',
                                lambda line: line.split()[-1] == 'ENABLED'),
            ':locust_bplasm_cache : file  :':
                AttributeParser('bplasma_file', lambda line: line.split()[-1]),
            ':locust_info : TF_Ncoil':
                AttributeParser('ncoil', lambda line: int(line.split()[-1])),
            ':locust_info : TF_Rcoil =':
                AttributeParser('rcoil', lambda line: float(line.split()[-1])),
            ':locust_info : TF_Rcoil_inner':
                AttributeParser('rcoil_inner', lambda line: float(line.split()[-1])),
            " :write_vtk   : Integrated power to PFCs           :":
                AttributeParser('total_stopped_power_string',
                                lambda line: " ".join(line.split()[-3:])),
            "Ensemble power":
                AttributeParser('pinj',
                                lambda line: float(line.split()[-1][:-2])),
            "....... finishing ....... : time":
                AttributeParser('simulation_time',
                                lambda line: float(line.split()[-5][:-1])),
            "locust_geqdsk : Attempt to open":
                AttributeParser('eqdsk_fname', lambda line: line.split('/')[-1][:-1]),
            ":locust_info : 3D field in":
                AttributeParser('i3dr',
                                lambda line: int('LOCUST' in line) -
                                             int('MARS-F' in line)),
        }

        self.log_path = log_path
        self.parse_log(log_path, parsing_map)
        self.post_parse()

    def parse_log(self, log_path, parsing_map):
        """
        Open and read log file and update the Run object with the parsed
        information.
        """
        with open(log_path, 'r', encoding='utf-8') as file:
            for line in file:
                for distinctive_string, parser in parsing_map.items():
                    if distinctive_string in line:
                        parsed_value = parser.parse_action(line)
                        setattr(self, parser.attribute_name, parsed_value)

    def post_parse(self):
        """
        After reading the log file, we can determine additional attributes
        based on the parsed information.
        """
        if self.bplasma_file is not None:
            self.bplasma = True
        else:
            self.bplasma = False
        if self.bplasma is False and self.analytic_ripple is False:
            self.axisymmetric = True
        else:
            self.axisymmetric = False
        if self.bplasma_file is not None:
            self.bplasma_n = int(self.bplasma_file[-1])
            if 'efcc' in self.bplasma_file:
                self.coil_set = 'exterior_rmp'
            elif 'rwm' in self.bplasma_file:
                self.coil_set = 'interior_rmp'
            elif self.bplasma_n == 1:
                self.coil_set = 'rwm_acc'
            else:
                raise ValueError("Unknown coil set encountered.")
            if 'current=' in self.bplasma_file:
                self.rmp_current = float(self.bplasma_file.split('current=')[1].split('_')[0])
            if 'phase=' in self.bplasma_file:
                self.rmp_phase = float(self.bplasma_file.split('phase=')[1].split('_')[0])
            if 'response=' in self.bplasma_file:
                self.rmp_response = int(self.bplasma_file.split('response=')[1].split('_')[0]) == 1
            if 'G=' in self.bplasma_file:
                self.rwm_gain = float(self.bplasma_file.split('G=')[1].split('_')[0])
            if 'bscale=' in self.bplasma_file:
                self.rwm_bscale = float(self.bplasma_file.split('bscale=')[1].split('_')[0])
        if self.total_stopped_power_string is None:
            self.total_stopped_power = None
        elif self.total_stopped_power_string.split()[1] == '+-':
            self.total_stopped_power = float(self.total_stopped_power_string.split()[0])
            self.total_stopped_power_error = float(self.total_stopped_power_string.split()[2][:-2])
        else:
            self.total_stopped_power = float(self.total_stopped_power_string.split()[2][:-2])
