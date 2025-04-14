#unit: Angstrom for coordinates, nm for bonds, angles, dihedrals and impropers
from math import sqrt
from math import pi
import numpy as np
from collections import namedtuple
from collections import defaultdict

TEMPERATURE = 299.15
Nm_L = 4
rho_L = 3
a_ww = 102.067126235101
delta_aij_lambda = 0.28042708
gamma_dict = {20: 19.628722739512}
gamma_mix_dict = {20: 10.610342649291}

if False:
    rho_L = 20
    a_ww = 14.07677776458
    delta_aij_lambda = 3.75488652
    gamma_mix_dict = {20: 19.628722739512}


const_Rc_A = ((Nm_L*rho_L*30)**(1/3))
const_Rc_nm = const_Rc_A / 10
reciprocal_const_Rc_nm = 1 / const_Rc_nm
reciprocal_const_Rc_A = 1 / const_Rc_A

CONST_Z = 8.843629103332

RTNA_KJpermol = 1.38064852 * TEMPERATURE * 6.0221409 / 1000
k_constraints = 30000
#unit of charge
const_e = 1.60217662*10**(-19) / sqrt(8.85418782*10**(-12) * 4 * pi * const_Rc_A * 10**(-10) * 1.38064852*10**(-23) * TEMPERATURE)
martini_table = {'Qda': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qd': {'Qda': 5.6, 'Qd': 5.0, 'Qa': 5.6, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 4.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Qa': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.0, 'P1': 5.0, 'Nda': 5.0, 'Nd': 5.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'Q0': {'Qda': 4.5, 'Qd': 4.5, 'Qa': 4.5, 'Q0': 3.5, 'P5': 5.0, 'P4': 5.6, 'P3': 5.0, 'P2': 4.5, 'P1': 4.0, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.3, 'C2': 2.0, 'C1': 2.0}, 'P5': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.6, 'P3': 5.6, 'P2': 5.6, 'P1': 5.6, 'Nda': 5.0, 'Nd': 5.0, 'Na': 5.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P4': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.6, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.0, 'Nd': 4.0, 'Na': 4.0, 'N0': 3.5, 'C5': 3.1, 'C4': 2.7, 'C3': 2.7, 'C2': 2.3, 'C1': 2.0}, 'P3': {'Qda': 5.6, 'Qd': 5.6, 'Qa': 5.6, 'Q0': 5.0, 'P5': 5.6, 'P4': 5.0, 'P3': 5.0, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P2': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.5, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.1, 'C2': 2.7, 'C1': 2.3}, 'P1': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.6, 'P4': 4.5, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 4.0, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'Nda': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Nd': {'Qda': 5.0, 'Qd': 4.0, 'Qa': 5.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.0, 'Na': 4.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'Na': {'Qda': 5.0, 'Qd': 5.0, 'Qa': 4.0, 'Q0': 4.0, 'P5': 5.0, 'P4': 4.0, 'P3': 4.5, 'P2': 4.5, 'P1': 4.5, 'Nda': 4.5, 'Nd': 4.5, 'Na': 4.0, 'N0': 3.5, 'C5': 3.5, 'C4': 3.1, 'C3': 2.7, 'C2': 2.7, 'C1': 2.7}, 'N0': {'Qda': 3.5, 'Qd': 3.5, 'Qa': 3.5, 'Q0': 3.5, 'P5': 3.5, 'P4': 3.5, 'P3': 3.5, 'P2': 4.0, 'P1': 4.0, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 2.7}, 'C5': {'Qda': 3.1, 'Qd': 3.1, 'Qa': 3.1, 'Q0': 3.1, 'P5': 3.1, 'P4': 3.1, 'P3': 3.5, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.5, 'Nd': 3.5, 'Na': 3.5, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C4': {'Qda': 2.7, 'Qd': 2.7, 'Qa': 2.7, 'Q0': 2.7, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.5, 'P1': 3.5, 'Nda': 3.1, 'Nd': 3.1, 'Na': 3.1, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.1, 'C1': 3.1}, 'C3': {'Qda': 2.3, 'Qd': 2.3, 'Qa': 2.3, 'Q0': 2.3, 'P5': 2.7, 'P4': 2.7, 'P3': 3.1, 'P2': 3.1, 'P1': 3.5, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.5, 'C5': 3.5, 'C4': 3.5, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C2': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.3, 'P4': 2.3, 'P3': 2.7, 'P2': 2.7, 'P1': 3.1, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 3.1, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}, 'C1': {'Qda': 2.0, 'Qd': 2.0, 'Qa': 2.0, 'Q0': 2.0, 'P5': 2.0, 'P4': 2.0, 'P3': 2.3, 'P2': 2.3, 'P1': 2.7, 'Nda': 2.7, 'Nd': 2.7, 'Na': 2.7, 'N0': 2.7, 'C5': 3.1, 'C4': 3.1, 'C3': 3.5, 'C2': 3.5, 'C1': 3.5}}

class Particle(namedtuple("Particle", ["ID", "molecule_ID", "type", "x", "y", "z", "charge"])):
    def set_charge(self, charge):
        return self._replace(charge = charge)
    def set_type(self, type):
        return self._replace(type = type)
    def shift(self, dx, dy, dz):
        return self._replace(x=self.x + dx)._replace(y=self.y + dy)._replace(z=self.z + dz)
    def distance(self, i, j, k):
        return np.sqrt((self.x-i)**2 + (self.y-j)**2 + (self.z-k)**2)
    def rotate(self, yaw=0, pitch=0, roll=0, shift_x=0, shift_y=0, shift_z=0):
        Rx = np.asmatrix([[1,0,0],[0,np.cos(yaw),-np.sin(yaw)],[0,np.sin(yaw),np.cos(yaw)]])
        Ry = np.asmatrix([[np.cos(pitch),0,np.sin(pitch)],[0,1,0],[-np.sin(pitch),0,np.cos(pitch)]])
        Rz = np.asmatrix([[np.cos(roll),-np.sin(roll),0],[np.sin(roll),np.cos(roll),0],[0,0,1]])
        [_x_new], [_y_new], [_z_new] = np.array(Rz * Ry * Rx * np.asmatrix([[self.x - shift_x],[self.y - shift_y],[self.z - shift_z]]))
        return self._replace(x=_x_new + shift_x)._replace(y=_y_new + shift_y)._replace(z=_z_new + shift_z)
    def __str__(self):
        return ' '.join([str(self.ID), str(self.molecule_ID), str(self.type), str(self.charge * const_e), str(self.x * reciprocal_const_Rc_A), str(self.y * reciprocal_const_Rc_A), str(self.z * reciprocal_const_Rc_A)])
    def output_real_unit(self):
        return ' '.join([str(self.ID), str(self.molecule_ID), str(self.type), str(self.charge), str(self.x), str(self.y), str(self.z)])
class Bond(namedtuple("Bond", ["ID", "type", "p1", "p2"])):
    def __str__(self):
        return ' '.join([str(self.ID), str(self.type), str(self.p1.ID), str(self.p2.ID)])
class Angle(namedtuple("Angle", ["ID", "type", "p1", "p2", "p3"])):
    def __str__(self):
        return ' '.join([str(self.ID), str(self.type), str(self.p1.ID), str(self.p2.ID), str(self.p3.ID)])
class Dihedral(namedtuple("Dihedral", ["ID", "type", "p1", "p2", "p3", "p4"])):
    def __str__(self):
        return ' '.join([str(self.ID), str(self.type), str(self.p1.ID), str(self.p2.ID), str(self.p3.ID), str(self.p4.ID)])
class Improper(namedtuple("Improper", ["ID", "type", "p1", "p2", "p3", "p4"])):
    def __str__(self):
        return ' '.join([str(self.ID), str(self.type), str(self.p1.ID), str(self.p2.ID), str(self.p3.ID), str(self.p4.ID)])
class Types:
    class Atom_Type:
        def __init__(self, name, ID, **args):
            self.name = name
            self.ID = ID
            self.index = args.get('index', name)
            if self.index[0] == "S":
                self.index = self.index[1:]
            if self.index == "AC1":
                self.index = "C1"
            if self.index == "AC2":
                self.index = "C2"
            self.size = args.get('size', 3 if name[0] == 'S' else 4)
            self.mass = args.get('mass', (45 / 72) if name[0] == 'S' else 1)
            self.solid = args.get('solid', False)
            self.virtual = args.get('virtual', False)
        def __str__(self):
            return ' '.join([str(self.ID), str(self.mass), '#', str(self.name), '\n'])
        def output_real_unit(self):
            return ' '.join([str(self.ID), str(self.mass * 72), '#', str(self.name), '\n'])

        def a_ij(self, other):
            if self.virtual or other.virtual:
                return 0
            p1 = self.index
            p2 = other.index
            if not (p1 in martini_table and p2 in martini_table):
                raise ValueError(f'undefined particle: {p1}, {p2}')
            if p2 in martini_table[p1]:
                chi = (martini_table[p1][p2] - 0.5 * (martini_table[p1][p1] + martini_table[p2][p2])) * (-CONST_Z / RTNA_KJpermol)
            elif p1 in martini_table[p2]:
                chi = (martini_table[p2][p1] - 0.5 * (martini_table[p2][p2] + martini_table[p1][p1])) * (-CONST_Z / RTNA_KJpermol)
            else:
                raise ValueError(f'undefined particle: {p1}, {p2}')
            return (a_ww + chi * 1 / delta_aij_lambda) * self.Rc(other)**2 

        def gamma(self, other):
            if self.virtual or other.virtual:
                return 0
            elif self.size == other.size and self.size in gamma_dict:
                return gamma_dict[self.size]
            elif self.size != other.size and (self.size in gamma_mix_dict or other.size in gamma_mix_dict):
                return gamma_mix_dict[self.size] if self.size in gamma_mix_dict else gamma_mix_dict[other.size]
            else:
                return 4.5
            if self.size == 3 and other.size == 3:
                return 5.322375198853
            else:
                return 4.5
        def Rc(self, other):
            if (self.size == 3 and other.size == 3):
                return (3/4)**(1/3)
            elif (self.size == 3) or (other.size == 3):
                    _size1 = 4 if self.size == 3 else self.size
                    _size2 = 4 if other.size == 3 else other.size
                    return ((_size1/4)**(1/3) + (_size2/4)**(1/3)) / 2
            return ((self.size / 4)**(1/3) + (other.size / 4)**(1/3)) / 2

    class Bond_Type(namedtuple("Bond_Type", ["K", "r", "ID"])):
        def __str__(self):
            return ' '.join([str(self.ID), str(self.K * const_Rc_nm**2 / RTNA_KJpermol / 2), str(self.r * reciprocal_const_Rc_nm), '#', str(self.K), str(self.r), '\n'])        
        def output_real_unit(self):
            return ' '.join([str(self.ID), str(self.K / 10**2 * 0.239006 / 2), str(self.r * 10), '#', str(self.K), str(self.r), '\n'])
        def __hash__(self):
            return hash((self.K, self.r))
    class Angle_Type(namedtuple("Angle_Type", ["K", "theta", "ID"])):
        def __str__(self):
            #angle_style cosine/squared
            return ' '.join([str(self.ID), str(self.K / RTNA_KJpermol / 2), str(self.theta), '#', str(self.K), '\n'])        
        def output_real_unit(self):
            #angle_style cosine/squared
            return ' '.join([str(self.ID), str(self.K * 0.239006 / 2), str(self.theta), '#', str(self.K), '\n'])
        def __hash__(self):
            return hash((self.K, self.theta))
    class Dihedral_Type(namedtuple("Dihedral_Type", ["K", "theta", "n", "ID"])):
        def __str__(self):
            #dihedral_style charmm
            return ' '.join([str(self.ID), str(self.K / RTNA_KJpermol), str(self.n), str(self.theta), '0', '#', str(self.K), '\n'])        
        def output_real_unit(self):
            #dihedral_style charmm
            return ' '.join([str(self.ID), str(self.K * 0.239006), str(self.n), str(self.theta), '0', '#', str(self.K), '\n'])
        def __hash__(self):
            return hash((self.K, self.theta, self.n))
    class Improper_Type(namedtuple("Improper_Type", ["K", "theta", "ID"])):
        def __str__(self):
            return ' '.join([str(self.ID), str(self.K / RTNA_KJpermol / 2), str(self.theta), '#', str(self.K), '\n'])        
        def output_real_unit(self):
            return ' '.join([str(self.ID), str(self.K * 0.239006 / 2), str(self.theta), '#', str(self.K), '\n'])
        def __hash__(self):
            return hash((self.K, self.theta))
    def __init__(self):
        self.atoms = {}
        self.beads = []
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.impropers = {}

    def pairs(self):
        result = ""
        for i in range(1, len(self.beads) + 1):
            for j in range(i, len(self.beads) + 1):
                p1 = self.beads[i - 1]
                p2 = self.beads[j - 1]
                if (p1.solid or p2.solid) and (not (p1.virtual or p2.virtual)):
                    result += f"pair_coeff {i} {j} dpd {p1.a_ij(p2) + 10000} {0} {p1.Rc(p2)} # {p1.name} {p2.name}\n"
                else:
                    result += f"pair_coeff {i} {j} dpd {p1.a_ij(p2)} {p1.gamma(p2)} {p1.Rc(p2)} # {p1.name} {p2.name}\n"
        return result
    def pairs_MARTINI(self):
        result = ""
        for i in range(1, len(self.beads) + 1):
            for j in range(i, len(self.beads) + 1):
                p1 = self.beads[i - 1]
                p2 = self.beads[j - 1]
                if ((p1.index in ["Qa", "Qd", "Qda", "Q0"]) and (p2.index in ["C1", "C2"])) or ((p2.index in ["Qa", "Qd", "Qda", "Q0"]) and (p1.index in ["C1", "C2"])):
                    _sigma = 6.2
                else:
                    _sigma = 4.7
                _epsilon = martini_table[p1.index][p2.index] * 0.239006
                if (p1.size == 3 and p2.size == 3):
                    _sigma = 4.3
                    _epsilon = _epsilon * 0.75
                result += f"pair_coeff {i} {j} {_epsilon} {_sigma} # {p1.name} {p2.name}\n"
        return result

    def atom_id(self, name, **args):
        if name in self.atoms:
            return self.atoms[name].ID
        else:
            id = len(self.atoms) + 1
            self.atoms[name] = Types.Atom_Type(name=name, ID=id, **args)
            self.beads.append(self.atoms[name])
            assert self.beads[id - 1] == self.atoms[name]
            return id
    def bond_id(self, K, r):
        key = (K, r)
        if key in self.bonds:
            return self.bonds[key].ID
        else:
            id = len(self.bonds) + 1
            self.bonds[key] = Types.Bond_Type(K, r, id)
            return id
    def angle_id(self, K, theta):
        key = (K, theta)
        if key in self.angles:
            return self.angles[key].ID
        else:
            id = len(self.angles) + 1
            self.angles[key] = Types.Angle_Type(K, theta, id)
            return id
    def dihedral_id(self, K, theta, n):
        key = (K, theta, n)
        if key in self.dihedrals:
            return self.dihedrals[key].ID
        else:
            id = len(self.dihedrals) + 1
            self.dihedrals[key] = Types.Dihedral_Type(K, theta, n, id)
            return id
    def improper_id(self, K, theta):
        key = (K, theta)
        if key in self.impropers:
            return self.impropers[key].ID
        else:
            id = len(self.impropers) + 1
            self.impropers[key] = Types.Improper_Type(K, theta, id)
            return id

def Phrase_PDB_line(l):
    record_name = l[:6]
    if record_name != "ATOM  ":
        return None
    atom_name = l[12:16].replace(' ', '')
    resname = l[17:21].replace(' ', '')
    resid = int(l[22:27].replace(' ', ''))
    x = float(l[30:38].replace(' ', ''))
    y = float(l[38:46].replace(' ', ''))
    z = float(l[46:54].replace(' ', ''))
    if resname == "HSP" or resname == "HSD":
        resname = "HIS"
    return (atom_name, resname, resid, x, y, z)
def martini_top_split(l):
    if ';' in l:
        l = l[:l.find(';')]
    if '#' in l:
        l = l[:l.find('#')]
    result = []
    items = l.split()
    for i in items:
        try:
            int(i)
            result.append(int(i))
            continue
        except ValueError:
            pass
        try:
            float(i)
            result.append(float(i))
            continue
        except ValueError:
            pass
        result.append(i)
    return result
class Data:
    def __init__(self, x=None, y=None, z=None):
        self.molecules = defaultdict(list)
        self.groups = defaultdict(list)
        self.atoms = {}
        self.residues = {}
        self.bonds = {}
        self.angles = {}
        self.dihedrals = {}
        self.impropers = {}
        self.types = Types()

        self.backbone = []

        self.x_len = x
        self.y_len = y
        self.z_len = z
    def add_particle(self, molecule_name: str, type_name: str, x: float, y: float, z: float, charge: float = 0.0, group_name: str = None, residue_name: str = None, residue_id: int=None, bead_reference_pdb: str = None):
        type_id = self.types.atom_id(type_name)
        if molecule_name in self.molecules:
            molecule_id = self.atoms[self.molecules[molecule_name][0]].molecule_ID
        else:
            molecule_id = len(self.molecules) + 1
        particle = Particle(len(self.atoms) + 1, molecule_id, type_id, x, y, z, charge)
        self.atoms[particle.ID] = particle
        if bead_reference_pdb != "BB":
            if (residue_name, residue_id) not in self.residues:
                self.residues[(residue_name, residue_id)] = [particle]
            else:
                self.residues[(residue_name, residue_id)].append(particle)
        self.molecules[molecule_name].append(particle.ID)
        if group_name is not None:
            if isinstance(group_name, str):
                self.groups[group_name].append(particle.ID)
            elif isinstance(group_name, list):
                for g in group_name:
                    self.groups[g].append(particle.ID)
        return particle.ID
    def add_bond(self, p1, p2, K, r):
        if isinstance(p1, int):
            p1 = self.atoms[p1]
        if isinstance(p2, int):
            p2 = self.atoms[p2]
        bond_type_id = self.types.bond_id(K, r)
        current = Bond(len(self.bonds) + 1, bond_type_id, p1, p2)
        self.bonds[current.ID] = current
        return current.ID
    def add_angle(self, p1, p2, p3, K, theta):
        if isinstance(p1, int):
            p1 = self.atoms[p1]
        if isinstance(p2, int):
            p2 = self.atoms[p2]
        if isinstance(p3, int):
            p3 = self.atoms[p3]
        angle_type_id = self.types.angle_id(K, theta)
        current = Angle(len(self.angles) + 1, angle_type_id, p1, p2, p3)
        self.angles[current.ID] = current
    def add_dihedral(self, p1, p2, p3, p4, K, theta, n):
        if isinstance(p1, int):
            p1 = self.atoms[p1]
        if isinstance(p2, int):
            p2 = self.atoms[p2]
        if isinstance(p3, int):
            p3 = self.atoms[p3]
        if isinstance(p4, int):
            p4 = self.atoms[p4]
        dihedral_type_id = self.types.dihedral_id(K, theta, n)
        current = Dihedral(len(self.dihedrals) + 1, dihedral_type_id, p1, p2, p3, p4)
        self.dihedrals[current.ID] = current
    def add_improper(self, p1, p2, p3, p4, K, theta):
        if isinstance(p1, int):
            p1 = self.atoms[p1]
        if isinstance(p2, int):
            p2 = self.atoms[p2]
        if isinstance(p3, int):
            p3 = self.atoms[p3]
        if isinstance(p4, int):
            p4 = self.atoms[p4]
        improper_type_id = self.types.improper_id(K, theta)
        current = Improper(len(self.impropers) + 1, improper_type_id, p1, p2, p3, p4)
        self.impropers[current.ID] = current
    def center_box(self):
        x_max_box = float(self.x_max)
        x_min_box = float(self.x_min)
        y_max_box = float(self.y_max)
        y_min_box = float(self.y_min)
        z_max_box = float(self.z_max)
        z_min_box = float(self.z_min)
        dx = (x_max_box + x_min_box) / 2
        dy = (y_max_box + y_min_box) / 2
        dz = (z_max_box + z_min_box) / 2
        for i in self.atoms:
            self.atoms[i] = self.atoms[i].shift(-dx, -dy, -dz)
    def shift_groups(self, group_names, dx, dy, dz):
        #assert mol_name in self.molecules
        _to_shift = set()
        for g in group_names:
            for i in self.groups[g]:
                _to_shift.add(i)
        for i in _to_shift:
            self.atoms[i] = self.atoms[i].shift(dx, dy, dz)
    def rotate_groups(self, group_names, yaw=0, pitch=0, roll=0):
        _to_rotate = set()
        _x_list = []
        _y_list = []
        _z_list = []
        for g in group_names:
            for i in self.groups[g]:
                _to_rotate.add(i)
                _x_list.append(self.atoms[i].x)
                _y_list.append(self.atoms[i].y)
                _z_list.append(self.atoms[i].z)
        for i in _to_rotate:
            self.atoms[i] = self.atoms[i].rotate(yaw, pitch, roll, np.mean(_x_list), np.mean(_y_list), np.mean(_z_list))
    def center_groups(self, group_names):
        _to_shift = set()
        for g in group_names:
            for i in self.groups[g]:
                _to_shift.add(i)
        x_list = []
        y_list = []
        z_list = []
        for i in _to_shift:
            x_list.append(self.atoms[i].x)
            y_list.append(self.atoms[i].y)
            z_list.append(self.atoms[i].z)
        x_center = (max(x_list) + min(x_list)) / 2
        y_center = (max(y_list) + min(y_list)) / 2
        z_center = (max(z_list) + min(z_list)) / 2
        self.shift_groups(group_names, -x_center, -y_center, -z_center)
    def center_molecule(self, mol_name, shift_x = False, shift_y = False, shift_z = False):
        assert mol_name in self.molecules
        x_list = []
        y_list = []
        z_list = []
        for atom in self.molecules[mol_name]:
            x_list.append(atom.x)
            y_list.append(atom.y)
            z_list.append(atom.z)
        x_center = (max(x_list) + min(x_list)) / 2
        y_center = (max(y_list) + min(y_list)) / 2
        z_center = (max(z_list) + min(z_list)) / 2
        x_max_box = float(self.x_max)
        x_min_box = float(self.x_min)
        y_max_box = float(self.y_max)
        y_min_box = float(self.y_min)
        z_max_box = float(self.z_max)
        z_min_box = float(self.z_min)

        dx = -((x_max_box + x_min_box) / 2 - x_center)
        dy = -((y_max_box + y_min_box) / 2 - y_center)
        dz = -((z_max_box + z_min_box) / 2 - z_center)

        for atom in self.atoms:
            if shift_x:
                if x_min_box <= atom.x - dx <= x_max_box:
                    atom.x -= dx
                elif x_min_box > atom.x - dx:
                    atom.x = atom.x - dx + (x_max_box - x_min_box)
                elif x_max_box < atom.x - dx:
                    atom.x = atom.x - dx - (x_max_box - x_min_box)
            if shift_y:
                if y_min_box <= atom.y - dy <= y_max_box:
                    atom.y -= dy
                elif y_min_box > atom.y - dy:
                    atom.y = atom.y - dy + (y_max_box - y_min_box)
                elif y_max_box < atom.y - dy:
                    atom.y = atom.y - dy - (y_max_box - y_min_box)
            if shift_z:
                if z_min_box <= atom.z - dz <= z_max_box:
                    atom.z -= dz
                elif z_min_box > atom.z - dz:
                    atom.z = atom.z - dz + (z_max_box - z_min_box)
                elif z_max_box < atom.z - dz:
                    atom.z = atom.z - dz - (z_max_box - z_min_box)
    def volume(self, solid_volume = 0):
        total_size = 0
        for i in self.atoms:
            if self.types.beads[self.atoms[i].type - 1].solid == False:
                total_size += self.types.beads[self.atoms[i].type - 1].size
        return total_size * 30 + solid_volume
    @property
    def x_max(self):
        from math import ceil
        x_list = []
        for i in self.atoms:
            x_list.append(self.atoms[i].x)
        return ceil(max(x_list) * 100 * reciprocal_const_Rc_A) / 100
    @property
    def x_max_real(self):
        from math import ceil
        x_list = []
        for i in self.atoms:
            x_list.append(self.atoms[i].x)
        return ceil(max(x_list) * 100) / 100
    
    @property
    def x_min(self):
        from math import floor
        x_list = []
        for i in self.atoms:
            x_list.append(self.atoms[i].x)
        return floor(min(x_list) * 100 * reciprocal_const_Rc_A) / 100    
    @property
    def x_min_real(self):
        from math import floor
        x_list = []
        for i in self.atoms:
            x_list.append(self.atoms[i].x)
        return floor(min(x_list) * 100) / 100
    @property
    def y_max(self):
        from math import ceil
        y_list = []
        for i in self.atoms:
            y_list.append(self.atoms[i].y)
        return ceil(max(y_list) * 100 * reciprocal_const_Rc_A) / 100    
    @property
    def y_max_real(self):
        from math import ceil
        y_list = []
        for i in self.atoms:
            y_list.append(self.atoms[i].y)
        return ceil(max(y_list) * 100) / 100
    @property
    def y_min(self):
        from math import floor
        y_list = []
        for i in self.atoms:
            y_list.append(self.atoms[i].y)
        return floor(min(y_list) * 100 * reciprocal_const_Rc_A) / 100    
    @property
    def y_min_real(self):
        from math import floor
        y_list = []
        for i in self.atoms:
            y_list.append(self.atoms[i].y)
        return floor(min(y_list) * 100) / 100
    @property
    def z_max(self):
        from math import ceil
        z_list = []
        for i in self.atoms:
            z_list.append(self.atoms[i].z)
        return ceil(max(z_list) * 100 * reciprocal_const_Rc_A) / 100    
    @property
    def z_max_real(self):
        from math import ceil
        z_list = []
        for i in self.atoms:
            z_list.append(self.atoms[i].z)
        return ceil(max(z_list) * 100) / 100
    @property
    def z_min(self):
        from math import floor
        z_list = []
        for i in self.atoms:
            z_list.append(self.atoms[i].z)
        return floor(min(z_list) * 100 * reciprocal_const_Rc_A) / 100    
    @property
    def z_min_real(self):
        from math import floor
        z_list = []
        for i in self.atoms:
            z_list.append(self.atoms[i].z)
        return floor(min(z_list) * 100) / 100
    @property
    def total_charge(self):
        _result = 0
        for i in self.atoms:
            _result += self.atoms[i].charge
        return _result
    
    @property
    def group_boxes(self):
        boxes = {}
        for mol_name in self.groups:
            x_list = [self.atoms[i].x for i in self.groups[mol_name]]
            y_list = [self.atoms[i].y for i in self.groups[mol_name]]
            z_list = [self.atoms[i].z for i in self.groups[mol_name]]
            boxes[mol_name] = ((min(x_list, default=0), max(x_list, default=0)), (min(y_list, default=0), max(y_list, default=0)), (min(z_list, default=0), max(z_list, default=0)))
        return boxes
            
    def fill_water(self, ion_concentration=0.15, solid_volume = 0.0, factor = 1.0, factor_rho = 1.0, skip_box = True):
        #ion_concentration unit: mol/L
        N_ion_pairs = 6.0221408 * 10**23 * 10**(-27) * ion_concentration * self.x_len * self.y_len * self.z_len
        
        self.types.atom_id(name='W', index="P4", size=4 * factor, mass=1 * factor)
        W_id = self.types.atom_id('W')
        self.types.atom_id(name='K', index="Qd", size=4, mass=1)
        K_id = self.types.atom_id("K")
        self.types.atom_id(name='CL', index="Qa", size=4, mass=1)
        CL_id = self.types.atom_id("CL")
        boxes = self.group_boxes

        
        K_counter = 1
        CL_counter = 1
        W_counter = 1
        remaining_charge = int(self.total_charge)
        print("remaining charge: ", remaining_charge)
        #if abs(remaining_charge) > N_ion_pairs:
        #    raise("Ion concentration is not sufficient to neutralize the system")
        while remaining_charge > 0.001:
            skip = False
            x = np.random.uniform(-self.x_len / 2, self.x_len / 2)
            y = np.random.uniform(-self.y_len / 2, self.y_len / 2)
            z = np.random.uniform(-self.z_len / 2, self.z_len / 2)
            if skip_box == True:
                for mol_name in boxes:
                    ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = boxes[mol_name]
                    if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
                        skip = True
                        break
            if not skip:
                self.add_particle(f"CL{CL_counter}", "CL", x, y, z, -1, "solution")
                CL_counter += 1
                remaining_charge -= 1
                N_ion_pairs -= 1
        while remaining_charge < -0.001:
            skip = False
            x = np.random.uniform(-self.x_len / 2, self.x_len / 2)
            y = np.random.uniform(-self.y_len / 2, self.y_len / 2)
            z = np.random.uniform(-self.z_len / 2, self.z_len / 2)
            if skip_box == True:
                for mol_name in boxes:
                    ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = boxes[mol_name]
                    if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
                        skip = True
                        break
            if not skip:
                self.add_particle(f"K{K_counter}", "K", x, y, z, 1, "solution")
                K_counter += 1
                remaining_charge += 1
                N_ion_pairs -= 1
        print("#ion pairs: ", N_ion_pairs)

        for i in range(int(N_ion_pairs)):
            skip = True
            while skip:
                x = np.random.uniform(-self.x_len / 2, self.x_len / 2)
                y = np.random.uniform(-self.y_len / 2, self.y_len / 2)
                z = np.random.uniform(-self.z_len / 2, self.z_len / 2)
                if skip_box == True:
                    for mol_name in boxes:
                        ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = boxes[mol_name]
                        if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
                            skip = True
                        else:
                            skip = False
                else:
                    skip = False
            self.add_particle(f"K{K_counter}", "K", x, y, z, 1, "solution")
            K_counter += 1
        for i in range(int(N_ion_pairs)):
            skip = True
            while skip:
                x = np.random.uniform(-self.x_len / 2, self.x_len / 2)
                y = np.random.uniform(-self.y_len / 2, self.y_len / 2)
                z = np.random.uniform(-self.z_len / 2, self.z_len / 2)
                if skip_box == True:
                    for mol_name in boxes:
                        ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = boxes[mol_name]
                        if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
                            skip = True
                        else:
                            skip = False
                else:
                    skip = False
            self.add_particle(f"CL{CL_counter}", "CL", x, y, z, -1, "solution")
            CL_counter += 1
        sc_grid_dis_water = (1/3 * factor_rho)**(1/3) * 2 * const_Rc_A
        x_grid = np.arange(-self.x_len/2, self.x_len/2, sc_grid_dis_water)
        y_grid = np.arange(-self.y_len/2, self.y_len/2, sc_grid_dis_water)
        z_grid = np.arange(-self.z_len/2, self.z_len/2, sc_grid_dis_water)
        for _x in x_grid:
            for _y in y_grid:
                for _z in z_grid:
                    skip = False
                    if skip_box == True:
                        for mol_name in boxes:
                            if skip:
                                break
                            ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = boxes[mol_name]
                            #if xmin - sc_grid_dis_water <= _x <= xmax + sc_grid_dis_water and ymin - sc_grid_dis_water <= _y <= ymax + sc_grid_dis_water and zmin - sc_grid_dis_water <= _z <= zmax + sc_grid_dis_water:
                            if xmin <= _x <= xmax and ymin <= _y <= ymax and zmin <= _z <= zmax:
                                skip = True
                                break
                    if not skip:
                        self.add_particle("WATER" + str(W_counter), "W", _x, _y, _z, 0, "solution")
                        W_counter += 1
        remaining_volume = self.x_len * self.y_len * self.z_len - self.volume(solid_volume)
        print("remaining volume: ", remaining_volume)
        if remaining_volume < 0:
            raise("fill water remaining volume is negative!")
        while remaining_volume > 0:
            skip = False
            x = np.random.uniform(-self.x_len / 2, self.x_len / 2)
            y = np.random.uniform(-self.y_len / 2, self.y_len / 2)
            z = np.random.uniform(-self.z_len / 2, self.z_len / 2)
            if skip_box == True:
                for mol_name in boxes:
                    ((xmin, xmax), (ymin, ymax), (zmin, zmax)) = boxes[mol_name]
                    if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
                        skip = True
                        break
            if not skip:
                self.add_particle("WATER" + str(W_counter), "W", x, y, z, 0, "solution")
                W_counter += 1
                remaining_volume -= self.types.beads[W_id - 1].size * 30

    def adjust_ph(self, ph=7):
        for (residue_name, residue_id) in self.residues:
            #if residue_name in ["ASP", "GLU", "LYS", "ARG", "HIS"]:
            if ph <= 3.9 and residue_name == "ASP":
                #change ASP to neutral
                particle = self.residues[(residue_name, residue_id)][0]
                assert len(self.residues[(residue_name, residue_id)]) == 1
                assert particle.type == self.types.atom_id("Qa")
                new_id = self.types.atom_id("P3")
                self.atoms[particle.ID] = self.atoms[particle.ID].set_charge(0).set_type(new_id)
            if ph <= 4.2 and residue_name == "GLU":
                #change GLU to neutral
                particle = self.residues[(residue_name, residue_id)][0]
                assert len(self.residues[(residue_name, residue_id)]) == 1
                assert particle.type == self.types.atom_id("Qa")
                new_id = self.types.atom_id("P1")
                self.atoms[particle.ID] = self.atoms[particle.ID].set_charge(0).set_type(new_id)
            if ph >= 10.5 and residue_name == "LYS":
                #change LYS to neutral
                assert len(self.residues[(residue_name, residue_id)]) == 2
                if self.residues[(residue_name, residue_id)][0].type == self.types.atom_id("Qd"):
                    assert self.residues[(residue_name, residue_id)][1].type == self.types.atom_id("C3")
                    particle = self.residues[(residue_name, residue_id)][0]
                else:
                    assert self.residues[(residue_name, residue_id)][0].type == self.types.atom_id("C3")
                    particle = self.residues[(residue_name, residue_id)][1]
                new_id = self.types.atom_id("P1")
                self.atoms[particle.ID] = self.atoms[particle.ID].set_charge(0).set_type(new_id)
            if ph >= 12.5 and residue_name == "ARG":
                #change ARG to neutral
                assert len(self.residues[(residue_name, residue_id)]) == 2
                if self.residues[(residue_name, residue_id)][0].type == self.types.atom_id("Qd"):
                    assert self.residues[(residue_name, residue_id)][1].type == self.types.atom_id("N0")
                    particle = self.residues[(residue_name, residue_id)][0]
                else:
                    assert self.residues[(residue_name, residue_id)][0].type == self.types.atom_id("N0")
                    particle = self.residues[(residue_name, residue_id)][1]
                new_id = self.types.atom_id("P4")
                self.atoms[particle.ID] = self.atoms[particle.ID].set_charge(0).set_type(new_id)
            if ph <= 6 and residue_name == "HIS":
                #change HIS to charged
                assert len(self.residues[(residue_name, residue_id)]) == 3
                if self.residues[(residue_name, residue_id)][0].type == self.types.atom_id("SP1"):
                    assert (self.residues[(residue_name, residue_id)][1].type == self.types.atom_id("SP1") and self.residues[(residue_name, residue_id)][2].type == self.types.atom_id("SC4")) or (self.residues[(residue_name, residue_id)][1].type == self.types.atom_id("SC4") and self.residues[(residue_name, residue_id)][2].type == self.types.atom_id("SP1"))
                    particle = self.residues[(residue_name, residue_id)][0]
                elif self.residues[(residue_name, residue_id)][0].type == self.types.atom_id("SC4"):
                    assert self.residues[(residue_name, residue_id)][1].type == self.types.atom_id("SP1") and self.residues[(residue_name, residue_id)][2].type == self.types.atom_id("SP1")
                    particle = self.residues[(residue_name, residue_id)][1]
                new_id = self.types.atom_id("SQd")
                self.atoms[particle.ID] = self.atoms[particle.ID].set_charge(1).set_type(new_id)
    def add_molecule(self, pdb_file, top_file, mol_name, group_name=None):
        def get_section(topfile, keyword):
            pos = topfile.find(keyword) + len(keyword)
            content = topfile[pos:topfile.find("[", pos,)].splitlines()
            content = list(filter(None, content))
            result = []
            for l in content:
                if not l[0] == ';':
                    result.append(l)
            return result
        protein_data = []
        for l in pdb_file:
            if Phrase_PDB_line(l) == None:
                continue
            bead_reference_pdb, resname_pdb, resid_pdb, x_pdb, y_pdb, z_pdb = Phrase_PDB_line(l)
            if resname_pdb in ["W", 'NA', "CL"]:
                continue
            protein_data.append([bead_reference_pdb, resname_pdb, resid_pdb, x_pdb, y_pdb, z_pdb])
        top_reference = {}
        if "[ atoms ]" in top_file:
            section_atom = get_section(top_file, "[ atoms ]")
        elif "[atoms]" in top_file:
            section_atom = get_section(top_file, "[atoms]")
        else:
            section_atom = ""
        if "[ bonds ]" in top_file:
            section_bond = get_section(top_file, "[ bonds ]")
        elif "[bonds]" in top_file:
            section_bond = get_section(top_file, "[bonds]")
        else:
            section_bond = ""
        if "[ constraints ]" in top_file:
            section_constraints = get_section(top_file, "[ constraints ]")
        elif "[constraints]" in top_file:
            section_constraints = get_section(top_file, "[constraints]")
        else:
            section_constraints = ""
        if "[ angles ]" in top_file:
            section_angles = get_section(top_file, "[ angles ]")
        elif "[angles]" in top_file:
            section_angles = get_section(top_file, "[angles]")
        else:
            section_angles = ""
        if "[ dihedrals ]" in top_file:
            section_dihedrals = get_section(top_file, "[ dihedrals ]")
        elif "[dihedrals]" in top_file:
            section_dihedrals = get_section(top_file, "[dihedrals]")
        else:
            section_dihedrals = ""

        assert len(protein_data) == len(section_atom)
        for i in range(len(section_atom)):
            l = section_atom[i]
            if len(martini_top_split(l)) == 7:
                tmp, atom_type_martini, resid_martini, resname_martini, bead_reference_martini, atom_id_martini, charge_matini = martini_top_split(l)
                [bead_reference_pdb, resname_pdb, resid_pdb, x_pdb, y_pdb, z_pdb] = protein_data[i]
                current_atom_ID = self.add_particle(mol_name, atom_type_martini, x_pdb, y_pdb, z_pdb, charge_matini, group_name, resname_pdb, resid_pdb, bead_reference_pdb)
                top_reference[atom_id_martini] = current_atom_ID
        for l in section_bond:
            if len(martini_top_split(l)) == 5:
                atom_1_id, atom_2_id, _tmp, r_martini, k_martini = martini_top_split(l)
                assert _tmp == 1
                self.add_bond(top_reference[atom_1_id], top_reference[atom_2_id], k_martini, r_martini)
        for l in section_constraints:
            if len(martini_top_split(l)) == 4:
                atom_1_id, atom_2_id, _tmp, r_martini = martini_top_split(l)
                assert _tmp == 1
                self.add_bond(top_reference[atom_1_id], top_reference[atom_2_id], k_constraints, r_martini)
        for l in section_angles:
            if len(martini_top_split(l)) == 6:
                atom_1_id, atom_2_id, atom_3_id, _tmp, theta_martini, k_martini = martini_top_split(l)
                assert _tmp == 2
                self.add_angle(top_reference[atom_1_id], top_reference[atom_2_id], top_reference[atom_3_id], k_martini, theta_martini)
        for l in section_dihedrals:
            if len(martini_top_split(l)) == 8:
                atom_1_id, atom_2_id, atom_3_id, atom_4_id, _tmp, theta_martini, k_martini, n_martini = martini_top_split(l)
                assert _tmp == 1
                self.add_dihedral(top_reference[atom_1_id], top_reference[atom_2_id], top_reference[atom_3_id], top_reference[atom_4_id], k_martini, theta_martini, n_martini)
            if len(martini_top_split(l)) == 7:
                atom_1_id, atom_2_id, atom_3_id, atom_4_id, _tmp, theta_martini, k_martini = martini_top_split(l)
                assert _tmp == 2
                self.add_improper(top_reference[atom_1_id], top_reference[atom_2_id], top_reference[atom_3_id], top_reference[atom_4_id], k_martini, theta_martini)
        
    def add_water_from_pdb(self, water_pdb):
        import random
        import numpy as np
        self.types.atom_id(name='W', index="P4", size=4, mass=1)
        water_type_id = self.types.atom_id('W')
        self.types.atom_id(name='K', index="Qd", size=4, mass=1)
        Positive_ion_type_id = self.types.atom_id("K")
        self.types.atom_id(name='CL', index="Qa", size=4, mass=1)
        Negative_ion_type_id = self.types.atom_id("CL")
        _last_water_mol_id = None
        water_xyz = []
        POSITIVEION_xyz = []
        NEGATIVEION_xyz = []
        for l in water_pdb:
            if Phrase_PDB_line(l) == None:
                continue
            bead_reference, resname, resid, x, y, z = Phrase_PDB_line(l)
            if resname == "W":
                water_xyz.append([x, y, z])
            elif resname == "NA":
                POSITIVEION_xyz.append([x, y, z])
            elif resname == "CL":
                NEGATIVEION_xyz.append([x, y, z])
            else:
                continue
        while int(int(self.total_charge) + len(POSITIVEION_xyz) - len(NEGATIVEION_xyz)) > 0:
            NEGATIVEION_xyz.append(water_xyz[-1])
            del water_xyz[-1]
        while int(int(self.total_charge) + len(POSITIVEION_xyz) - len(NEGATIVEION_xyz)) < 0:
            POSITIVEION_xyz.append(water_xyz[-1])
            del water_xyz[-1]
        counter = 1
        for xyz in water_xyz:
            x, y, z = xyz
            self.add_particle('WATER' + str(counter), "W", x, y, z, 0, "solution")
            counter += 1
        for xyz in POSITIVEION_xyz:
            x, y, z = xyz
            self.add_particle('WATER' + str(counter), "Qd", x, y, z, 1, "solution")
            counter += 1
        for xyz in NEGATIVEION_xyz:
            x, y, z = xyz
            self.add_particle('WATER' + str(counter), "Qa", x, y, z, -1, "solution")
            counter += 1
        assert -0.001 < self.total_charge < 0.001
    
    def dump(self, outfile = "protein.data", additional_group_str = "", MARTINI = "False"):
        self.center_box()

        group_str = ""
        for name in self.groups:
            tmp = ' '.join(str(a) for a in self.groups[name])
            if len(self.groups[name]) > 0 and len(self.groups[name]) == max(self.groups[name]) - min(self.groups[name]) + 1:
                tmp = str(min(self.groups[name])) + ':' + str(max(self.groups[name]))
            group_str += ' '.join(["group", name, 'id', tmp, '\n'])
        with open("group.settings", 'w') as group_out:
            group_out.write(group_str + additional_group_str)
        if not MARTINI:
            with open(outfile, 'w') as output:
                output.write('\n')
                if len(self.atoms) != 0:
                    output.write(str(len(self.atoms)) + ' atoms\n')
                if len(self.bonds) != 0:
                    output.write(str(len(self.bonds)) + ' bonds\n')
                if len(self.angles) != 0:
                    output.write(str(len(self.angles)) + ' angles\n')
                if len(self.dihedrals) != 0:
                    output.write(str(len(self.dihedrals)) + ' dihedrals\n')
                if len(self.impropers) != 0:
                    output.write(str(len(self.impropers)) + ' impropers\n')
                output.write('\n')
                if len(self.types.atoms) != 0:
                    output.write(str(len(self.types.atoms)) + ' atom types\n')
                if len(self.types.bonds) != 0:
                    output.write(str(len(self.types.bonds)) + ' bond types\n')
                if len(self.types.angles) != 0:
                    output.write(str(len(self.types.angles)) + ' angle types\n')
                if len(self.types.dihedrals) != 0:
                    output.write(str(len(self.types.dihedrals)) + ' dihedral types\n')
                if len(self.types.impropers) != 0:
                    output.write(str(len(self.types.impropers)) + ' improper types\n')
                output.write('\n')
                output.write(str(self.x_min) + ' ' + str(self.x_max) + ' xlo xhi\n')
                output.write(str(self.y_min) + ' ' + str(self.y_max) + ' ylo yhi\n')
                output.write(str(self.z_min) + ' ' + str(self.z_max) + ' zlo zhi\n')
                
                if len(self.types.atoms) > 0:
                    output.write('\nMasses\n\n')
                    for name in self.types.atoms:
                        output.write(str(self.types.atoms[name]))
                if len(self.atoms) != 0:
                    output.write('\nAtoms\n\n')
                    for i in self.atoms:
                        output.write(str(self.atoms[i]) + '\n')
                    output.write('\n')
                if len(self.bonds) != 0:
                    output.write('Bonds\n')
                    output.write('\n')
                    for b in self.bonds:
                        output.write(str(self.bonds[b]) + '\n')
                    output.write('\n')
                    output.write('Bond Coeffs\n')
                    output.write('\n')
                    for i in self.types.bonds:
                        output.write(str(self.types.bonds[i]))
                    output.write('\n')
                if len(self.angles) != 0:
                    output.write('Angles\n')
                    output.write('\n')
                    for i in self.angles:
                        output.write(str(self.angles[i]) + '\n')
                    output.write('\n')
                    output.write('Angle Coeffs\n')
                    output.write('\n')
                    for i in self.types.angles:
                        output.write(str(self.types.angles[i]))
                    output.write('\n')
                if len(self.dihedrals) != 0:
                    output.write('Dihedrals\n')
                    output.write('\n')
                    for i in self.dihedrals:
                        output.write(str(self.dihedrals[i]) + '\n')
                    output.write('\n')
                    output.write('Dihedral Coeffs\n')
                    output.write('\n')
                    for i in self.types.dihedrals:
                        output.write(str(self.types.dihedrals[i]))
                    output.write('\n')
                if len(self.impropers) != 0:
                    output.write('Impropers\n')
                    output.write('\n')
                    for i in self.impropers:
                        output.write(str(self.impropers[i]) + '\n')
                    output.write('\n')
                    output.write('Improper Coeffs\n')
                    output.write('\n')
                    for i in self.types.impropers:
                        output.write(str(self.types.impropers[i]))
                    output.write('\n')   
            with open("pair.settings", 'w') as outpair:
                if len(self.types.atoms) > 0:
                    #output.write('\nPairIJ Coeffs\n\n')
                    outpair.write(str(self.types.pairs()))
                    outpair.write(f"pair_coeff *{len(self.types.atoms)} *{len(self.types.atoms)} coul/slater/long\n")
                    outpair.write("dielectric 78.3\n")     
        else:
            with open(outfile, 'w') as output:
                output.write('\n')
                if len(self.atoms) != 0:
                    output.write(str(len(self.atoms)) + ' atoms\n')
                if len(self.bonds) != 0:
                    output.write(str(len(self.bonds)) + ' bonds\n')
                if len(self.angles) != 0:
                    output.write(str(len(self.angles)) + ' angles\n')
                if len(self.dihedrals) != 0:
                    output.write(str(len(self.dihedrals)) + ' dihedrals\n')
                if len(self.impropers) != 0:
                    output.write(str(len(self.impropers)) + ' impropers\n')
                output.write('\n')
                if len(self.types.atoms) != 0:
                    output.write(str(len(self.types.atoms)) + ' atom types\n')
                if len(self.types.bonds) != 0:
                    output.write(str(len(self.types.bonds)) + ' bond types\n')
                if len(self.types.angles) != 0:
                    output.write(str(len(self.types.angles)) + ' angle types\n')
                if len(self.types.dihedrals) != 0:
                    output.write(str(len(self.types.dihedrals)) + ' dihedral types\n')
                if len(self.types.impropers) != 0:
                    output.write(str(len(self.types.impropers)) + ' improper types\n')
                output.write('\n')
                output.write(str(self.x_min_real) + ' ' + str(self.x_max_real) + ' xlo xhi\n')
                output.write(str(self.y_min_real) + ' ' + str(self.y_max_real) + ' ylo yhi\n')
                output.write(str(self.z_min_real) + ' ' + str(self.z_max_real) + ' zlo zhi\n')
                
                if len(self.types.atoms) > 0:
                    output.write('\nMasses\n\n')
                    for name in self.types.atoms:
                        output.write(str(self.types.atoms[name].output_real_unit()))
                if len(self.atoms) != 0:
                    output.write('\nAtoms\n\n')
                    for i in self.atoms:
                        output.write(str(self.atoms[i].output_real_unit()) + '\n')
                    output.write('\n')
                if len(self.bonds) != 0:
                    output.write('Bonds\n')
                    output.write('\n')
                    for b in self.bonds:
                        output.write(str(self.bonds[b]) + '\n')
                    output.write('\n')
                    output.write('Bond Coeffs\n')
                    output.write('\n')
                    for i in self.types.bonds:
                        output.write(str(self.types.bonds[i].output_real_unit()))
                    output.write('\n')
                if len(self.angles) != 0:
                    output.write('Angles\n')
                    output.write('\n')
                    for i in self.angles:
                        output.write(str(self.angles[i]) + '\n')
                    output.write('\n')
                    output.write('Angle Coeffs\n')
                    output.write('\n')
                    for i in self.types.angles:
                        output.write(str(self.types.angles[i].output_real_unit()))
                    output.write('\n')
                if len(self.dihedrals) != 0:
                    output.write('Dihedrals\n')
                    output.write('\n')
                    for i in self.dihedrals:
                        output.write(str(self.dihedrals[i]) + '\n')
                    output.write('\n')
                    output.write('Dihedral Coeffs\n')
                    output.write('\n')
                    for i in self.types.dihedrals:
                        output.write(str(self.types.dihedrals[i].output_real_unit()))
                    output.write('\n')
                if len(self.impropers) != 0:
                    output.write('Impropers\n')
                    output.write('\n')
                    for i in self.impropers:
                        output.write(str(self.impropers[i]) + '\n')
                    output.write('\n')
                    output.write('Improper Coeffs\n')
                    output.write('\n')
                    for i in self.types.impropers:
                        output.write(str(self.types.impropers[i].output_real_unit()))
                    output.write('\n')
            with open("pair.MARTINI.settings", 'w') as outpair:
                if len(self.types.atoms) > 0:
                    outpair.write(str(self.types.pairs_MARTINI()))
                    outpair.write("dielectric 15\n")
                
                