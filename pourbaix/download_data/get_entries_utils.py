import qmpy_rester as qr
import numpy as np
from pymatgen.core.structure import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor

adaptor = AseAtomsAdaptor()
from ase.db import connect


def get_all_entries(maximum_pages=np.inf, page_limit=100):
    """
    Queries all entries from OQMD that meet the following criteria:
        - Formation energy <= 0
        - Energy above hull <= 0.1 eV/atom above hull
        - Is unary, binary, or ternary

    Args:
        max_pages: The max number of pages to query from the OQMD RESTful API
        page_limit: The number of entries that will be returned with each query
    """
    continue_bool = True
    batch_num = 0
    cumulative_data = []
    while continue_bool:
        with qr.QMPYRester() as q:
            kwargs = {
                "limit": page_limit,
                "offset": batch_num * page_limit,
                "filter": "stability <= 0.1 AND delta_e <= 0 AND ntypes < 4",
            }
            data = q.get_oqmd_phases(verbose=False, **kwargs)
        cumulative_data = [*cumulative_data, *data["data"]]
        batch_num += 1
        if data["links"]["next"] and batch_num <= maximum_pages:
            continue_bool = True
        else:
            continue_bool = False
    return cumulative_data


def get_icsd_entries(maximum_pages=np.inf, page_limit=100):
    """
    Queries all entries from OQMD from ICSD for data fitting purposes

    Args:
        max_pages: The max number of pages to query from the OQMD RESTful API
        page_limit: The number of entries that will be returned with each query
    """
    continue_bool = True
    batch_num = 0
    cumulative_data = []
    with qr.QMPYRester() as q:
        kwargs = {
            "limit": page_limit,
            "offset": batch_num * page_limit,
            "icsd": "T",
        }
        data = q.get_oqmd_phases(verbose=False, **kwargs)
        cumulative_data = [*cumulative_data, *data["data"]]
        batch_num += 1
        if data["links"]["next"] and batch_num <= maximum_pages:
            continue_bool = True
        else:
            continue_bool = False
    return cumulative_data


def get_ase_object(entry):
    """
    Converts the unit cell specifications queried from OQMD into an ASE atoms object

    Args:
        entry: dict-like object of an entry queried from OQMD
    """
    cell = entry["unit_cell"]
    elems, coords = [], []
    for atom_info in entry["sites"]:
        info = atom_info.split(" ")
        elems.append(info[0])
        coords.append([float(x) for x in info[2:]])
    lattice = Lattice(cell)
    s = Structure(lattice, elems, coords)
    atoms_obj = adaptor.get_atoms(s)
    entry["atoms_obj"] = atoms_obj
    entry["pmg_struct"] = s
    return atoms_obj


def write_data(entries, path):
    """
    Writes OQMD data to an ASE db.

    Args:
        entries: A list of entries queried from OQMD
        path: desired path and filename for the db file
    """
    with connect(path) as db:
        for entry in entries:
            atoms = get_ase_object(entry)
            oqmd_id = "oqmd-" + str(entry["entry_id"])
            db.write(
                atoms, source="oqmd", bulk_id=oqmd_id, formation_energy=entry["delta_e"]
            )
