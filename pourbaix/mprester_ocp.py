from pymatgen.ext.matproj import MPRester
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import IonEntry, PourbaixEntry
from pymatgen.core.ion import Ion
from pymatgen.entries.compatibility import (
    Compatibility,
    MaterialsProject2020Compatibility,
    MaterialsProjectAqueousCompatibility,
    MaterialsProjectCompatibility,
)


class MPRester_x_OCP(MPRester):
    """
    This class allows additional structures to be passed as an argumentt to MPRester.
    This
    """

    def get_pourbaix_entries_w_add_entries(
        self,
        chemsys,
        solid_compat="MaterialsProject2020Compatibility",
        additional_entries=[],
    ):
        """
        A helper function to get all entries necessary to generate
        a pourbaix diagram from the rest interface.
        Args:
            chemsys (str or [str]): Chemical system string comprising element
                symbols separated by dashes, e.g., "Li-Fe-O" or List of element
                symbols, e.g., ["Li", "Fe", "O"].
            solid_compat: Compatibility scheme used to pre-process solid DFT energies prior to applying aqueous
                energy adjustments. May be passed as a class (e.g. MaterialsProject2020Compatibility) or an instance
                (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies are used as-is.
                Default: MaterialsProject2020Compatibility
        """

        if solid_compat == "MaterialsProjectCompatibility":
            self.solid_compat = MaterialsProjectCompatibility()
        elif solid_compat == "MaterialsProject2020Compatibility":
            self.solid_compat = MaterialsProject2020Compatibility()
        elif isinstance(solid_compat, Compatibility):
            self.solid_compat = solid_compat
        else:
            raise ValueError(
                "Solid compatibility can only be 'MaterialsProjectCompatibility', "
                "'MaterialsProject2020Compatibility', or an instance of a Compatibility class"
            )

        pbx_entries = []

        if isinstance(chemsys, str):
            chemsys = chemsys.split("-")

        # Get ion entries first, because certain ions have reference
        # solids that aren't necessarily in the chemsys (Na2SO4)
        url = "/pourbaix_diagram/reference_data/" + "-".join(chemsys)
        ion_data = self._make_request(url)
        ion_ref_comps = [Composition(d["Reference Solid"]) for d in ion_data]
        ion_ref_elts = list(
            itertools.chain.from_iterable(i.elements for i in ion_ref_comps)
        )
        ion_ref_entries = self.get_entries_in_chemsys(
            list(set([str(e) for e in ion_ref_elts] + ["O", "H"])),
            property_data=["e_above_hull"],
            compatible_only=False,
        )

        # suppress the warning about supplying the required energies; they will be calculated from the
        # entries we get from MPRester
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="You did not provide the required O2 and H2O energies.",
            )
            compat = MaterialsProjectAqueousCompatibility(
                solid_compat=self.solid_compat
            )
        # suppress the warning about missing oxidation states
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Failed to guess oxidation states.*"
            )
            ion_ref_entries = compat.process_entries(ion_ref_entries)
        ion_ref_pd = PhaseDiagram(ion_ref_entries)

        # position the ion energies relative to most stable reference state
        for n, i_d in enumerate(ion_data):
            ion = Ion.from_formula(i_d["Name"])
            refs = [
                e
                for e in ion_ref_entries
                if e.composition.reduced_formula == i_d["Reference Solid"]
            ]
            if not refs:
                raise ValueError("Reference solid not contained in entry list")
            stable_ref = sorted(refs, key=lambda x: x.data["e_above_hull"])[0]
            rf = stable_ref.composition.get_reduced_composition_and_factor()[1]

            solid_diff = (
                ion_ref_pd.get_form_energy(stable_ref)
                - i_d["Reference solid energy"] * rf
            )
            elt = i_d["Major_Elements"][0]
            correction_factor = ion.composition[elt] / stable_ref.composition[elt]
            energy = i_d["Energy"] + solid_diff * correction_factor
            ion_entry = IonEntry(ion, energy)
            pbx_entries.append(PourbaixEntry(ion_entry, f"ion-{n}"))

        # Construct the solid pourbaix entries from filtered ion_ref entries

        extra_elts = (
            set(ion_ref_elts)
            - {Element(s) for s in chemsys}
            - {Element("H"), Element("O")}
        )
        ion_ref_entries = [*ion_ref_entries, *additional_entries]
        for entry in ion_ref_entries:
            entry_elts = set(entry.composition.elements)
            # Ensure no OH chemsys or extraneous elements from ion references
            if not (
                entry_elts <= {Element("H"), Element("O")}
                or extra_elts.intersection(entry_elts)
            ):
                # Create new computed entry
                form_e = ion_ref_pd.get_form_energy(entry)
                new_entry = ComputedEntry(
                    entry.composition, form_e, entry_id=entry.entry_id
                )
                pbx_entry = PourbaixEntry(new_entry)
                pbx_entries.append(pbx_entry)

        return pbx_entries
