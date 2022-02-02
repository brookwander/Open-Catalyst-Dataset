from get_entries_utils import get_icsd_entries, write_data

if __name__ == "__main__":
    entries = get_icsd_entries()
    write_data(entries, "oqmd_icsd_entries.db")
