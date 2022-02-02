from get_entries_utils import get_all_entries, write_data

if __name__ == '__main__':
    entries = get_all_entries()
    write_data(entries, 'oqmd_entries.db')