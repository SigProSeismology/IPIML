import os
import glob
from obspy import read
from obspy import UTCDateTime



import os

def find_files_and_merge(directory_path, start_time_str, end_time_str, output_directory_path=None):
    start_time = UTCDateTime(start_time_str)
    end_time = UTCDateTime(end_time_str)

    current_day = UTCDateTime(start_time.year, start_time.month, start_time.day)
    end_day = UTCDateTime(end_time.year, end_time.month, end_time.day)

    if output_directory_path is not None:
        if not os.path.exists(output_directory_path):
            os.makedirs(output_directory_path)

    num_files_merged = 0
    while current_day <= end_day:
        current_day_str = current_day.strftime("%Y%m%d")
        file_pattern = f"*__{current_day_str}T*.mseed"
        file_path = os.path.join(directory_path, file_pattern)
        files = glob.glob(file_path)

        components = set([os.path.basename(file).split("..")[1][-45:-42] for file in files])

        for component in components:
            component_files = [file for file in files if os.path.basename(file).split("..")[1][-45:-42] == component]
            st = read(component_files[0])

            for file in component_files[1:]:
                st += read(file)

            st.merge(method=1, fill_value='latest')
            st.trim(start_time, end_time)

            num_files_merged += 1
            if output_directory_path is not None:
                output_file_name = f"{st[0].stats.network}.{st[0].stats.station}..{component}__{start_time_str}_{end_time_str}.mseed"
                output_file_path = os.path.join(output_directory_path, output_file_name)
                st.write(output_file_path, format="MSEED")

        current_day += 86400

    return num_files_merged

directory_path = "./AAA/ori_data/"
output_directory_path = "./AAA/primary_events/"
start_time_str = "20121020T235800Z"
end_time_str = "20121021T025900Z"
# output_directory_path = "./AAA/primary_events/AKOS"
# find_files_and_merge(directory_path, start_time_str, end_time_str, output_directory_path)

# Get a list of all subdirectories in the directory_path
subdirectories = [os.path.join(directory_path, d) for d in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, d))]

for subdirectory in subdirectories:
    subdirectory_name = os.path.basename(subdirectory)
    subdirectory_output_path = os.path.join(output_directory_path, subdirectory_name)

    if not os.path.exists(subdirectory_output_path):
        os.makedirs(subdirectory_output_path)

    find_files_and_merge(subdirectory, start_time_str, end_time_str, subdirectory_output_path)



# Example usage: Merge all files in directory "/path/to/files" that contain data within the time window 20121020T103000Z - 20121020T113200Z
# directory_path = "./AAA/ori_data/AKOS"

