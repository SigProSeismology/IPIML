import glob
from obspy import read
from obspy import UTCDateTime

def find_files_and_trim(directory_path, start_time_str, end_time_str):
    # Parse the start and end times as UTCDateTime objects
    start_time = UTCDateTime(start_time_str)
    end_time = UTCDateTime(end_time_str)

    # Find files that contain the date of interest
    pattern = os.path.join(directory_path, f"*__{start_time.date.strftime('%Y%m%d')}*.mseed")
    files = glob.glob(pattern)

    # Loop over the files and trim the data within the desired time window
    for file_path in files:
        # Read the file using ObsPy
        st = read(file_path)

        # Check if the file contains data within the desired time window
        if start_time.date == end_time.date:
            st_trimmed = st.trim(starttime=start_time, endtime=end_time)
        else:
            st_trimmed = st.slice(starttime=start_time, endtime=end_time)

        # Write the trimmed data to a new file
        file_name, file_extension = os.path.splitext(file_path)
        output_file_path = f"{file_name}_{start_time_str}_{end_time_str}{file_extension}"
        st_trimmed.write(output_file_path, format="MSEED")

    return len(files)

# Example usage: Trim all files in directory "/path/to/files" that contain data within the time window 20121020T103000Z - 20121020T113200Z
directory_path = "/path/to/files"
start_time_str = "20121020T103000Z"
end_time_str = "20121020T113200Z"
num_files_trimmed = find_files_and_trim(directory_path, start_time_str, end_time_str)
print(f"Trimmed {num_files_trimmed} files.")