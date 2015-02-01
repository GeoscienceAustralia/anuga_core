"""Function to return the length of a file.

Intended to be used for simplfying the creation of boundary_tags in
run_model.py

Input: filename
Returns: number of lines in file
"""
def file_length(in_file):
    fid = open(in_file)
    data = fid.readlines()
    fid.close()
    return len(data)

