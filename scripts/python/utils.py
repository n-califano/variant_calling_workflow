import os.path

def remove_all_extensions(filename):
    while True:
        filename_new, ext = os.path.splitext(filename)
        if not ext:
            return filename_new
        filename = filename_new