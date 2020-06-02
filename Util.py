import os

def FileValidator(file):
    if os.path.exists(file) == False:
        raise Exception("\nERROR:File: " + file + " Doesn't exists")
    if os.path.isfile(file) == False:
        raise Exception('\nERROR:File expected:' + file + ':is not a file')
    if os.stat(file)[6] == 0:
        raise Exception('\nERROR:File:' + file + ':is empty')
    else:
        lines = 0
        for i in open(file):
            if i.startswith("#"):
                continue
            lines = 1
            break
        if lines == 0:
            raise Exception('\nERROR:File:' + file + ':is empty')

    if os.access(file, os.R_OK) == False:
        raise Exception('\nERROR:File:\n' + file + ':no read access ')
    return file


def DirectoryValidator(file, read_write="NO"):
    if os.path.exists(os.path.dirname(file)) == False:
        raise Exception('\nERROR:Path:' + file + '\n:Does not exist')
    if read_write == "YES":
        if os.access(os.path.dirname(file), os.W_OK) == False:
            raise Exception('\nERROR:File:\n' + os.path.dirname(file) + ':no write access ')
    if os.path.isfile(file) == True:
        raise Exception('\nERROR:Directory expected:\n' + file + ':is a file')
    return file



def ensure_dir(f):
    """ Checks if a directory is present or not and creates one if not present """
    if not os.path.exists(f):
        os.makedirs(f)
