import os

import openbabel


class OBConvError(Exception):
    def __init__(self, filename, msg=None):
        if msg is None:
            msg = "Error with file {}".format(filename)
        super(OBConvError, self).__init__(msg)
        self.filename = filename

class OBReadError(OBConvError):
    def __init__(self, filename):
        super(OBReadError, self).__init__(
            filename, msg="OpenBabel could not read {}. Make sure it is a"
                          " valid file".format(filename))


class OBWriteError(OBConvError):
    def __init__(self, filename):
        super(OBWriteError, self).__init__(
            filename, msg="OpenBabel can not write to {}. Make sure you have"
                          " write permissions and there is"
                          " space left".format(filename))


class OBUnrecognizedFormat(OBConvError):
    def __init__(self, extension):
        super(OBUnrecognizedFormat, self).__init__(
            None, msg =  "The extension {} is unrecognized or unsupported."
                         " Check OpenBabel website for a list of supported"
                         " formats.".format(extension))


def readFile(filename):
    """Returns an OBMol with the content of filename.

    If openbabel can not read the file, raise OBReadError. Raise
    OBUnrecognizedFormat if the file extension is unsupported."""
    input_ext = _file_extension(filename)
    in_conv = openbabel.OBConversion()
    format_status = in_conv.SetInFormat(input_ext)
    if not format_status:
        raise OBUnrecognizedFormat(input_ext)
    mol = openbabel.OBMol()
    read_status = in_conv.ReadFile(mol, filename)
    if not read_status:
        raise OBReadError(filename)
    return mol

def writeFile(filename, molecule):
    """Write molecule to file.

    Raise OBWriteError if openbabel fails to write the file; raise
    OBUnrecognizedFormat if the file extension is unsupported."""
    out_ext = _file_extension(filename)
    out_conv = openbabel.OBConversion()
    format_status = out_conv.SetOutFormat(out_ext)
    if not format_status:
        raise OBUnrecognizedFormat(out_ext)
    write_status = out_conv.WriteFile(molecule , filename)
    if not write_status:
        raise OBWriteError(filename)

def _file_extension(filename):
    """Return file extension without the dot"""
    # openbabel expects the extension without the dot, but os.path.splitext
    #  returns the extension with it
    dotext = os.path.splitext(filename)[1]
    return dotext[1:]