import logging
import numpy as np
import mrcfile
import struct

from pyworkflow.utils import cyanStr, redStr

logger = logging.getLogger(__name__)


class EMFileHandler:
    """
    Class to handle basic reading and writing of .em format files.
    Based on the standard TOM toolbox header structure.
    """
    HEADER_SIZE = 512

    # Data type codes in .em format
    EM_BYTE = 1
    EM_SHORT = 2
    EM_LONG = 4
    EM_FLOAT = 5
    EM_DOUBLE = 9

    @staticmethod
    def read(filepath: str):
        """Reads an .em file and returns a numpy array."""
        with open(filepath, 'rb') as f:
            header = f.read(EMFileHandler.HEADER_SIZE)
            if len(header) < EMFileHandler.HEADER_SIZE:
                raise ValueError("File is too small to be a valid .em file")

            # Dimensions (x, y, z) are stored at bytes 4, 8, 12
            # These are 32-bit integers
            dims = struct.unpack('<3i', header[4:16])
            width, height, depth = dims

            # Data type code is at byte 3
            data_type_code = header[3]
            if data_type_code == EMFileHandler.EM_BYTE:
                dtype = np.int8
            elif data_type_code == EMFileHandler.EM_SHORT:
                dtype = np.int16
            elif data_type_code == EMFileHandler.EM_LONG:
                dtype = np.int32
            elif data_type_code == EMFileHandler.EM_FLOAT:
                dtype = np.float32
            elif data_type_code == EMFileHandler.EM_DOUBLE:
                dtype = np.float64
            else:
                # Fallback: defaults to float32 if type is unknown
                logger.info(cyanStr(f"Warning: Unknown data type code {data_type_code}. Assuming float32."))
                dtype = np.float32

            # Read remaining binary data
            data = np.fromfile(f, dtype=dtype)

            # Verify size consistency
            expected_len = width * height * depth
            if data.size != expected_len:
                raise ValueError(f"Incorrect data size. Expected {expected_len}, got {data.size}")

            # Reshape to correct volume dimensions.
            # .em uses Fortran order (X, Y, Z), numpy uses C order (Z, Y, X).
            # We read in Fortran order and then transpose to (Z, Y, X).
            volume = data.reshape((width, height, depth), order='F')

            # Transpose to match Python/MRC convention (Z, Y, X)
            volume = volume.transpose((2, 1, 0))

            return volume

    @staticmethod
    def write(filepath, data):
        """Writes a numpy array to .em format."""
        # Ensure array is (Z, Y, X) -> Transpose back to (X, Y, Z) for .em
        data_to_write = data.transpose((2, 1, 0))

        # Get dimensions
        width, height, depth = data_to_write.shape  # X, Y, Z

        header = bytearray(EMFileHandler.HEADER_SIZE)

        # Byte 0: Machine coding (6 = PC)
        header[0] = 6
        # Byte 3: Data type
        if data.dtype == np.float32:
            header[3] = EMFileHandler.EM_FLOAT
        elif data.dtype == np.int16:
            header[3] = EMFileHandler.EM_SHORT
        elif data.dtype == np.int8:
            header[3] = EMFileHandler.EM_BYTE
        else:
            # Convert to float32 by default for other types
            data_to_write = data_to_write.astype(np.float32)
            header[3] = EMFileHandler.EM_FLOAT

        # Pack dimensions (x, y, z) into bytes 4, 8, 12
        struct.pack_into('<3i', header, 4, width, height, depth)

        with open(filepath, 'wb') as f:
            f.write(header)
            # Write data in Fortran order (column-major)
            data_to_write.tofile(f, sep="", format="")


# --- Conversion Functions ---

def em_to_mrc(em_path, mrc_path, pixel_size=None):
    """Converts .em file to .mrc"""
    try:
        volume_data = EMFileHandler.read(em_path)
        with mrcfile.new(mrc_path, overwrite=True) as mrc:
            mrc.set_data(volume_data)
            if pixel_size:
                mrc.voxel_size = pixel_size

    except Exception as e:
        logger.error(redStr(f"Error in em_to_mrc: {e}"))


def mrc_to_em(mrc_path, em_path):
    """Converts .mrc file to .em"""
    try:
        # Use permissive=True for MRC files that might have slightly corrupt headers
        with mrcfile.open(mrc_path, permissive=True) as mrc:
            volume_data = mrc.data

        EMFileHandler.write(em_path, volume_data)

    except Exception as e:
        logger.error(redStr(f"Error in mrc_to_em: {e}"))
