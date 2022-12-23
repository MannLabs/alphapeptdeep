import os
import numpy as np
import pandas as pd
from alphabase.io.hdf import HDF_File

try:
    # should be replaced by AlphaRaw in the near future
    from peptdeep.legacy.thermo_raw.pyrawfilereader import RawFileReader
except (ImportError,AttributeError):
    RawFileReader = None

class MSReaderBase:
    def __init__(self):
        self.spectrum_df:pd.DataFrame = pd.DataFrame()
        self.peak_df:pd.DataFrame = pd.DataFrame()
        # self.mzs: np.ndarray = np.array([])
        # self.intensities: np.ndarray = np.array([])

    def load(self, file_path):
        raise NotImplementedError('load()')

    def build_spectrum_df(self,
        scan_list:list,
        scan_indices:np.ndarray,
        rt_list:list,
        mobility_list:list = None,
    ):
        """Build spectrum_df by the given information

        Parameters
        ----------
        scan_list : list
            scan number list

        scan_indices : np.array
            starts and end positions of ms2
            peaks for each scan

        rt_list : list
            retention time (minutes) for each scan 

        mobility_list : list, optional
            mobility for each scan. Defaults to None.
            
        """
        def set_col(col, indexes, values, dtype, na_value):
            self.spectrum_df.loc[indexes, col] = values
            self.spectrum_df[col].fillna(na_value, inplace=True)
            self.spectrum_df[col] = self.spectrum_df[col].astype(dtype)

        scan_list = np.array(scan_list, dtype=np.int32)
        if scan_list.min() > 0:
            # thermo scan >= 1
            scan_list -= 1
        idx_len = np.max(scan_list)+1
        self.spectrum_df = pd.DataFrame(index=np.arange(idx_len, dtype=np.int64))
        self.spectrum_df['spec_idx'] = self.spectrum_df.index.values
        set_col('peak_start_idx', scan_list, scan_indices[:-1], np.int64, -1)
        set_col('peak_end_idx', scan_list, scan_indices[1:], np.int64, -1)
        set_col('rt', scan_list, rt_list, np.float64, np.nan)
        if mobility_list is not None:
            set_col('mobility', scan_list, mobility_list, np.float64, np.nan)

    def get_peaks(self, spec_idx:int):
        """Get peak (mz and intensity) values by `spec_idx`

        Parameters
        ----------
        spec_idx : int
            indicator for a spectrum, 
            could be scan_num-1 for thermo data.

        Returns
        -------
        np.array
            mz values for the given spec_idx

        np.array
            intensity values for the given spec_idx

        """
        if spec_idx not in self.spectrum_df.index:
            return None, None
        start_idx, end_idx = self.spectrum_df.loc[
            spec_idx, ['peak_start_idx','peak_end_idx']
        ].values.astype(np.int64)
        return (
            self.peak_df.mz.values[start_idx:end_idx],
            self.peak_df.intensity.values[start_idx:end_idx]
        )

    def get_peaks_by_scan_num(self, scan_num:int):
        """Get peak (mz and intensity) values by `spec_idx`

        Parameters
        ----------
        scan_num : int
            scan_num of thermodata

        Returns
        -------
        np.array
            mz values for the given spec_idx (scan)

        np.array
            intensity values for the given spec_idx
            
        """
        return self.get_peaks(scan_num-1)

class AlphaPept_HDF_MS1_Reader(MSReaderBase):
    """MS1 from AlphaPept HDF"""
    def load(self, file_path):
        hdf = HDF_File(file_path)
        self.peak_df['mz'] = hdf.Raw.MS1_scans.mass_list_ms1.values
        self.peak_df['intensity'] = hdf.Raw.MS1_scans.int_list_ms1.values
        self.build_spectrum_df(
            scan_list=hdf.Raw.MS1_scans.scan_list_ms1.values, 
            scan_indices=hdf.Raw.MS1_scans.indices_ms1.values,
            rt_list=hdf.Raw.MS1_scans.rt_list_ms1.values,
            mobility_list=hdf.Raw.MS1_scans.mobility.values 
            if hasattr(hdf.Raw.MS1_scans, 'mobility') else None,
        )

class AlphaPept_HDF_MS2_Reader(MSReaderBase):
    """MS2 from AlphaPept HDF"""
    def load(self, file_path):
        hdf = HDF_File(file_path)
        self.peak_df['mz'] = hdf.Raw.MS2_scans.mass_list_ms2.values
        self.peak_df['intensity'] = hdf.Raw.MS2_scans.int_list_ms2.values
        if hasattr(hdf.Raw.MS2_scans, 'mobility2'):
            scan_list = np.arange(len(hdf.Raw.MS2_scans.rt_list_ms2))
        else:
            scan_list = hdf.Raw.MS2_scans.scan_list_ms2.values
        self.build_spectrum_df(
            scan_list=scan_list, 
            scan_indices=hdf.Raw.MS2_scans.indices_ms2.values,
            rt_list=hdf.Raw.MS2_scans.rt_list_ms2.values,
            mobility_list=hdf.Raw.MS2_scans.mobility2.values 
            if hasattr(hdf.Raw.MS2_scans, 'mobility2') else None,
        )

def read_until(file, until):
    lines = []
    while True:
        line = file.readline().strip()
        if line == "": break
        elif line.startswith(until):
            break
        else:
            lines.append(line)
    return lines

def find_line(lines, start):
    for line in lines:
        if line.startswith(start):
            return line
    return None

def parse_pfind_scan_from_TITLE(pfind_title):
    return int(pfind_title.split('.')[-4])

def is_pfind_mgf(mgf):
    return mgf.upper().endswith('_HCDFT.MGF')

def index_ragged_list(ragged_list: list)  -> np.ndarray:
    """Create lookup indices for a list of arrays for concatenation.

    Parameters
    ----------
    value : list
        Input list of arrays.

    Returns
    -------
    indices
        A numpy array with indices.
        
    """
    indices = np.zeros(len(ragged_list) + 1, np.int64)
    indices[1:] = [len(i) for i in ragged_list]
    indices = np.cumsum(indices)

    return indices

class MGFReader(MSReaderBase):
    """MGF Reader (MS2)"""
    def load(self, mgf):
        if isinstance(mgf, str):
            f = open(mgf)
        else:
            f = mgf
        scanset = set()
        masses_list = []
        intens_list = []
        scan_list = []
        rt_list = []
        while True:
            line = f.readline()
            if not line: break
            if line.startswith('BEGIN IONS'):
                lines = read_until(f, 'END IONS')
                masses = []
                intens = []
                scan = None
                RT = 0
                for line in lines:
                    if line[0].isdigit():
                        mass,inten = [float(i) for i in line.strip().split()]
                        masses.append(mass)
                        intens.append(inten)
                    elif line.startswith('SCAN='):
                        scan = int(line.split('=')[1])
                    elif line.startswith('RTINSECOND'):
                        RT = float(line.split('=')[1])/60
                if not scan:
                    title = find_line(lines, 'TITLE=')
                    scan = parse_pfind_scan_from_TITLE(title)
                if scan in scanset: continue
                scanset.add(scan)
                scan_list.append(scan)
                rt_list.append(RT)
                masses_list.append(np.array(masses))
                intens_list.append(np.array(intens))
        if isinstance(mgf, str): 
            f.close()
        self.build_spectrum_df(
            scan_list, 
            index_ragged_list(masses_list), 
            rt_list
        )
        self.peak_df['mz'] = np.concatenate(masses_list)
        self.peak_df['intensity'] = np.concatenate(intens_list)

class MSReaderProvider:
    """Factory class to register and get MS Readers"""
    def __init__(self):
        self.reader_dict = {}
    def register_reader(self, ms2_type, reader_class):
        self.reader_dict[ms2_type.lower()] = reader_class

    def get_reader(self, file_type)->MSReaderBase:
        if file_type not in self.reader_dict: return None
        else: return self.reader_dict[file_type.lower()]()

ms2_reader_provider = MSReaderProvider()
ms2_reader_provider.register_reader('mgf', MGFReader)
ms2_reader_provider.register_reader('alphapept', AlphaPept_HDF_MS2_Reader)
ms2_reader_provider.register_reader('alphapept_hdf', AlphaPept_HDF_MS2_Reader)

ms1_reader_provider = MSReaderProvider()
ms1_reader_provider.register_reader('alphapept', AlphaPept_HDF_MS1_Reader)
ms1_reader_provider.register_reader('alphapept_hdf', AlphaPept_HDF_MS1_Reader)

if RawFileReader is None:
    class ThermoRawMS1Reader:
        def __init__(self):
            raise NotImplementedError("RawFileReader is not available")
    
    class ThermoRawMS2Reader:
        def __init__(self):
            raise NotImplementedError("RawFileReader is not available")
else:
    class ThermoRawMS1Reader(MSReaderBase):
        """Thermo Raw MS1 Reader"""
        def __init__(self):
            super().__init__()
            self.profile_mode = False

        def load(self, raw_path):
            rawfile = RawFileReader(raw_path)
 
            spec_indices = np.array(
                range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1)
            )
            scan_list = []
            rt_list = []
            masses_list = []
            intens_list = []
            for i in spec_indices:
                try:
                    ms_order = rawfile.GetMSOrderForScanNum(i)

                    if ms_order == 1:
                        if self.profile_mode:
                            masses, intens = rawfile.GetProfileMassListFromScanNum(i)
                        else:
                            masses, intens = rawfile.GetCentroidMassListFromScanNum(i)
                        scan_list.append(i)
                        rt_list.append(rawfile.RTInSecondsFromScanNum(i))
                        masses_list.append(masses)
                        intens_list.append(intens)

                except KeyboardInterrupt as e:
                    raise e
                except SystemExit as e:
                    raise e
                except Exception as e:
                    print(f"Bad scan={i} in raw file '{raw_path}'")
            
            self.build_spectrum_df(
                scan_list,
                index_ragged_list(masses_list),
                rt_list,
            )
            self.peak_df['mz'] = np.concatenate(masses_list)
            self.peak_df['intensity'] = np.concatenate(intens_list)
            rawfile.Close()

    class ThermoRawMS2Reader(MSReaderBase):
        """Thermo RAW MS2 Reader"""
        def __init__(self):
            super().__init__()
            self.profile_mode = False

        def load(self, raw_path):
            rawfile = RawFileReader(raw_path)
 
            spec_indices = np.array(
                range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1)
            )
            scan_list = []
            rt_list = []
            masses_list = []
            intens_list = []
            for i in spec_indices:
                try:
                    ms_order = rawfile.GetMSOrderForScanNum(i)

                    if ms_order == 2:
                        if self.profile_mode:
                            masses, intens = rawfile.GetProfileMassListFromScanNum(i)
                        else:
                            masses, intens = rawfile.GetCentroidMassListFromScanNum(i)
                        scan_list.append(i)
                        rt_list.append(rawfile.RTFromScanNum(i))
                        masses_list.append(masses)
                        intens_list.append(intens)

                except KeyboardInterrupt as e:
                    raise e
                except SystemExit as e:
                    raise e
                # except Exception as e:
                #     print(f"Bad scan={i} in raw file '{raw_path}'")
            
            self.build_spectrum_df(
                scan_list,
                index_ragged_list(masses_list),
                rt_list,
            )
            self.peak_df['mz'] = np.concatenate(masses_list)
            self.peak_df['intensity'] = np.concatenate(intens_list)
            rawfile.Close()
    
    ms2_reader_provider.register_reader('thermo', ThermoRawMS2Reader)
    ms2_reader_provider.register_reader('thermo_raw', ThermoRawMS2Reader)
    ms1_reader_provider.register_reader('thermo', ThermoRawMS1Reader)
    ms1_reader_provider.register_reader('thermo_raw', ThermoRawMS1Reader)
