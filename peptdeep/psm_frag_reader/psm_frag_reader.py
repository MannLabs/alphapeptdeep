import typing
import numpy as np
import pandas as pd

from alphabase.peptide.fragment import get_charged_frag_types
from alphabase.io.psm_reader.psm_reader import PSMReaderBase

class PSMReader_w_FragBase:
    '''
    Read PSMs and fragments
    '''
    def __init__(self,
        *,
        frag_types=['b','y','b_modloss','y_modloss'], 
        max_frag_charge=2,
        **kwargs,
    ):

        self.charged_frag_types = get_charged_frag_types(
            frag_types, max_frag_charge
        )
        self._fragment_intensity_df:pd.DataFrame = pd.DataFrame(
            columns=self.charged_frag_types
        )
    
    @property
    def fragment_intensity_df(self):
        return self._fragment_intensity_df


class PSM_w_FragReaderProvider:
    def __init__(self):
        self.reader_dict = {}

    def register_reader(self, reader_name, reader_class):
        # for example, we can register the AlphaPept reader 
        self.reader_dict[reader_name.lower()] = reader_class

    def get_reader(self, reader_name, 
        frag_types=['b','y','b_modloss','y_modloss'], 
        max_frag_charge=2, **kwargs,
    )->PSMReader_w_FragBase:
        return self.reader_dict[reader_name.lower()](
            frag_types=frag_types, max_frag_charge=max_frag_charge, 
            **kwargs
        )

psm_w_frag_reader_provider = PSM_w_FragReaderProvider()
