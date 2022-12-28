import pandas as pd
import numpy as np
import typing
from tqdm import tqdm

import alphabase.constants.modification as ap_mod
from alphabase.peptide.fragment import (
    concat_precursor_fragment_dataframes,
    init_fragment_by_precursor_dataframe
)

from alphabase.io.psm_reader.pfind_reader import (
    pFindReader, get_pFind_mods, translate_pFind_mod
)

from peptdeep.psm_frag_reader.psm_frag_reader import (
    PSMReader_w_FragBase, psm_w_frag_reader_provider
)


class PSMLabelReader(pFindReader,PSMReader_w_FragBase):
    def __init__(self, 
        frag_types=['b','y','b_modloss','y_modloss'], 
        max_frag_charge=2,
        **kwargs
    ):
        PSMReader_w_FragBase.__init__(self,
            frag_types=frag_types,
            max_frag_charge=max_frag_charge,
            **kwargs
        )
        pFindReader.__init__(self)

        psmlabel_columns = 'b,b-NH3,b-H20,b-ModLoss,y,y-HN3,y-H20,y-ModLoss'.split(',')
        self.psmlabel_frag_columns = []
        self.frag_df_columns = {}
        for _type in psmlabel_columns:
            frag_idxes = [
                i for i,_t in enumerate(
                    self.charged_frag_types
                ) if _t.startswith(_type.replace('-','_').lower()+'_')
            ]
            if frag_idxes:
                self.psmlabel_frag_columns.append(_type)
                self.frag_df_columns[_type] = np.array(
                    frag_idxes, dtype=int
                )

    def _init_column_mapping(self):
        self.column_mapping = {
            'sequence': 'peptide',
            'charge': 'charge',
            'rt': 'RT',
            'raw_name': 'raw_name',
            'query_id': 'spec',
            'scan_num': 'scan_num',
        }

    def _load_file(self, filename):
        psmlabel_df = pd.read_csv(filename, sep="\t")
        psmlabel_df.fillna('', inplace=True)

        if psmlabel_df['spec'].values[0].count('.')>=4: #pfind
            psmlabel_df['raw_name'], psmlabel_df['scan_num'], psmlabel_df['charge'] = zip(
                *psmlabel_df['spec'].str.split('.').apply(lambda x: (x[0], x[-4], x[-3]))
            )
        else:
            psmlabel_df['raw_name'], psmlabel_df['scan_num'] = zip(
                *psmlabel_df['spec'].str.split('.').apply(lambda x: (x[0], x[1]))
            )
        psmlabel_df['scan_num'] = psmlabel_df['scan_num'].astype(int)
        psmlabel_df['charge'] = psmlabel_df['charge'].astype(int)
        return psmlabel_df

    def _load_modifications(self, psmlabel_df: pd.DataFrame):
        (
            self._psm_df['mods'], self._psm_df['mod_sites'] 
        ) = zip(*psmlabel_df['modinfo'].apply(get_pFind_mods))

    def _translate_decoy(self, df):
        pass
    def _translate_score(self, df):
        pass

    def _translate_modifications(self):
        self._psm_df['mods'] = self._psm_df['mods'].apply(translate_pFind_mod)

    def _post_process(self, psmlabel_df: pd.DataFrame):
        psmlabel_df['nAA'] = psmlabel_df.peptide.str.len()
        self._psm_df['nAA'] = psmlabel_df.nAA
        psmlabel_df = psmlabel_df[
            ~self._psm_df['mods'].isna()
        ].reset_index(drop=True)
        
        self._psm_df = self._psm_df[
            ~self._psm_df['mods'].isna()
        ].reset_index(drop=True)

        self._fragment_intensity_df = init_fragment_by_precursor_dataframe(
            psmlabel_df, self.charged_frag_types
        )

        for ith_psm, (nAA, start,end) in enumerate(
            psmlabel_df[['nAA','frag_start_idx','frag_stop_idx']].values
        ):
            intens = np.zeros((nAA-1, len(self.charged_frag_types)))
            for ion_type in self.psmlabel_frag_columns:
                if ion_type not in psmlabel_df.columns: continue

                pos_end = ion_type.find('-')-len(ion_type)-2 if '-' in ion_type else -2
                typed_frags = psmlabel_df.loc[ith_psm,ion_type]
                if not typed_frags: continue
                typed_frags = typed_frags.strip(';').split(';')
                frag_pos = []
                frag_charge = []
                frag_inten = []

                for frag in typed_frags:
                    frag, inten = frag.split(',')
                    frag_pos.append(int(frag[1:pos_end]))
                    frag_charge.append(int(frag[-1]))
                    frag_inten.append(float(inten))
                if not frag_inten: continue
                
                frag_pos = np.array(frag_pos, dtype=int)
                frag_col = np.array(frag_charge, dtype=int)-1
                
                if ion_type[0] in 'xyz':
                    frag_pos = nAA - frag_pos -1
                else:
                    frag_pos -= 1
                intens[frag_pos,self.frag_df_columns[ion_type][frag_col]] = frag_inten
            if np.any(intens>0):
                intens /= np.max(intens)
            self._fragment_intensity_df.iloc[
                start:end,:
            ] = intens
        
        self._psm_df[
            ['frag_start_idx','frag_stop_idx']
        ] = psmlabel_df[['frag_start_idx','frag_stop_idx']]

psm_w_frag_reader_provider.register_reader(
    'psmlabel', PSMLabelReader
)

def load_psmlabel_list(
    psmlabel_list,
    nce_list,
    instrument_list,
    frag_types=['b','y','b_modloss','y_modloss'], 
    frag_charge=2,
    include_mod_list=[
        'Oxidation@M','Phospho@S','Phospho@T','Phospho@Y','Acetyl@Protein N-term'
    ]
):
    psm_df_list = []
    fragment_inten_df_list = []
    for i,psmlabel in tqdm(enumerate(psmlabel_list)):
        psm_reader = PSMLabelReader(
            frag_types=frag_types,
            max_frag_charge = frag_charge
        )
        psm_reader.import_file(psmlabel)
        psm_reader.filter_psm_by_modifications(include_mod_list)
        psm_reader.psm_df['nce'] = nce_list[i]
        psm_reader.psm_df['instrument'] = instrument_list[i]
        psm_df_list.append(psm_reader.psm_df)
        fragment_inten_df_list.append(psm_reader.fragment_intensity_df)
    return concat_precursor_fragment_dataframes(psm_df_list, fragment_inten_df_list)

