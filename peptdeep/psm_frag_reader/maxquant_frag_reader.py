import pandas as pd
import numpy as np
import numba

from alphabase.peptide.fragment import (
    init_fragment_by_precursor_dataframe
)
from alphabase.io.psm_reader.maxquant_reader import (
    MaxQuantReader
)

from peptdeep.psm_frag_reader.psm_frag_reader import (
    PSMReader_w_FragBase,
    psm_w_frag_reader_provider
)

def parse_phos_probs(mods, prob_seq, prob):
    if mods == 'Unmodified': return '', ''
    idx = mods.find('Phospho')
    if idx == -1: return '', ''
    elif idx == 0 or mods[idx-1]==',':
        num_phos = 1
    elif mods[idx-2].isdigit():
        num_phos = int(mods[idx-2])
    else:
        return 'x', 'x'
    
    idx = prob_seq.find('(')
    keep_probs = []
    keep_sites = []
    while idx != -1:
        end = prob_seq.find(')',idx+2)
        if prob_seq[idx-1] in 'STY':
            if float(prob_seq[idx+1:end])>=prob:
                keep_probs.append(prob_seq[idx+1:end])
                keep_sites.append(str(idx))
        prob_seq = prob_seq[:idx]+prob_seq[end+1:]
        idx = prob_seq.find('(',idx)
    if len(keep_probs) >= num_phos: 
        return ';'.join(keep_probs),';'.join(keep_sites)
    else: return 'x','x'

def filter_phos(mq_df, prob):
    if 'Phospho (STY) Probabilities' not in mq_df.columns:
        mq_df['PhosProbs'] = ''
        mq_df['PhosSites'] = ''
        return mq_df

    (
        mq_df['PhosProbs'], mq_df['PhosSites']
    ) = zip(*mq_df[['Modifications','Phospho (STY) Probabilities']].apply(
        lambda x: parse_phos_probs(x[0],x[1],prob), axis=1
    ))
    return mq_df[mq_df['PhosProbs']!='x']
    

class MaxQuantMSMSReader(MaxQuantReader, PSMReader_w_FragBase):
    def __init__(self,
        frag_types=['b','y','b_modloss','y_modloss'], 
        max_frag_charge=2,
        score_threshold=100,
        rt_unit='minute',
        **kwargs
    ):
        PSMReader_w_FragBase.__init__(self,
            frag_types = frag_types,
            max_frag_charge = max_frag_charge,
            **kwargs
        )

        MaxQuantReader.__init__(self, rt_unit=rt_unit)

        self.column_mapping['phos_probs'] = 'PhosProbs'
        self.column_mapping['phos_sites'] = 'PhosSites'
        self._score_thres = score_threshold
        self._phos_prob = 0.75
    
    @property
    def fragment_intensity_df(self):
        return self._fragment_intensity_df

    def _load_file(self, filename):
        df = MaxQuantReader._load_file(self, filename)
        df = filter_phos(df, self._phos_prob)
        df = df[
            (df.Score >= self._score_thres)|
            ((df.Score>=60)&(df.PhosProbs!=''))]
        df.reset_index(drop=True, inplace=True)
        return df

    def _post_process(self, 
        mq_df
    ):  
        self._psm_df['nAA'] = self._psm_df.sequence.str.len()
        mq_df['nAA'] = self._psm_df.nAA
        self.normalize_rt_by_raw_name()

        self._fragment_intensity_df = init_fragment_by_precursor_dataframe(
            mq_df, self.charged_frag_types
        )

        frag_col_dict = dict(zip(
            self.charged_frag_types, 
            range(len(self.charged_frag_types))
        ))

        for ith_psm, (nAA, start,end) in enumerate(
            mq_df[['nAA','frag_start_idx','frag_stop_idx']].values
        ):
            intens = np.zeros((nAA-1, len(self.charged_frag_types)))

            frag_types = mq_df['Matches'].values[ith_psm]
            frag_intens = mq_df['Intensities'].values[ith_psm]
            for frag_type, frag_inten in zip(
                frag_types.split(';'), frag_intens.split(';')
            ):
                if '-' in frag_type: continue
                if any(_.isupper() for _ in frag_type): continue
                idx = frag_type.find('(')
                charge = '1'
                if idx > 0:
                    frag_type, charge = frag_type[:idx], frag_type[idx+1:-2]
                if not frag_type[1].isdigit(): continue # no H2O or NH3 loss
                frag_type, frag_pos = frag_type[0], frag_type[1:]
                if frag_pos.endswith('*'):
                    frag_pos = int(frag_pos[:-1])
                    modloss=True
                else: 
                    frag_pos = int(frag_pos)
                    modloss=False
                if frag_type in 'xyz':
                    frag_pos = nAA - frag_pos -1
                else:
                    frag_pos -= 1 
                frag_type += ('_modloss_z' if modloss else '_z') +charge
                if frag_type not in frag_col_dict: continue
                frag_col = frag_col_dict[frag_type]
                intens[frag_pos,frag_col] = float(frag_inten)

            if np.any(intens>0):
                intens /= np.max(intens)
            self._fragment_intensity_df.iloc[
                start:end,:
            ] = intens
        
        self._psm_df[
            ['frag_start_idx','frag_stop_idx']
        ] = mq_df[['frag_start_idx','frag_stop_idx']]

        self._psm_df = self._psm_df[~self._psm_df.mods.isna()]


psm_w_frag_reader_provider.register_reader('maxquant', MaxQuantMSMSReader)
