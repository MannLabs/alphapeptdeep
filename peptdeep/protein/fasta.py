from alphabase.protein.fasta import SpecLibFasta
from peptdeep.spec_lib.predict_lib import PredictSpecLib
from peptdeep.pretrained_models import ModelManager


class PredictSpecLibFasta(SpecLibFasta, PredictSpecLib):
    """
    Predicted spec lib from fasta files or other peptide files.
    """
    def __init__(self,
        model_manager:ModelManager = None,
        *,
        charged_frag_types:list = ['b_z1','b_z2','y_z1','y_z2'],
        protease:str = 'trypsin',
        max_missed_cleavages:int = 2,
        peptide_length_min:int = 7,
        peptide_length_max:int = 35,
        precursor_charge_min:int = 2,
        precursor_charge_max:int = 4,
        precursor_mz_min:float = 400.0, 
        precursor_mz_max:float = 1800.0,
        var_mods:list = ['Acetyl@Protein N-term','Oxidation@M'],
        min_var_mod_num:int = 0,
        max_var_mod_num:int = 2,
        fix_mods:list = ['Carbamidomethyl@C'],
        labeling_channels:dict = None,
        special_mods:list = [],
        min_special_mod_num:int = 0,
        max_special_mod_num:int = 1,
        special_mods_cannot_modify_pep_n_term:bool=False,
        special_mods_cannot_modify_pep_c_term:bool=False,
        decoy: str = None, # or pseudo_reverse or diann
        include_contaminants: bool=False,
        I_to_L=False,
        generate_precursor_isotope:bool = False,
        rt_to_irt:bool = False,
    ):
        """
        Parameters
        ----------
        model_manager : ModelManager, optional
            ModelManager of MS2/RT/CCS... models, by default None

        charged_frag_types : list, optional
            Fragment types with charge, 
            by default [ 'b_z1','b_z2','y_z1', 'y_z2' ]

        protease : str, optional
            Could be pre-defined protease name defined in :data:`protease_dict`,
            or a regular expression. 
            By default 'trypsin'

        max_missed_cleavages : int, optional
            Maximal missed cleavages, by default 2
            
        peptide_length_min : int, optional
            Minimal cleaved peptide length, by default 7

        peptide_length_max : int, optional
            Maximal cleaved peptide length, by default 35

        precursor_charge_min : int, optional
            Minimal precursor charge, by default 2

        precursor_charge_max : int, optional
            Maximal precursor charge, by default 4

        precursor_mz_min : float, optional
            Minimal precursor mz, by default 200.0

        precursor_mz_max : float, optional
            Maximal precursor mz, by default 2000.0

        var_mods : list, optional
            list of variable modifications, 
            by default ['Acetyl@Protein N-term','Oxidation@M']

        max_var_mod_num : int, optional
            Minimal number of variable modifications on a peptide sequence, 
            by default 0

        max_var_mod_num : int, optional
            Maximal number of variable modifications on a peptide sequence, 
            by default 2

        fix_mods : list, optional
            list of fixed modifications, by default ['Carbamidomethyl@C']

        labeling_channels : dict, optional
            Add isotope labeling with different channels, 
            see :meth:`add_peptide_labeling()`. 
            Defaults to None

        special_mods : list, optional
            Special modifications.
            It is useful for modificaitons like Phospho which may largely 
            explode the number of candidate modified peptides.
            The number of special_mods per peptide 
            is controlled by `max_append_mod_num`.
            Defaults to [].

        min_special_mod_num : int, optional
            Control the min number of special_mods per peptide, by default 0.

        max_special_mod_num : int, optional
            Control the max number of special_mods per peptide, by default 1.

        special_mods_cannot_modify_pep_c_term : bool, optional
            Some modifications cannot modify the peptide C-term, 
            this will be useful for GlyGly@K as if C-term is di-Glyed, 
            it cannot be cleaved/digested. 
            Defaults to False.

        special_mods_cannot_modify_pep_n_term : bool, optional
            Similar to `special_mods_cannot_modify_pep_c_term`, but at N-term.
            Defaults to False.

        decoy : str, optional
            Decoy type, see `alphabase.spectral_library.decoy_library`,
            by default None

        include_contaminants : bool, optional
            If include contaminants.fasta, by default False

        generate_precursor_isotope : bool, optional
            If :meth:`peptdeep.spec_lib.predict_lib.PredictSpecLib.predict_all()` 
            includes :meth:`peptdeep.spec_lib.predict_lib.PredictSpecLib.calc_precursor_isotope()`.
            Defaults to False
        
        rt_to_irt : bool, optional
            If convert predicted RT to iRT values
        """
        SpecLibFasta.__init__(self,
            charged_frag_types=charged_frag_types,
            protease=protease,
            max_missed_cleavages=max_missed_cleavages,
            peptide_length_min=peptide_length_min,
            peptide_length_max=peptide_length_max,
            precursor_charge_min=precursor_charge_min,
            precursor_charge_max=precursor_charge_max,
            precursor_mz_min=precursor_mz_min, 
            precursor_mz_max=precursor_mz_max,
            var_mods=var_mods,
            min_var_mod_num=min_var_mod_num,
            max_var_mod_num=max_var_mod_num,
            fix_mods=fix_mods,
            labeling_channels=labeling_channels,
            special_mods=special_mods,
            min_special_mod_num=min_special_mod_num,
            max_special_mod_num=max_special_mod_num,
            special_mods_cannot_modify_pep_n_term=special_mods_cannot_modify_pep_n_term,
            special_mods_cannot_modify_pep_c_term=special_mods_cannot_modify_pep_c_term,
            decoy=decoy,
            include_contaminants=include_contaminants,
            I_to_L=I_to_L,
        )

        PredictSpecLib.__init__(self,
            model_manager=model_manager,
            charged_frag_types=self.charged_frag_types,
            precursor_mz_min=self.min_precursor_mz,
            precursor_mz_max=self.max_precursor_mz,
            decoy=self.decoy,
            generate_precursor_isotope=generate_precursor_isotope,
            rt_to_irt=rt_to_irt,
        )
