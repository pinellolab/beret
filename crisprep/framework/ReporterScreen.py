from typing import List, Union, Collection
from functools import reduce
from copy import deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
import anndata as ad
from perturb_tools import Screen
from .Edit import Edit, Allele
from ..framework._supporting_fn import get_aa_alleles, filter_allele_by_pos


def _get_counts(filename_pattern, guides, reps, conditions):
    for rep in reps:
        for cond in conditions:
            res = pd.read_csv(filename_pattern.format(r=rep, c=cond), delimiter="\t")
            res = res.set_index("guide_name")[["read_counts"]]
            res.columns = res.columns.map(lambda x: "{r}_{c}".format(r=rep, c=cond))
            guides = guides.join(res, how="left")
        else:
            guides = guides.fillna(0)
            return guides


def _get_edits(filename_pattern, guide_info, reps, conditions, count_exact=True):
    edits = pd.DataFrame(index=(guide_info.index))
    for rep in reps:
        for cond in conditions:
            res = pd.read_csv(filename_pattern.format(r=rep, c=cond), delimiter="\t")
            res = res.set_index("name")
            if count_exact:
                res_info = guide_info.join(res, how="left").reset_index()
                this_edit = res_info.loc[
                    (
                        (res_info["Target base position in gRNA"] == res_info.pos + 1)
                        & (res_info.ref_base == "A")
                        & (res_info.sample_base == "G")
                    )
                ][["name", "count"]].set_index("name", drop=True)["count"]
            else:
                this_edit = (
                    res.loc[(res.ref_base == "A") & (res.sample_base == "G"), :]
                    .groupby("name")["count"]
                    .sum()
                )
            this_edit.name = "{r}_{c}".format(r=rep, c=cond)
            edits = edits.join(this_edit, how="left")
        else:
            edits = edits.fillna(0)
            return edits


class ReporterScreen(Screen):
    def __init__(self, X=None, X_edit=None, X_bcmatch=None, *args, **kwargs):
        (super().__init__)(X, *args, **kwargs)
        if (X_edit is None) != (X_bcmatch is None):
            raise ValueError("Only one of number of edits or barcode matched guide counts is specified.")
        if not X_edit is None:
            self.layers["edits"] = X_edit
        if not X_bcmatch is None:
            self.layers["X_bcmatch"] = X_bcmatch
        # if "edit_counts" in self.uns.keys():
        #     if "edit" in self.uns["edit_counts"].columns:
        #         self.uns["edit_counts"].edit = self.uns["edit_counts"].edit.map(lambda s: Edit.from_str(s))
        # if "allele_counts" in self.uns.keys():
        #     if "allele" in self.uns["allele_counts"].columns:
        #         self.uns["allele_counts"].allele = self.uns["allele_counts"].allele.map(lambda s: Allele.from_str(s))
        # if "guide_reporter_allele_counts" in self.uns.keys():
        #     self.uns["guide_reporter_allele_counts"].reporter_allele = \
        #         self.uns["guide_reporter_allele_counts"].reporter_allele.map(lambda s: Allele.from_str(s))            
        #     self.uns["guide_reporter_allele_counts"].guide_allele = \
        #         self.uns["guide_reporter_allele_counts"].guide_allele.map(lambda s: Allele.from_str(s))
        
    @classmethod
    def from_file_paths(
        cls,
        reps: List[str] = None,
        conditions: List[str] = None,
        guide_info_file_name: str = None,
        guide_count_filenames: str = None,
        guide_bcmatched_count_filenames: str = None,
        edit_count_filenames: str = None,
    ):
        guide_info = pd.read_csv(guide_info_file_name).set_index("name")
        guides = pd.DataFrame(index=(pd.read_csv(guide_info_file_name)["name"]))
        guides_lenient = _get_counts(guide_count_filenames, guides, reps, conditions)
        edits_ag = _get_edits(
            edit_count_filenames, guide_info, reps, conditions, count_exact=False
        )
        edits_exact = _get_edits(edit_count_filenames, guide_info, reps, conditions)
        if guide_bcmatched_count_filenames is not None:
            guides_bcmatch = _get_counts(
                guide_bcmatched_count_filenames, guides, reps, conditions
            )
            repscreen = cls(
                guides_lenient, edits_exact, X_bcmatch=guides_bcmatch, guides=guide_info
            )
        else:
            repscreen = cls(guides_lenient, edits_exact, guides=guide_info)
        repscreen.layers["edits_ag"] = edits_ag
        repscreen.condit["replicate"] = np.repeat(reps, len(conditions))
        repscreen.condit["sort"] = np.tile(conditions, len(reps))
        repscreen.uns["replicates"] = reps
        repscreen.uns["conditions"] = conditions
        return repscreen

    @classmethod
    def from_adata(cls, adata):
        if "X_bcmatch" in adata.layers:
            X_bcmatch = adata.layers["X_bcmatch"]
        else:
            X_bcmatch = None
        repscreen = cls(
            (adata.X),
            (adata.layers["edits"]),
            X_bcmatch=X_bcmatch,
            guides=(adata.obs),
            condit=(adata.var),
            uns=(adata.uns),
            layers=(adata.layers)
        )
        return repscreen

    def copy(self):
        adata = super().copy()
        return(type(self).from_adata(adata))

    def __add__(self, other):
        if all(self.guides.index == other.guides.index):
            if all(self.condit.index == other.condit.index):
                added_uns = dict()
                for k in self.uns.keys():
                    if k not in other.uns.keys(): continue
                    if k == "guide_edit_counts":
                        index_pair = ["guide", "edit"]
                    elif k == "edit_counts":
                        index_pair = ["guide", "edit"]
                    elif k == "allele_counts":
                        index_pair = ["guide", "allele"]
                    elif k == "guide_reporter_allele_counts":
                        index_pair = ["guide", "reporter_allele", "guide_allele"]
                    else:
                        continue
                    self_df = self.uns[k].set_index(index_pair, drop = True)
                    add_df = other.uns[k].set_index(index_pair, drop = True)
                    add_df = add_df.loc[:,self_df.columns]
                    added_uns[k] = self_df.add(add_df, fill_value = 0).astype(int).reset_index()

                if "X_bcmatch" in self.layers and "X_bcmatch" in other.layers:
                
                    added = ReporterScreen(
                        (self.X + other.X),
                        (self.layers["edits"] + other.layers["edits"]),
                        (self.layers["X_bcmatch"] + other.layers["X_bcmatch"]),
                        guides=(self.guides),
                        condit=(self.condit),
                        uns = added_uns
                    )
                else:
                    added = ReporterScreen(
                        (self.X + other.X),
                        (self.layers["edits"] + self.layers["edits"]),
                        guides=(self.guides),
                        condit=(self.condit),
                        uns = added_uns
                    )
                return added
        raise ValueError("Guides/sample description mismatch")

    def __getitem__(self, index):
        ''' TODO: currently the condit names are in ['index'] column. Making it to be the idnex will 
        allow the subsetting by condition names.
        '''
        oidx, vidx = self._normalize_indices(index)
        if isinstance(oidx, (int, np.integer)):
            if not (-self.n_obs <= oidx < self.n_obs):
                raise IndexError(f"Observation index `{oidx}` is out of range.")
            oidx += self.n_obs * (oidx < 0)
            oidx = slice(oidx, oidx + 1, 1)
        if isinstance(vidx, (int, np.integer)):
            if not (-self.n_vars <= vidx < self.n_vars):
                raise IndexError(f"Variable index `{vidx}` is out of range.")
            vidx += self.n_vars * (vidx < 0)
            vidx = slice(vidx, vidx + 1, 1)
        
        guides_include = self.guides.iloc[oidx]["name"]
        condit_include = self.condit.iloc[vidx]['index'].tolist()
        print(condit_include)
        adata = super().__getitem__(index)
        new_uns = deepcopy(adata.uns)
        for k, df in adata.uns.items():
            if "guide" in df.columns:
                if "allele" in df.columns: key_col = ["guide", "allele"]
                elif "edit" in df.columns: key_col = ["guide", "edit"]
                else: continue
                df_new = df.loc[df.guide.isin(guides_include), key_col + condit_include]
                new_uns[k] = df_new
        adata.uns = new_uns
        return(type(self).from_adata(adata))



    def get_edit_mat_from_uns(
        self, 
        ref_base, 
        alt_base, 
        match_target_position=True,
        rel_pos_start = 0,
        rel_pos_end = np.Inf,
        rel_pos_is_reporter = True,
        target_pos_col = "target_pos"
        ):
        edits = self.uns["edit_counts"].copy()
        if not self.layers["edits"] is None:
            old_edits = self.layers["edits"].copy()
            self.layers["edits"] = np.zeros_like(
                self.X,
            )
        else:
            old_edits = None
        edits["ref_base"] = edits.edit.map(lambda e: e.ref_base)
        edits["alt_base"] = edits.edit.map(lambda e: e.alt_base)
        edits = edits.loc[(edits.ref_base == ref_base) & (edits.alt_base == alt_base),:].reset_index()
        guide_len = self.guides.sequence.map(len)
        guide_name_to_idx = self.guides[["name"]].reset_index().set_index("name")
        edits["guide_idx"] = guide_name_to_idx.loc[edits.guide, "index"].reset_index(drop=True)
        edits["guide_start_pos"] = 32 - 6 - guide_len[edits.guide_idx].reset_index(drop=True)
        if not match_target_position:
            edits["rel_pos"] = edits.edit.map(lambda e: e.rel_pos)
            edits["in_rel_pos_range"] = (edits.rel_pos >= rel_pos_start + edits.guide_start_pos) & (edits.rel_pos < rel_pos_end + edits.guide_start_pos)
            good_edits = edits.loc[edits.in_rel_pos_range, ["guide", "edit"] + self.condit["index"].tolist()]
        else:
            edits["rel_pos"] = edits.edit.map(lambda e: e.rel_pos)
            edits["target_pos_guides"] = self.guides.loc[edits.guide_idx, target_pos_col].reset_index(drop=True)
            edits["target_pos_matches"] = edits.rel_pos == edits.target_pos_guides
            good_edits = edits.loc[edits.target_pos_matches, ["guide", "edit"] + self.condit["index"].tolist()]
        
        good_guide_idx = guide_name_to_idx.loc[good_edits.guide, "index"].astype(int)
        try:
            for gidx, eidx in zip(good_guide_idx, good_edits.index):
                self.layers["edits"][gidx, :] += good_edits.loc[eidx, self.condit["index"].tolist()].astype(int)
        except:
            return((good_guide_idx, good_edits))
        print("New edit matrix saved in .layers['edits']. Returning old edits.")
        return old_edits

    def get_edit_rate(
        self, 
        normalize_by_editable_base = False, 
        edited_base = None,
        editable_base_start = 3,
        editable_base_end = 8,
        bcmatch_thres = 1,
        prior_weight: float = None,
        return_result = False,
        count_layer = "X_bcmatch",
        edit_layer = "edits"
        ):
        """
        prior_weight: 
        Considering the edit rate to have prior of beta distribution with mean 0.5,
        prior weight to use when calculating posterior edit rate.
        """
        if self.layers[count_layer] is not None and self.layers[edit_layer] is not None:
            num_targetable_sites = 1.0
            if normalize_by_editable_base:
                if not edited_base in ["A", "C", "T", "G"]: 
                    raise ValueError("Specify the correct edited_base")
                num_targetable_sites = self.guides.sequence.map(
                    lambda s: s[editable_base_start:editable_base_end].count(edited_base))
            bulk_idx = np.where(
                self.condit.reset_index()["index"].map(lambda s: "bulk" in s)
            )[0]

            if prior_weight is None:
                prior_weight = 1
            n_edits = self.layers[edit_layer][:, bulk_idx].sum(axis=1)
            n_counts = self.layers[count_layer][:, bulk_idx].sum(axis=1)
            self.guides["edit_rate"] = (n_edits + prior_weight / 2) / \
                ((n_counts*num_targetable_sites) + prior_weight/2)
            self.guides["edit_rate"][n_counts < bcmatch_thres] = np.nan
            if normalize_by_editable_base:
                self.guides["edit_rate"][num_targetable_sites == 0] = np.nan
            if return_result:
                return(n_edits, n_counts)
        else:
            raise ValueError("edits or barcode matched guide counts not available.")

    def get_edit_from_allele(self):
        df = self.uns["allele_counts"].copy()
        df["edits"] = df.allele.map(lambda a: str(a).split(","))
        df = df.explode("edits").groupby(["guide", "edits"]).sum()
        df = df.reset_index().rename(columns={"edits":"edit"})
        df["edit"] = df.edit.map(lambda s: Edit.from_str(s))
        return(df)

    def _get_allele_norm(self, allele_count_df = None, thres = 10):
        '''
        Get bcmatched count to normalize allele counts
        thres: normalizing count below this number would be disregarded.
        '''
        if allele_count_df is None:
            allele_count_df = self.uns["allele_counts"]
        guide_to_idx = self.guides[["name"]].reset_index().set_index("name")
        guide_idx = guide_to_idx["index"][allele_count_df.guide]
        norm_counts = self.layers["X_bcmatch"][guide_idx.values.astype(int), :]
        norm_counts[norm_counts < thres] = np.nan
        return(norm_counts)

    def get_normalized_allele_counts(self, allele_count_df = None, norm_thres = 10):
        if allele_count_df is None:
            allele_count_df = self.uns["allele_counts"]
        alleles_norm = allele_count_df.copy()
        count_columns = self.condit["index"]
        norm_counts = self._get_allele_norm(allele_count_df = allele_count_df, thres = norm_thres)
        alleles_norm.loc[:,count_columns] = alleles_norm.loc[:,count_columns]/norm_counts
        return(alleles_norm)

    def filter_allele_counts_by_pos(self, 
        rel_pos_start = 0, rel_pos_end = 32,
        allele_uns_key = "allele_counts", filter_rel_pos = True):
        '''
        Filter alleles based on barcode matched counts, allele counts, 
        or proportion
        '''
        allele_count_df = self.uns[allele_uns_key]
        filtered_allele, filtered_edits = \
            zip(*allele_count_df.allele.map(lambda a:
                filter_allele_by_pos(a, rel_pos_start, rel_pos_end, filter_rel_pos)))
        self.uns[allele_uns_key].loc[:, "allele"] = filtered_allele
        # Hashing on Allele object messes up the order. Converting it to str and back to allele for groupby.
        allele_count_df["str_allele"] = allele_count_df.allele.map(str)
        allele_count_df = allele_count_df.groupby(["guide", "str_allele"]).sum().reset_index()
        allele_count_df.insert(1, "allele", allele_count_df.str_allele.map(lambda s: Allele.from_str(s)))
        allele_count_df.drop("str_allele", axis=1)
        print("{} edits filtered from {} alleles.".format(sum(filtered_edits), len(filtered_edits)))
        return(allele_count_df)
        

    def collapse_allele_by_target(self, ref_base, alt_base, rel_target_pos_column = "target_pos"):
        if not rel_target_pos_column in self.guides.columns:
            raise ValueError("The .guides have to have 'target_pos' specifying the relative position of target edit.")
        df = self.uns["allele_counts"].copy().reset_index(drop=True)
        df["target_pos"] = self.guides.set_index("name").loc[df.guide, rel_target_pos_column].reset_index(drop=True)
        df["has_target"] = df.apply(lambda row: row.allele.has_edit(ref_base, alt_base, rel_pos = row.target_pos), axis = 1)
        df["has_nontarget"] = df.apply(lambda row: row.allele.has_other_edit(ref_base, alt_base, rel_pos = row.target_pos), axis = 1)
        df = df.drop("target_pos", axis=1)
        res = df.groupby(["guide", "has_target", "has_nontarget"]).sum()
        return(res)


    def translate_alleles(self):
        
        if self.uns["allele_counts"] is None:
            print("No allele information. Run crisrpep-count with -a option.")
            return
        self.uns["allele_counts"]["aa_allele"] = self.uns["allele_counts"].allele.map(
            lambda s: get_aa_alleles(s)
        )

    def log_norm(self):
        super().log_norm()
        super().log_norm(output_layer = "lognorm_edits", read_count_layer="edits")

    def log_fold_changes(self, cond1, cond2, return_result=False):
        guide_fc = super().log_fold_change(cond1, cond2, return_result=return_result)
        edit_fc = super().log_fold_change(
            cond1,
            cond2,
            lognorm_counts_key="lognorm_edits",
            out_guides_suffix="edit_lfc",
            return_result=return_result,
        )
        if return_result:
            return (guide_fc, edit_fc)

    def log_fold_change_aggregates(
        self,
        cond1,
        cond2,
        aggregate_condit="replicate",
        compare_condit="sort",
        aggregate_fn="median",
        return_result=False,
        keep_per_replicate=False,
    ):
        guide_fc_agg = super().log_fold_change_aggregate(
            cond1,
            cond2,
            lognorm_counts_key="lognorm_counts",
            aggregate_condit=aggregate_condit,
            compare_condit=compare_condit,
            out_guides_suffix="lfc",
            aggregate_fn=aggregate_fn,
            return_result=return_result,
            keep_per_replicate=keep_per_replicate,
        )
        edit_fc_agg = super().log_fold_change_aggregate(
            cond1,
            cond2,
            lognorm_counts_key="lognorm_edits",
            aggregate_condit=aggregate_condit,
            compare_condit=compare_condit,
            out_guides_suffix="edit_lfc",
            aggregate_fn=aggregate_fn,
            return_result=return_result,
            keep_per_replicate=keep_per_replicate,
        )
        if return_result:
            return (guide_fc_agg, edit_fc_agg)

    def allele_log_fold_change_aggregate(
        self,
        cond1,
        cond2,
        aggregate_condit="replicate",
        compare_condit="sort",
        aggregate_fn="median",
        return_result=False,
        keep_per_replicate=False,
    ):
        pass

    # def write(self, out_path):
    #     """
    #     Write .h5ad
    #     """
    #     adata = self.copy()
    #     for k in adata.uns.keys():
    #         if 'edit' in adata.uns[k].columns:
    #             adata.uns[k].edit = adata.uns[k].edit.map(str)
    #         for c in [colname for colname in adata.uns[k].columns if "allele" in colname]:
    #             adata.uns[k].loc[:,c] = adata.uns[k][c].map(str)
    #     super(ReporterScreen, adata).write(out_path)




def concat(screens: Collection[ReporterScreen], *args, axis = 1, **kwargs):
    # TODO: var/obs info not concated if doesn't overlap
    # TODO: concat 2 times: allele_counts merging doesn't work?
    if axis == 1:
        if not all(screen.guides.index.equals(screens[0].guides.index) for screen in screens):
            raise ValueError("Guide index doesn't match.")

    adata = ad.concat(screens, *args, axis = axis, **kwargs)
    adata.obs = screens[0].guides
    keys = set()

    # Merging uns
    if axis == 1:
        for screen in screens:
            keys.update(screen.uns.keys())

        for k in keys:
            if k == "guide_edit_counts":
                merge_on = ["guide", "edit"]
            elif k == "edit_counts":
                merge_on = ["guide", "edit"]
            elif k == "allele_counts":
                merge_on = ["guide", "allele"]
            elif k == "guide_reporter_allele_counts":
                merge_on = ["guide", "reporter_allele", "guide_allele"]
            else:
                print("uns '{}' ignored during concat.".format(k))
                continue
            for i, screen in enumerate(screens):
                if i == 0:
                    adata.uns[k] = screen.uns[k]
                else:
                    try:
                        adata.uns[k] = pd.merge(
                            adata.uns[k], screen.uns[k], on=merge_on, how="outer"
                        )
                    except:
                        print("Error occurred concatenating '{}' for {}:".format(k, screen))
                        if not k in screen.uns.keys(): print("Warning: '{}' not in this screen. Skipping...".format(k))
                        else:
                            print(k)
                            print(screen.uns[k])
            adata.uns[k] = adata.uns[k].fillna(0)
            float_col = adata.uns[k].select_dtypes(include=["float64"])
            for col in float_col.columns.values:
                adata.uns[k][col] = adata.uns[k][col].astype("int64")

    return ReporterScreen.from_adata(adata)

def read_h5ad(filename):
    adata = ad.read_h5ad(filename)
    return ReporterScreen.from_adata(adata)

