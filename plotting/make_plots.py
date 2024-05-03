# Make plots for SUEP analysis. Reads in hdf5 files and outputs to root files
import argparse
import getpass
import logging
import os
import subprocess

import fill_utils
import numpy as np
import plot_utils
import uproot

# Import our own functions
from CMS_corrections import (
    GNN_syst,
    higgs_reweight,
    pileup_weight,
    track_killing,
    triggerSF,
)
from hist import Hist
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description="Famous Submitter")
# Name of the output file tag
parser.add_argument("-o", "--output", type=str, help="output tag", required=True)
# Where the files are coming from. Either provide:
# 1. filepath with -f
# 2. tag and dataset for something in dataDirLocal.format(tag, dataset)
# 3. tag, dataset and xrootd = 1 for something in dataDirXRootD.format(tag, dataset)
parser.add_argument(
    "-dataset",
    "--dataset",
    type=str,
    default="QCD",
    help="dataset name",
    required=False,
)
parser.add_argument(
    "-t", "--tag", type=str, default="IronMan", help="production tag", required=False
)
parser.add_argument(
    "-f", "--file", type=str, default="", help="Use specific input file"
)
parser.add_argument(
    "-s",
    "--save",
    type=str,
    help="Use specific output directory. Overrides MIT-specific paths.",
    required=False,
)
parser.add_argument(
    "--xrootd",
    type=int,
    default=0,
    help="Local data or xrdcp from hadoop (default=False)",
)
# optional: call it with --merged = 1 to append a /merged/ to the paths in options 2 and 3
parser.add_argument("--merged", type=int, default=1, help="Use merged files")
# some info about the files, highly encouraged to specify every time
parser.add_argument("-e", "--era", type=str, help="era", required=True)
parser.add_argument("--isMC", type=int, help="Is this MC or data", required=True)
parser.add_argument("--isSignal", type=int, help="Is this signal sample or not",default=0)
parser.add_argument("--scouting", type=int, default=0, help="Is this scouting or no")
# some parameters you can toggle freely
parser.add_argument("--doInf", type=int, default=0, help="make GNN plots")
parser.add_argument(
    "--doABCD", type=int, default=0, help="make plots for each ABCD+ region"
)
parser.add_argument("--doSyst", type=int, default=0, help="make systematic plots")
parser.add_argument(
    "--predictSR", type=int, default=0, help="Predict SR using ABCD method."
)
parser.add_argument(
    "--blind", type=int, default=1, help="Blind the data (default=True)"
)
parser.add_argument(
    "--weights",
    default="None",
    help="Pass the filename of the weights, e.g. --weights weights.npy",
)
options = parser.parse_args()

###################################################################################################################
# Script Parameters
###################################################################################################################

outDir = f"/data/submit/{getpass.getuser()}/SUEP/outputs/"
if options.save is not None and options.save != "None" and options != "none":
    outDir = options.save
redirector = "root://submit50.mit.edu/"
username = getpass.getuser()
if os.path.isdir("/data/submit/cms/store/user/" + username):
    # define these if --xrootd 0
    dataDirLocal = "/data/submit//cms/store/user/{}/SUEP/{}/{}/".format(
        getpass.getuser(), options.tag, options.dataset
    )
    # and these if --xrootd 1
    dataDirXRootD = "/cms/store/user/{}/SUEP/{}/{}/".format(
        getpass.getuser(), options.tag, options.dataset
    )
elif os.path.isdir("/data/submit/" + username):
    # define these if --xrootd 0
    dataDirLocal = "/data/submit/{}/SUEP/{}/{}/".format(
#        getpass.getuser(), options.tag, options.dataset
        "akanugan", options.tag, options.dataset
    )
    # and these if --xrootd 1
    dataDirXRootD = "/{}/SUEP/{}/{}/".format(
        getpass.getuser(), options.tag, options.dataset
    )
"""
Define output plotting methods, each draws from an input_method (outputs of SUEPCoffea),
and can have its own selections, ABCD regions, and signal region.
Multiple plotting methods can be defined for the same input method, as different
selections and ABCD methods can be applied.
N.B.: Include lower and upper bounds for all ABCD regions.
"""

if options.scouting:
    config = {
        "Cluster": {
            "input_method": "CL",
            "xvar": "SUEP_S1_CL",
            "xvar_regions": [0.3, 0.34, 0.5, 2.0],
            "yvar": "SUEP_nconst_CL",
            "yvar_regions": [0, 18, 50, 1000],
            "SR": [["SUEP_S1_CL", ">=", 0.5], ["SUEP_nconst_CL", ">=", 50]],
            "selections": [["ht", ">", 560], ["ntracks", ">", 0]],
        },
       "ClusterInverted": {
            "input_method": "CL",
            "xvar": "ISR_S1_CL",
            "xvar_regions": [0.3, 0.34, 0.5, 2.0],
            "yvar": "ISR_nconst_CL",
            "yvar_regions": [0, 18, 50, 1000],
            "SR": [["SUEP_S1_CL", ">=", 0.5], ["SUEP_nconst_CL", ">=", 50]],
            "selections": [["ht", ">", 560], ["ntracks", ">", 0]],
        },
    }
else:
    config = {
        "Cluster70": {
            "input_method": "CL",
            "xvar": "SUEP_S1_CL",
            "xvar_regions": [0.3, 0.4, 0.5, 2.0],
            "yvar": "SUEP_nconst_CL",
            "yvar_regions": [30, 50, 70, 1000],
            "SR": [["SUEP_S1_CL", ">=", 0.5], ["SUEP_nconst_CL", ">=", 70]],
            "selections": [["ht_JEC", ">", 1200], ["ntracks", ">", 0]],
        },
        "ClusterInverted": {
            "input_method": "CL",
            "xvar": "ISR_S1_CL",
            "xvar_regions": [0.3, 0.4, 0.5, 2.0],
            "yvar": "ISR_nconst_CL",
            "yvar_regions": [30, 50, 70, 1000],
            # "SR": [["SUEP_S1_CL", ">=", 0.5], ["SUEP_nconst_CL", ">=", 75]],
            "SR": [["SUEP_S1_CL", ">=", 0.5], ["SUEP_nconst_CL", ">=", 70]],
            "selections": [["ht_JEC", ">", 1200], ["ntracks", ">", 0]],
        },
    }

if options.doInf:
    config.update(
        {
            "GNN": {
                "input_method": "GNN",
                "xvar": "SUEP_S1_GNN",
                "xvar_regions": [0.3, 0.4, 0.5, 1.0],
                "yvar": "single_l5_bPfcand_S1_SUEPtracks_GNN",
                "yvar_regions": [0.0, 0.5, 1.0],
                "SR": [
                    ["SUEP_S1_GNN", ">=", 0.5],
                    ["single_l5_bPfcand_S1_SUEPtracks_GNN", ">=", 0.5],
                ],
                "SR2": [
                    ["SUEP_S1_CL", ">=", 0.5],
                    ["SUEP_nconst_CL", ">=", 80],
                ],  # both are blinded
                "selections": [["ht_JEC", ">", 1200], ["ntracks", ">", 40]],
                "models": ["single_l5_bPfcand_S1_SUEPtracks"],
                "fGNNsyst": "../data/GNN/GNNsyst.json",
                "GNNsyst_bins": [0.0j, 0.25j, 0.5j, 0.75j, 1.0j],
            },
            "GNNInverted": {
                "input_method": "GNNInverted",
                "xvar": "ISR_S1_GNNInverted",
                "xvar_regions": [0.0, 1.5, 2.0],
                "yvar": "single_l5_bPfcand_S1_SUEPtracks_GNNInverted",
                "yvar_regions": [0.0, 1.5, 2.0],
                "SR": [
                    ["ISR_S1_GNNInverted", ">=", 10.0],
                    ["single_l5_bPfcand_S1_SUEPtracks_GNNInverted", ">=", 10.0],
                ],
                # 'SR2': [['ISR_S1_CL', '>=', 0.5], ['ISR_nconst_CL', '>=', 80]], # both are blinded
                "selections": [["ht_JEC", ">", 1200], ["ntracks", ">", 40]],
                "models": ["single_l5_bPfcand_S1_SUEPtracks"],
            },
        }
    )


def open_file(options, redirector, ifile):
    if not options.xrootd:
        return fill_utils.h5load(ifile, "vars")
    if os.path.exists(options.dataset + ".hdf5"):
        os.system("rm " + options.dataset + ".hdf5")
    xrd_file = redirector + ifile
    os.system(f"xrdcp -s {xrd_file} {options.dataset}.hdf5")
    return fill_utils.h5load(options.dataset + ".hdf5", "vars")


def create_output_clusterInverted(output, label, regions_list):
    output.update(
        {
            # 2D histograms
            f"2D_ISR_S1_vs_ntracks_{label}": Hist.new.Reg(
                100, 0, 1.0, name=f"ISR_S1_{label}", label="$Sph_1$"
            )
            .Reg(200, 0, 500, name=f"ntracks_{label}", label="# Tracks")
            .Weight(),
            f"2D_ISR_S1_vs_ISR_nconst_{label}": Hist.new.Reg(
                100, 0, 1.0, name=f"ISR_S1_{label}", label="$Sph_1$"
            )
            .Reg(501, 0, 500, name=f"nconst_{label}", label="# Constituents")
            .Weight(),
            f"2D_ISR_nconst_vs_ISR_pt_avg_{label}": Hist.new.Reg(
                501, 0, 500, name=f"ISR_nconst_{label}"
            )
            .Reg(500, 0, 500, name=f"ISR_pt_avg_{label}")
            .Weight(),
        }
    )
    # variables from the dataframe for all the events, and those in A, B, C regions
    for r in regions_list:
        output.update(
            {
                f"{r}ISR_nconst_{label}": Hist.new.Reg(
                    501,
                    0,
                    500,
                    name=f"{r}ISR_nconst_{label}",
                    label="# Tracks in ISR",
                ).Weight(),
                f"{r}ISR_S1_{label}": Hist.new.Reg(
                    100, 0, 1, name=f"{r}ISR_S1_{label}", label="$Sph_1$"
                ).Weight(),
            }
        )
        if (
            "up" not in label and "down" not in label and r == ""
        ):  # don't care for systematics for these
            output.update(
                {
                    f"ISR_pt_{label}": Hist.new.Reg(
                        100,
                        0,
                        2000,
                        name=f"ISR_pt_{label}",
                        label=r"ISR $p_T$ [GeV]",
                    ).Weight(),
                    f"ISR_pt_avg_{label}": Hist.new.Reg(
                        500,
                        0,
                        500,
                        name=f"ISR_pt_avg_{label}",
                        label=r"ISR Components $p_T$ Avg.",
                    ).Weight(),
                    f"ISR_eta_{label}": Hist.new.Reg(
                        100,
                        -5,
                        5,
                        name=f"ISR_eta_{label}",
                        label=r"ISR $\eta$",
                    ).Weight(),
                    f"ISR_phi_{label}": Hist.new.Reg(
                        100,
                        -6.5,
                        6.5,
                        name=f"ISR_phi_{label}",
                        label=r"ISR $\phi$",
                    ).Weight(),
                    f"ISR_mass_{label}": Hist.new.Reg(
                        150,
                        0,
                        4000,
                        name=f"ISR_mass_{label}",
                        label="ISR Mass [GeV]",
                    ).Weight(),
                }
            )


def create_output_GNN(abcd, output, label, regions_list):
    # 2D histograms
    for model in abcd["models"]:
        output.update(
            {
                f"2D_SUEP_S1_vs_{model}_{label}": Hist.new.Reg(
                    100, 0, 1.0, name=f"SUEP_S1_{label}", label="$Sph_1$"
                )
                .Reg(100, 0, 1, name=f"{model}_{label}", label="GNN Output")
                .Weight(),
                f"2D_SUEP_nconst_vs_{model}_{label}": Hist.new.Reg(
                    501,
                    0,
                    500,
                    name=f"SUEP_nconst_{label}",
                    label="# Const",
                )
                .Reg(100, 0, 1, name=f"{model}_{label}", label="GNN Output")
                .Weight(),
            }
        )

    output.update(
        {
            f"2D_SUEP_nconst_vs_SUEP_S1_{label}": Hist.new.Reg(
                501, 0, 500, name=f"SUEP_nconst_{label}", label="# Const"
            )
            .Reg(100, 0, 1, name=f"SUEP_S1_{label}", label="$Sph_1$")
            .Weight(),
        }
    )

    for r in regions_list:
        output.update(
            {
                f"{r}SUEP_nconst_{label}": Hist.new.Reg(
                    501,
                    0,
                    500,
                    name=f"{r}SUEP_nconst{label}",
                    label="# Constituents",
                ).Weight(),
                f"{r}SUEP_S1_{label}": Hist.new.Reg(
                    100,
                    -1,
                    2,
                    name=f"{r}SUEP_S1_{label}",
                    label="$Sph_1$",
                ).Weight(),
            }
        )
        for model in abcd["models"]:
            output.update(
                {
                    f"{r}{model}_{label}": Hist.new.Reg(
                        100,
                        0,
                        1,
                        name=f"{r}{model}_{label}",
                        label="GNN Output",
                    ).Weight()
                }
            )


def create_output_GNNInverted(abcd, output, label, regions_list):
    # 2D histograms
    for model in abcd["models"]:
        output.update(
            {
                f"2D_ISR_S1_vs_{model}_{label}": Hist.new.Reg(
                    100, 0, 1.0, name=f"ISR_S1_{label}", label="$Sph_1$"
                )
                .Reg(100, 0, 1, name=f"{model}_{label}", label="GNN Output")
                .Weight(),
                f"2D_ISR_nconst_vs_{model}_{label}": Hist.new.Reg(
                    501, 0, 500, name=f"ISR_nconst_{label}", label="# Const"
                )
                .Reg(100, 0, 1, name=f"{model}_{label}", label="GNN Output")
                .Weight(),
            }
        )
    output.update(
        {
            f"2D_ISR_nconst_vs_ISR_S1_{label}": Hist.new.Reg(
                501, 0, 500, name=f"ISR_nconst_{label}", label="# Const"
            )
            .Reg(100, 0, 1, name=f"ISR_S1_{label}", label="$Sph_1$")
            .Weight()
        }
    )

    for r in regions_list:
        output.update(
            {
                f"{r}ISR_nconst_{label}": Hist.new.Reg(
                    501,
                    0,
                    500,
                    name=f"{r}ISR_nconst{label}",
                    label="# Tracks in ISR",
                ).Weight(),
                f"{r}ISR_S1_{label}": Hist.new.Reg(
                    100, -1, 2, name=f"{r}ISR_S1_{label}", label="$Sph_1$"
                ).Weight(),
            }
        )
        for model in abcd["models"]:
            output.update(
                {
                    f"{r}{model}_{label}": Hist.new.Reg(
                        100,
                        0,
                        1,
                        name=f"{r}{model}_{label}",
                        label="GNN Output",
                    ).Weight()
                }
            )


# output histos
def create_output_file(label, abcd, options):
    # don't recreate histograms if called multiple times with the same output label
    if label in output["labels"]:
        return output
    else:
        output["labels"].append(label)

    # ABCD histogram
    xvar = abcd["xvar"]
    yvar = abcd["yvar"]
    xvar_regions = abcd["xvar_regions"]
    yvar_regions = abcd["yvar_regions"]
    output.update(
        {
            f"ABCDvars_{label}": Hist.new.Reg(
                100, yvar_regions[0], yvar_regions[-1], name=xvar
            )
            .Reg(100, xvar_regions[0], xvar_regions[-1], name=yvar)
            .Weight()
        }
    )

    # defnie all the regions, will be used to make historgams for each region
    regions = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    n_regions = (len(xvar_regions) - 1) * (len(yvar_regions) - 1)
    regions_list = [""]
    if options.doABCD:
        regions_list += [regions[i] + "_" for i in range(n_regions)]

    ###########################################################################################################################
    # variables from the dataframe for all the events, and those in A, B, C regions
    runList_2018 = [315259, 315264, 315265, 315267, 315270, 315322, 315339, 315357, 315361, 315363, 315366, 315420, 315488, 315489, 315490, 315506, 315510, 315512, 315543, 315555, 315556, 315557, 315640, 315721, 315641, 315642, 315644, 315645, 315646, 315647, 315648, 315689, 315690, 315702, 315703, 315704, 315705, 315713, 315741, 315764, 315770, 315784, 315785, 315786, 315790, 315800, 315801, 315840, 315973, 316061, 315974, 316058, 316059, 316060, 316062, 316082, 316110, 316111, 316113, 316114, 316153, 316186, 316187, 316199, 316200, 316201, 316202, 316216, 316217, 316218, 316219, 316239, 316240, 316241, 316271, 316361, 316362, 316363, 316377, 316378, 316379, 316380, 316455, 316457, 316469, 316470, 316472, 316505, 316569, 316590, 316613, 316615, 316666, 316667, 316700, 316701, 316702, 316715, 316716, 316717, 316718, 316719, 316720, 316721, 316722, 316723, 316758, 316766, 316876, 316877, 316879, 316928, 316985, 316993, 316994, 316995, 317080, 317087, 317089, 317182, 317212, 317213, 317279, 317291, 317292, 317297, 317319, 317320, 317338, 317339, 317340, 317382, 317383, 317391, 317392, 317435, 317438, 317475, 317478, 317484, 317488, 317527, 317591, 317626, 317640, 317641, 317648, 317649, 317650, 317661, 317663, 317683, 317696, 318733, 318828, 318872, 318874, 318876, 318877, 319077, 319337, 319347, 319348, 319349, 319449, 319450, 319456, 319459, 319486, 319503, 319524, 319526, 319528, 319579, 319625, 319639, 319656, 319657, 319658, 319659, 319678, 319687, 319697, 319698, 319756, 319840, 319841, 319847, 319848, 319849, 319851, 319853, 319854, 319908, 319909, 319910, 319912, 319913, 319914, 319915, 319941, 319942, 319950, 319991, 319992, 319993, 320002, 320006, 320010, 320011, 320012, 320023, 320024, 320025, 320026, 320038, 320039, 320040, 320059, 320060, 320061, 320062, 320063, 320064, 320065, 320673, 320674, 320688, 320712, 320757, 320804, 320807, 320809, 320821, 320822, 320823, 320824, 320838, 320840, 320841, 320853, 320854, 320855, 320856, 320857, 320858, 320859, 320887, 320888, 320916, 320917, 320920, 320933, 320934, 320936, 320941, 320980, 320995, 321012, 320996, 321004, 321005, 321006, 321007, 321009, 321010, 321011, 321051, 321055, 321067, 321068, 321069, 321119, 321121, 321122, 321124, 321126, 321134, 321138, 321140, 321149, 321165, 321166, 321167, 321177, 321178, 321218, 321219, 321221, 321230, 321231, 321232, 321233, 321262, 321283, 321294, 321295, 321296, 321305, 321393, 321311, 321312, 321313, 321396, 321397, 321414, 321415, 321431, 321432, 321433, 321434, 321436, 321457, 321461, 321475, 321710, 321712, 321730, 321732, 321735, 321755, 321758, 321760, 321773, 321774, 321775, 321776, 321777, 321778, 321780, 321781, 321813, 321815, 321817, 321818, 321820, 321831, 321832, 321833, 321834, 321879, 321880, 321887, 321908, 321909, 321917, 321919, 321933, 321960, 321961, 321973, 321975, 321988, 321990, 322013, 322014, 322022, 322040, 322057, 322068, 322079, 322106, 322113, 322118, 322179, 322201, 322204, 322625, 322222, 322252, 322317, 322319, 322322, 322324, 322332, 322348, 322355, 322356, 322381, 322407, 322430, 322431, 322480, 322492, 322510, 322599, 322602, 322603, 322605, 322617, 322633, 323414, 323423, 323470, 323471, 323472, 323473, 323474, 323475, 323487, 323488, 323492, 323493, 323495, 323524, 323525, 323526, 323693, 323696, 323702, 323725, 323726, 323727, 323755, 323775, 323778, 323790, 323794, 323841, 323857, 323940, 323954, 323976, 323978, 323980, 323983, 323997, 324021, 324022, 324077, 324201, 324202, 324205, 324206, 324207, 324209, 324237, 324245, 324293, 324315, 324318, 324420, 324729, 324747, 324764, 324765, 324769, 324772, 324785, 324791, 324835, 324840, 324841, 324846, 324878, 324897, 324970, 324980, 324997, 324998, 324999, 325000, 325001, 325022, 325057, 325097, 325098, 325099, 325100, 325101, 325110, 325117, 325159, 325168, 325169, 325170, 325172]
    output.update(
        {
            f"ht_{label}": Hist.new.Reg(
                50, 50, 3500, name=f"ht_{label}", label="HT"
            ).Weight(),
            f"ht_JEC_{label}": Hist.new.Reg(
                50, 50, 3500, name=f"ht_JEC_{label}", label="HT JEC"
            ).Weight(),
            f"ht_JEC_JER_up_{label}": Hist.new.Reg(
                200,
                0,
                7000,
                name=f"ht_JEC_JER_up_{label}",
                label="HT JEC up",
            ).Weight(),
            f"ht_JEC_JER_down_{label}": Hist.new.Reg(
                200,
                0,
                7000,
                name=f"ht_JEC_JER_down_{label}",
                label="HT JEC JER down",
            ).Weight(),
            f"ht_JEC_JES_up_{label}": Hist.new.Reg(
                200,
                0,
                7000,
                name=f"ht_JEC_JES_up_{label}",
                label="HT JEC JES up",
            ).Weight(),
            f"ht_JEC_JES_down_{label}": Hist.new.Reg(
                200,
                0,
                7000,
                name=f"ht_JEC_JES_down_{label}",
                label="HT JEC JES down",
            ).Weight(),
            f"ntracks_{label}": Hist.new.Reg(
                101,
                0,
                500,
                name=f"ntracks_{label}",
                label="# Tracks in Event",
            ).Weight(),
            f"ngood_fastjets_{label}": Hist.new.Reg(
                9,
                0,
                10,
                name=f"ngood_fastjets_{label}",
                label="# FastJets in Event",
            ).Weight(),
            f"PV_npvs_{label}": Hist.new.Reg(
                199,
                0,
                200,
                name=f"PV_npvs_{label}",
                label="# PVs in Event ",
            ).Weight(),
            f"Pileup_nTrueInt_{label}": Hist.new.Reg(
                199,
                0,
                200,
                name=f"Pileup_nTrueInt_{label}",
                label="# True Interactions in Event ",
            ).Weight(),
            f"ngood_ak4jets_{label}": Hist.new.Reg(
                19,
                0,
                20,
                name=f"ngood_ak4jets_{label}",
                label="# ak4jets in Event",
            ).Weight(),
            f"run_{label}": Hist.new.Reg(
                #runList_2018,
                9914,
                315258.5,
                325172.5,
                #9411,
                #297049.5,
                #306460.5,
                name=f"run_{label}",
                label="# Events in run",
            ).Weight(),
            f"lumi_{label}": Hist.new.Reg(
                2500,
                0,
                2500,
                name=f"lumi_{label}",
                label="lumisection",
            ).Weight(),
            
        }
    )

    ###########################################################################################################################
    if any([lbl in label for lbl in ["ISRRemoval", "Cluster", "Cone"]]):
        # 2D histograms
        output.update(
            {
                f"2D_SUEP_S1_vs_SUEP_nconst_{label}": Hist.new.Reg(
                    100, 0, 1.0, name=f"SUEP_S1_{label}", label="$Sph_1$"
                )
                .Reg(501, 0, 500, name=f"nconst_{label}", label="# Constituents")
                .Weight(),
            }
        )
        if (
            "up" not in label and "down" not in label
        ):  # don't care for systematics for these
            output.update(
                {
                    f"2D_SUEP_S1_vs_ntracks_{label}": Hist.new.Reg(
                        100, 0, 1.0, name=f"SUEP_S1_{label}", label="$Sph_1$"
                    )
                    .Reg(100, 0, 500, name=f"ntracks_{label}", label="# Tracks")
                    .Weight(),
                    f"2D_SUEP_nconst_vs_SUEP_pt_avg_{label}": Hist.new.Reg(
                        501, 0, 500, name=f"SUEP_nconst_{label}", label="# Const"
                    )
                    .Reg(200, 0, 500, name=f"SUEP_pt_avg_{label}", label="$p_T Avg$")
                    .Weight(),
                    f"2D_SUEP_eta_vs_SUEP_nconst_{label}": Hist.new.Reg(
                        100, -5, 5, name=f"SUEP_eta_{label}", label=r"$\eta$"
                    )
                    .Reg(501, 0, 500, name=f"nconst_{label}", label="# Constituents")
                    .Weight(),
                    f"2D_SUEP_pt_vs_SUEP_nconst_{label}": Hist.new.Reg(
                        100, 0, 2000, name=f"SUEP_pt_{label}", label="SUEP $p_T$"
                    )
                    .Reg(501, 0, 500, name=f"nconst_{label}", label="# Constituents")
                    .Weight(),
                }
            )

        # variables from the dataframe for all the events, and those in A, B, C regions
        for r in regions_list:
            output.update(
                {
                    f"{r}SUEP_nconst_{label}": Hist.new.Reg(
                        501,
                        0,
                        500,
                        name=f"{r}SUEP_nconst_{label}",
                        label="# Constituents",
                    ).Weight(),
                    f"{r}SUEP_S1_{label}": Hist.new.Reg(
                        100, 0, 1, name=f"{r}SUEP_S1_{label}", label="$Sph_1$"
                    ).Weight(),
                }
            )
            if (
                r == "" and "up" not in label and "down" not in label
            ):  # don't care for systematics for these
                output.update(
                    {
                        f"{r}SUEP_genMass_{label}": Hist.new.Reg(
                            100,
                            0,
                            1200,
                            name=f"{r}SUEP_genMass_{label}",
                            label="Gen Mass of SUEP ($m_S$) [GeV]",
                        ).Weight(),
                        f"{r}SUEP_pt_{label}": Hist.new.Reg(
                            100,
                            0,
                            2000,
                            name=f"{r}SUEP_pt_{label}",
                            label=r"SUEP $p_T$ [GeV]",
                        ).Weight(),
                        f"{r}SUEP_delta_pt_genPt_{label}": Hist.new.Reg(
                            400,
                            -2000,
                            2000,
                            name=f"{r}SUEP_delta_pt_genPt_{label}",
                            label="SUEP $p_T$ - genSUEP $p_T$ [GeV]",
                        ).Weight(),
                        f"{r}SUEP_pt_avg_{label}": Hist.new.Reg(
                            200,
                            0,
                            500,
                            name=f"{r}SUEP_pt_avg_{label}",
                            label=r"SUEP Components $p_T$ Avg.",
                        ).Weight(),
                        f"{r}SUEP_eta_{label}": Hist.new.Reg(
                            100,
                            -5,
                            5,
                            name=f"{r}SUEP_eta_{label}",
                            label=r"SUEP $\eta$",
                        ).Weight(),
                        f"{r}SUEP_phi_{label}": Hist.new.Reg(
                            100,
                            -6.5,
                            6.5,
                            name=f"{r}SUEP_phi_{label}",
                            label=r"SUEP $\phi$",
                        ).Weight(),
                        f"{r}SUEP_mass_{label}": Hist.new.Reg(
                            150,
                            0,
                            2000,
                            name=f"{r}SUEP_mass_{label}",
                            label="SUEP Mass [GeV]",
                        ).Weight(),
                        f"{r}SUEP_delta_mass_genMass_{label}": Hist.new.Reg(
                            400,
                            -2000,
                            2000,
                            name=f"{r}SUEP_delta_mass_genMass_{label}",
                            label="SUEP Mass - genSUEP Mass [GeV]",
                        ).Weight(),
                    }
                )
    ###########################################################################################################################
    if "ClusterInverted" in label:
        create_output_clusterInverted(output, label, regions_list)

    ###########################################################################################################################
    if label == "GNN" and options.doInf:
        create_output_GNN(abcd, output, label, regions_list)

    ###########################################################################################################################
    if label == "GNNInverted" and options.doInf:
        create_output_GNNInverted(abcd, output, label, regions_list)

    ###########################################################################################################################

    return output


def calculate_systematic(
    df,
    config,
    syst,
    options,
):
    # prepare new event weight
    if options.isMC:
        df["event_weight"] = df["genweight"].to_numpy()
    else:
        df["event_weight"] = np.ones(df.shape[0])

    if options.isMC == 1:
        if options.scouting != 1:
            # 1) pileup weights
            puweights, puweights_up, puweights_down = pileup_weight.pileup_weight(
                options.era
            )
            pu = pileup_weight.get_pileup_weights(
                df, syst, puweights, puweights_up, puweights_down
            )
            df["event_weight"] *= pu

            # 2) TriggerSF weights
            (
                trig_bins,
                trig_weights,
                trig_weights_up,
                trig_weights_down,
            ) = triggerSF.triggerSF(options.era)
            trigSF = triggerSF.get_trigSF_weight(
                df,
                syst,
                trig_bins,
                trig_weights,
                trig_weights_up,
                trig_weights_down,
            )
            df["event_weight"] *= trigSF

            # 3) PS weights
            if "PSWeight" in syst and syst in df.keys():
                df["event_weight"] *= df[syst]

            # 3) prefire weights
            if options.era == "2016" or options.era == "2017":
                if "prefire" in syst and syst in df.keys():
                    df["event_weight"] *= df[syst]
                else:
                    df["event_weight"] *= df["prefire_nom"]

        else:
            # 1) pileup weights
            puweights, puweights_up, puweights_down = pileup_weight.pileup_weight(
                options.era
            )
            pu = pileup_weight.get_pileup_weights(
                df, syst, puweights, puweights_up, puweights_down
            )
            df["event_weight"] *= pu

            # 2) TriggerSF weights
            trigSF = triggerSF.get_scout_trigSF_weight(np.array(df["ht"]).astype(int), syst, options.era)
            df["event_weight"] *= trigSF

            # 3) PS weights
            if "PSWeight" in syst and syst in df.keys():
                df["event_weight"] *= df[syst]


        # 5) Higgs_pt weights
        if "mS125" in options.dataset:
            (
                higgs_bins,
                higgs_weights,
                higgs_weights_up,
                higgs_weights_down,
            ) = higgs_reweight.higgs_reweight(df["SUEP_genPt"])
            higgs_weight = higgs_reweight.get_higgs_weight(
                df,
                syst,
                higgs_bins,
                higgs_weights,
                higgs_weights_up,
                higgs_weights_down,
            )
            df["event_weight"] *= higgs_weight

    # 6) scaling weights
    # N.B.: these aren't part of the systematics, just an optional scaling
    if scaling_weights is not None:
        df = fill_utils.apply_scaling_weights(
            df.copy(),
            scaling_weights,
            config["Cluster"],
            regions="ABCDEFGHI",
            x_var="SUEP_S1_CL",
            y_var="SUEP_nconst_CL",
            z_var="ht",
        )

    for label_out, config_out in config.items():
        if "track_down" in label_out and syst != "":
            continue  # don't run other systematics when doing track killing systematic
        if options.isMC and syst != "":
            if any([j in label_out for j in jet_corrections]):
                continue  # don't run other systematics when doing jet systematics

        # rename if we have applied a systematic
        if len(syst) > 0:
            label_out = label_out + "_" + syst

        # initialize new hists, if needed
        output.update(create_output_file(label_out, config_out, options))

        # prepare the DataFrame for plotting: blind, selections
        df_plot = fill_utils.prepareDataFrame(
            df.copy(), config_out, label_out, isMC=options.isMC, blind=options.blind
        )

        # auto fill all histograms
        fill_utils.auto_fill(
            df_plot,
            output,
            config_out,
            label_out,
            isMC=options.isMC,
            do_abcd=options.doABCD,
        )


#############################################################################################################

# variables that will be filled
nfailed = 0
total_gensumweight = 0
output = {"labels": []}

# get list of files
username = getpass.getuser()
if options.file:
    files = [options.file]
elif options.xrootd:
    dataDir = dataDirXRootD
    if options.merged:
        dataDir += "merged/"
    result = subprocess.check_output(["xrdfs", redirector, "ls", dataDir])
    result = result.decode("utf-8")
    files = result.split("\n")
    files = [f for f in files if len(f) > 0]
else:
    dataDir = dataDirLocal
    if options.merged:
        dataDir += "merged/"
    files = [dataDir + f for f in os.listdir(dataDir)]

# get cross section
xsection = 1.0
if options.isMC:
    xsection = fill_utils.getXSection(options.dataset, options.era, SUEP=bool(options.isSignal))

# custom per region weights
scaling_weights = None
if options.weights is not None and options.weights != "None":
    scaling_weights = fill_utils.read_in_weights(options.weights)

# add new output methods for track killing and jet energy systematics
if options.isMC and options.doSyst:
    # track systematics: need to use the track_down version of the variables
    new_config_track_killing = fill_utils.get_track_killing_config(config)

    # jet systematics: just change ht to ht_SYS (e.g. ht -> ht_JEC_JES_up)
    jet_corrections = [
        "JER_up",
        "JER_down",
        "JES_up",
        "JES_down",
    ]
    new_config_jet_corrections = fill_utils.get_jet_corrections_config(
        config, jet_corrections
    )

    config = config | new_config_jet_corrections
    config = config | new_config_track_killing

logging.info("Setup ready, filling histograms now.")

# Plotting loop #######################################################################
for ifile in tqdm(files):
    #####################################################################################
    # ---- Load file
    #####################################################################################

    # get the file
    df, metadata = open_file(options, redirector, ifile)

    # check if file is corrupted
    if type(df) == int:
        nfailed += 1
        continue

    # update the gensumweight
    if options.isMC and metadata != 0:
        total_gensumweight += metadata["gensumweight"]

    # check if file is empty
    if "empty" in list(df.keys()):
        continue
    if df.shape[0] == 0:
        continue

    if options.isMC and options.doSyst:
        sys_loop = [
            "",
            "puweights_up",
            "puweights_down",
            "trigSF_up",
            "trigSF_down",
            "PSWeight_ISR_up",
            "PSWeight_ISR_down",
            "PSWeight_FSR_up",
            "PSWeight_FSR_down",
            "prefire_up",
            "prefire_down",
        ]
        if "mS125" in options.dataset:
            sys_loop += [
                "higgs_weights_up",
                "higgs_weights_down",
            ]
    else:
        sys_loop = [""]

    for syst in sys_loop:
        # prepare new event weight
        calculate_systematic(df, config, syst, options)

    #####################################################################################
    # ---- End
    #####################################################################################

    # remove file at the end of loop
    if options.xrootd:
        os.system("rm " + options.dataset + ".hdf5")

logging.warning("Number of files that failed to be read: " + str(nfailed))
# End plotting loop ###################################################################

# not needed anymore
output.pop("labels")

logging.info("Applying symmetric systematics and normalization.")

# do the systematics
if options.isMC and options.doSyst:
    # do the tracks UP systematic
    output = track_killing.generate_up_histograms(config.keys(), output)

    if options.doInf:
        # do the GNN systematic
        GNN_syst.apply_GNN_syst(
            output,
            config["GNN"]["fGNNsyst"],
            config["GNN"]["models"],
            config["GNN"]["GNNsyst_bins"],
            options.era,
            out_label="GNN",
        )

# apply normalization
if options.isMC:
    print("xsection",xsection, "total_gensumweight", total_gensumweight)
    output = fill_utils.apply_normalization(output, xsection / total_gensumweight)

# Make ABCD expected histogram for signal region

if options.doABCD and options.blind and options.predictSR:
    logging.info("Predicting SR using ABCD method.")

    # Loop through every configuration
    for label_out, config_out in config.items():
        xregions = np.array(config_out["xvar_regions"]) * 1.0j
        yregions = np.array(config_out["yvar_regions"]) * 1.0j
        xvar = config_out["xvar"].replace("_" + config_out["input_method"], "")
        yvar = config_out["yvar"].replace("_" + config_out["input_method"], "")

        hist_name = f"2D_{xvar}_vs_{yvar}_{label_out}"

        # Check if histogram exists
        if hist_name not in output.keys():
            logging.warning(f"{hist_name} has not been created.")
            continue

        # Only calculate predicted for 9 region ABCD
        if (len(xregions) != 4) or (len(yregions) != 4):
            logging.warning(
                f"Can only calculate SR for 9 region ABCD, skipping {label_out}"
            )
            continue

        # Calculate SR from ABCD method
        # sum_var = 'x' corresponds to scaling F histogram
        SR, SR_exp = plot_utils.ABCD_9regions_errorProp(
            output[hist_name], xregions, yregions, sum_var="x"
        )

        output[f"I_{yvar}_{label_out}_exp"] = SR_exp

if options.dataset:
    outFile = outDir + "/" + options.dataset + "_" + options.output
else:
    outFile = os.path.join(outDir, options.output)
logging.info("Saving outputs to " + outFile + ".")
with uproot.recreate(outFile + ".root") as froot:
    for h, hist in output.items():
        froot[h] = hist
