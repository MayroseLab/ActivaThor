"""Available on-target scoring functions"""
import subprocess
import os

from DeepHF.scripts import prediction_util

PATH = os.getcwd()
deephf_config = None
deephf_loaded_model = None


def deephf(target_lst: list[str]) -> list[float]:
    """
    This function use the model of deephf that was improved by Yaron Orenstein`s lab

    :param target_lst: list of targets with PAM
    :return: list of on-target scores
    """
    # take 21 nt from targets
    targets = [target[0:21] for target in target_lst]
    # get deephf scores
    scores = prediction_util.get_predictions(deephf_loaded_model, deephf_config, targets)
    return list(scores)


def ucrispr(sg_seq_list: list[str]) -> list[float]:
    """
    This function will run the uCRISPR algorithm for a list of targets and will return a list of the on-target scores
    IMPORTANT: before running you need to give the path to data tables that are part of uCRISPR
    enter this to the .sh file or to the bashrc file
    'export DATAPATH=<path to folder>/uCRISPR/RNAstructure/data_tables/'
    Also make sure the uCRISPR file inside the uCRISPR folder has exe permission

    :param sg_seq_list: list of sgrnas (with PAM)
    :return: a list of ucrispr on-target scores
    """
    # make a file with guides for inputs to ucrispr
    with open(f"{PATH}/targets.txt", "w") as f:
        for sg in sg_seq_list:
            f.write(f"{sg}\n")

    # run ucrispr in terminal
    p = subprocess.run([f"{PATH}/uCRISPR/uCRISPR", "-on", f"{PATH}/targets.txt"],
                       stdout=subprocess.PIPE)

    # pares the results
    res_lst = p.stdout.decode('utf-8').split("\n")
    # delete input file
    subprocess.run(["rm", f"{PATH}/targets.txt"])
    # return a list of the results
    return [float(i.split(" ")[1]) for i in res_lst[1:len(res_lst) - 1]]
