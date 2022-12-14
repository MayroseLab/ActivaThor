"""DeepHF imports and model loading"""
import pickle
from keras.models import load_model

import ScoringFunctions


def load_deephf():
    """
    Function to load DeepHF algorithm model
    """
    with open("DeepHF/models/parallel/config.pkl", "rb") as fp:
        config = pickle.load(fp)
    model = load_model("DeepHF/models/parallel/model")
    ScoringFunctions.deephf_loaded_model = model
    ScoringFunctions.deephf_config = config
    return ScoringFunctions.deephf
