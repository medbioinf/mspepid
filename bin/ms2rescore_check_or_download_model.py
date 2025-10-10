#!/usr/bin/env python

import argparse
import logging

from ms2pip.constants import MODELS
from ms2pip._utils.xgb_models import validate_requested_xgb_model


def argparse_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ms2pip_model", help="Model for MS2PIP", default="HCD", type=str)
    parser.add_argument("-model_dir", help="Directory to store/find MS2PIP model", default="./ms2pip-model", type=str)

    return parser.parse_args()

if __name__ == "__main__":
    args = argparse_setup()
    logging.basicConfig(level=logging.INFO)

    ms2pip_model = args.ms2pip_model
    model_dir = args.model_dir

    # Validate / download requested model
    if ms2pip_model in MODELS.keys():
        print(f"Checking {ms2pip_model} model")
        if "xgboost_model_files" in MODELS[ms2pip_model].keys():
            validate_requested_xgb_model(
                MODELS[ms2pip_model]["xgboost_model_files"],
                MODELS[ms2pip_model]["model_hash"],
                model_dir,
            )
    