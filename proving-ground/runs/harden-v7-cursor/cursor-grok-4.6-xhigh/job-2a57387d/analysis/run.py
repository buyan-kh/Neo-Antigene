"""Driver: load the screen data and print the significant features."""

import csv
import os

from stats import screen

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data", "features.csv")


def load():
    features = {}
    with open(DATA, newline="") as handle:
        for row in csv.DictReader(handle):
            feature_id = row["feature_id"]
            features.setdefault(feature_id, {"A": [], "B": []})
            features[feature_id][row["group"]].append(float(row["value"]))
    return features


def main():
    features = load()
    _, significant = screen(features)
    print("{} features screened".format(len(features)))
    print("{} significant".format(len(significant)))
    for feature_id in significant:
        print(feature_id)


if __name__ == "__main__":
    main()
