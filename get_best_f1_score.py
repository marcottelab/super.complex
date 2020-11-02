from glob import glob
from argparse import ArgumentParser as argparse_ArgumentParser
parser = argparse_ArgumentParser("Input parameters")
parser.add_argument("--main_folder", default="/*", help="Input parameters file name")
args = parser.parse_args()


def get_f1_score(file_name):
    res = 0
    with open(file_name) as f:
        lines = f.readlines()
        for line in lines[::-1]:  # Reverse order
            words = line.strip().split("=")
            if words[0] == "Prediction F1 score ":
                if float(words[1]) > 0:
                    res = float(words[1])
                    break
    return res


# level1subd = './humap/*/res_metrics*'
allsubd = './humap' + args.main_folder + '*/res_metrics*'
# fname = "./humap/results_73_neg_same_size_distmetropolis/res_metrics.out"
max_f1_score = 0
max_fname = ""
all_sets = []
for fname in glob(allsubd, recursive=True):
    f1score = get_f1_score(fname)
    all_sets.append((fname, f1score))
    if f1score > max_f1_score:
        max_f1_score = f1score
        max_fname = fname
all_sets = sorted(all_sets, key = lambda x: x[1])
for item in all_sets:
    print(item[0]," ", item[1])
