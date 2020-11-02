from argparse import ArgumentParser as argparse_ArgumentParser
import matplotlib.pyplot as plt


def sizewise_scores(filename, figname):
    with open(filename) as f:
        rawlines = f.readlines()

    scores = {}
    for rawline in rawlines:
        words = rawline.split()
        size = len(words) - 1
        score = words[-1]
        if size not in scores:
            scores[size] = [float(score), 1]
        else:
            scores[size][0] += float(score)
            scores[size][1] += 1

    sizes = []
    counts = []
    avg_scores = []
    for size in scores:
        sizes.append(size)
        avg_scores.append(float(scores[size][0])/scores[size][1])
        counts.append(scores[size][1])

    plt.figure()
    plt.plot(sizes, avg_scores, 'b.')
    plt.title("Classifier score variation with size")
    plt.xlabel("Size")
    plt.ylabel("Classifier score")
    for i, txt in enumerate(counts):
        plt.annotate(txt, (sizes[i], avg_scores[i]))
    plt.savefig(figname)
    plt.close()


def main():
    parser = argparse_ArgumentParser("Input parameters")
    parser.add_argument("--main_folder", help="Folder of the results")

    args = parser.parse_args()
    main_folder="./humap/results_73_neg_same_size_distmetropolis/"
    #main_folder="./humap/results_73_neg_same_size_distisa/"
    #main_folder="./humap/results_73_neg_same_size_distcliques/"
    #main_folder="./humap/results_73_neg_same_size_distsearch_top_neigs/"

    if args.main_folder:
        main_folder = args.main_folder
    filename = main_folder + "res_pred.out"
    figname = main_folder + "res_sizewise_scores_pred.png"
    sizewise_scores(filename, figname)
if __name__ == '__main__':
    main()
