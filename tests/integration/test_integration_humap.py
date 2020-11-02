from os import system as os_system, path as os_path, chdir as os_chdir

os_chdir(os_path.dirname(os_path.abspath(__file__)))

pythonCommand="python"
#pythonCommand="python3"
inputFile = "input_humap_test.yaml"
def test_read():
    result = os_system(pythonCommand + " ../../main_read.py --input_file_name " + inputFile)
    assert result == 0

def test_train():
    result = os_system(pythonCommand + " ../../train.py --input_file_name " + inputfile)
    assert result == 0

def test_partition():
    result = os_system(pythonCommand + " ../../partition_search_seeds.py --input_file_name " + inputfile)
    assert result == 0

def test_sample():
    result = os_system(pythonCommand + " ../../main_sample.py --input_file_name " + inputfile)
    assert result == 0

def test_postprocess():
    result = os_system(pythonCommand + " ../../main_postprocess.py --input_file_name " + inputfile)
    assert result == 0

def test_eval():
    result = os_system(pythonCommand + " ../../main_eval.py --input_file_name " + inputfile)
    assert result == 0

def test_check_results():
    res = 0
    with open("../../humap/results/res_metrics.out") as f:
        lines = f.readlines()
        for line in lines[::-1]: # Reverse order
            words = line.strip().split("=")
            if words[0] == "Prediction F1 score ":
                if float(words[1]) > 0:
                    res = 1
                    break
    assert res == 1

