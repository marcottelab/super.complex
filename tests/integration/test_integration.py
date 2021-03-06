from os import system as os_system, path as os_path, chdir as os_chdir

os_chdir(os_path.dirname(os_path.abspath(__file__)))

pythonCommand="python"
#pythonCommand="python3"
def test_read():
    result = os_system(pythonCommand + " ../../main1_read.py")
    assert result == 0

def test_train():
    result = os_system(pythonCommand + " ../../main2_train.py")
    assert result == 0

def test_partition():
    result = os_system(pythonCommand + " ../../main3_partition_search_seeds.py")
    assert result == 0

def test_sample():
    result = os_system(pythonCommand + " ../../main4_sample.py")
    assert result == 0

def test_postprocess():
    result = os_system(pythonCommand + " ../../main5_postprocess.py")
    assert result == 0

def test_eval():
    result = os_system(pythonCommand + " ../../main6_eval.py")
    assert result == 0

def test_check_results():
    res = 0
    with open("../../toy_network/results/res_metrics.out") as f:
        lines = f.readlines()
        for line in lines[::-1]: # Reverse order
            words = line.strip().split("=")
            if words[0] == "Prediction F1 score ":
                if float(words[1]) == 1:
                    res = 1
                    break
    assert res == 1

