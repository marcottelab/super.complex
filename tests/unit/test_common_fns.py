from os import path as os_path, chdir as os_chdir

os_chdir(os_path.dirname(os_path.abspath(__file__)))
from sys import path as sys_path
# insert at 1, 0 is the script path (or '' in REPL)
sys_path.insert(1, '../../functions_py3/')

from jaccard_coeff import jaccard_coeff
from testClassi import generate_freq_dict

def test_jaccard_coefficient():
    assert jaccard_coeff({1,2,3,4},{3,4,5,6}) == float(1)/3

def test_generate_freq_dict():
    X_pos_test = {'nodes': [3.0,3.0,3.0,3.0,4.0]}
    res_pos = [1.0,1.0,0,0,1.0]
    pos_sizes, pos_accs, pos_tots = generate_freq_dict(X_pos_test, res_pos, 1)
    assert pos_sizes == [3,4]
    assert pos_accs == [0.5,1]
    assert pos_tots == [4,1]
