import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.tsa.vector_ar.var_model import VAR


def main():
    _left = [-4.6097261953412385, 0.36803839099457036, -2.3240113911468114, 0.6098348718318576, -1.39782506791029, 4.133510766269701, -1.8415295888375525, 3.261136085002043, 2.013894564773196, -3.487127094401316, -1.4213880007799276, -2.9863970254223204, 1.6991005723630366, 3.067525211229224, 1.8811229362805193, 2.928405738286969, -2.9690408032129056, -5.708114897645281, 4.157504326904196, -0.4721753526500123, -2.627663790427377, 0.08773962024166337, 3.067421349336243, -2.2172850166845706, 1.653825434591667, -1.9183571519984923, -0.33519611851281894, 2.8640967931921395, -1.6223947691050284, -4.312219043727005, 11.40489547666457, -5.346511561800327, 1.3517382413087944, 0.10388548057258618, -1.493413059508697, -3.8675388389206873, -0.794140025865886, -0.8906514704546353, -1.1513936195650132, 2.7444289280392553, -1.1460112890000005, 6.089052287581694, 2.1113112695016696, -2.1018486807851957, 4.769260777162998, -5.0454717318925475, 3.1072051742474116, -5.272324768184276, -2.623555547826683, -1.4445764922233764, 4.707826596896425, -1.6347450952851474, 1.8021259198691766, -1.8994276369583005, -2.4753151508159092, 4.157623601587872, 4.900096757133596, 1.030241111565175, -9.934396115055364, 5.408854635692876, -2.036332337048691, -3.063575730346699, 2.2830913484548745, -2.6717599118463937, 1.1408805662431263, -1.7441918280604725, 5.392354597876228, 0.977424590191859, -0.07110747854515864, 0.7723743359215405, -3.4023368951176707, 3.2075163398692865, 1.338643790849673, -3.1307189542483655, 4.314133986928098, -6.701982853487394, -2.3705145401530885, 1.9182004089979543, 1.6417177914110397, 2.5794213498270295, -4.5496852938504375, 1.801470588235297, 6.629901960784309, -7.809983860619916, 5.142956834315662, -2.9556961094014724, -2.724213506217332, -0.11098819467824939, 1.1388888888888857, 2.9647403523960776, 2.49250421277776, -2.893158341968288, -1.9139524058862207, 0.6764217435804198, 1.6901928051133268, 1.9717589699004492, -3.0530498195736584, 2.3575042535744757]
    _right = [-1.7351859419697604, -2.7444152637216206, 0.4661153045879871, 2.8704536166734727, -2.135355388696624, 0.11741391841555071, -1.6372549019607732, 3.34502064686383, 2.755918645008336, -1.26793395836782, 1.857762640351865, -7.036388928953954, 6.029832447895373, -2.4399224307568943, -1.356625663347927, 0.6919301946360861, 2.9861743805273164, 2.571078431372541, -3.8067810457516345, 5.757761437908499, -10.318627450980387, 0.9644604504381675, 5.405803024111165, -3.569201383046078, 0.24308777454587016, 1.3623043823168786, -4.901552287581694, 3.4775396911250454, 2.6911694592017454, -4.8473357863177995, 9.456125460749078, -3.9840621168778085, 2.3105843890478184, -3.813171423725464, 7.0750476161834825, -3.680572597137015, 1.0473324002962983, -7.978860760635698, 1.9927610059341418, 0.4562551103843049, 2.459851885316212, 0.7450154050177531, 3.9399049055666353, -6.813270918638537, 1.4262664233381628, -1.3306781880440468, 2.438653815309351, -1.0744652063417135, -0.9429543275745118, -2.6279792390899104, 11.625694498417587, -3.8515475799364935, 1.1326015477765878, -4.14444795300534, -7.1803330275033375, 3.619192456529788, 2.460496959081169, 0.5361667347772823, -4.772056764971595, 6.111247610905394, 0.9328023710505846, -0.3496053981690892, -2.762191660253791, -1.273785445224533, 1.8659391027643721, 0.9713935431140186, 6.121779269533718, -1.8880718954248294, 1.3047385620915009, 2.4052287581699403, -7.285095647689445, 6.596778654225389, 1.482380007217614, -2.3437859304369226, -4.265636560440825, -4.955450484920846, 1.1355975437443817, 0.9421036880274016, -1.8490245749916028, 6.315663351393098, -6.863506334286882, 6.7960768287699125, 2.2734746227874183, -2.142293584781683, 3.167074940496491, -6.008922100283996, -3.6272987331426236, 3.9247564711758542, -2.2909191193164276, -0.22950571728488, 0.4428281529203417, -1.760638610436672, -2.7377450980392126, 5.179626046032084, 3.324070290151198, 1.8111973845525142, -7.065145334673815, 0.32505542866645953]
    _data = [[_l, _r] for _l, _r in zip(_left, _right)]
    _df = pd.DataFrame(data=_data, columns=['left', 'right'])
    _model = VAR(_df)
    _lag_results = _model.select_order()
    _chosen_lag = 1
    _model_results = _model.fit(maxlags=1)
    _whiteness = _model_results.test_whiteness(nlags=min(2, _chosen_lag))
    _normality = _model_results.test_normality()
    _granger_left_causing = _model_results.test_causality(caused='right', causing='left')
    _granger_right_causing = _model_results.test_causality(caused='left', causing='right')
    _inst_granger_left_causing = _model_results.test_inst_causality(causing='left')
    _inst_granger_right_causing = _model_results.test_inst_causality(causing='right')


if __name__ == '__main__':
    main()