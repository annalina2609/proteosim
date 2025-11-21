from proteosim.liquid_chromatography import predict_lc_retention_times

def test_predict_lc_retention_times():
    peptides = ['MATSR','YEPVAEIGVGAYGTVYK','DPHSGHFVALK',]
    expected = {
        'MATSR': 8.9,''
        'YEPVAEIGVGAYGTVYK': 42.7, 
        'DPHSGHFVALK': 29.0
        }

    actual = predict_lc_retention_times(peptides)
    assert actual == expected

from proteosim.liquid_chromatography import select_retention_time_window

def test_select_retention_time_window():
    peptide_rt_map = {'LEHSNDPFNVYIESNAWQEK': 42.6,'AYVQAR': 15.6,'VLESYR': 18.4,'SCYVVENHLAIEQPNTHLPETKPSP': 54.2}
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=0, upper_ret_time=20)

    assert selected == {'AYVQAR': 15.6,'VLESYR': 18.4}