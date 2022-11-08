from catsd.parse.Keywords import KeywordsSet


class Config:

    n_cores = 8
    max_core = 2000
    keep_input_files = True

    class ORCA:
        path = None
        implicit_solvation_type = 'cpcm'
        keywords = KeywordsSet(opt=['Opt', 'PBE0', 'RIJCOSX', 'D3BJ',
                                    'def2-SVP', 'def2/J'],
                               opt_ts=['OptTS', 'Freq', 'PBE0', 'RIJCOSX',
                                       'D3BJ', 'def2-SVP', 'def2/J'],
                               hess=['Freq', 'PBE0', 'RIJCOSX', 'D3BJ',
                                     'def2-SVP', 'def2/J'],
                               optts_block=('%geom\n'
                                            'Calc_Hess true\n'
                                            'NumFreq true\n'
                                            'Recalc_Hess 20\n'
                                            'MaxIter 999\n'
                                            'end'))

    class XTB:
        path = None
        implicit_solvation_type = 'gbsa'
        keywords = KeywordsSet()