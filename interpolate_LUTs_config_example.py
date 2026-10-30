# helper config file to be used when running interpolate_LUTs.py
# 
# nbins: number of eta bins to interpolate over
# num_cols: purely cosmetic; number of columns output LUTs will have before going to new line
# factor_bits: dictionary with coeffs as keys telling how many bits to shift the raw parameters when storing them in LUTs
# lut_bits: dictionary with coeffs as keys telling how many bits of precision are to be saved in LUTs
# trash_eta_zones_by_layer: dictionary indexed by layer, then each resulting 2D list contains a sequence of pairs of 
#                           eta values for which the fits should be thrown out and interpolated to instead

nbins = 512
num_cols = 32

factor_bits = {'prop1_0':  9, 'prop1_1':  9, 'prop1_2':  9, 'prop1_3':  9, 'prop1_4':  9, 'prop1_5':  9, 'prop1_6':  9, 'prop1_7':  9, 'prop1_8':  9, 'prop1_9':  9,
               'prop2_0': 23, 'prop2_1': 23, 'prop2_2': 23, 'prop2_3': 23, 'prop2_4': 23, 'prop2_5': 23, 'prop2_6': 23, 'prop2_7': 23, 'prop2_8': 23, 'prop2_9': 23,
               'prop3_0': 22, 'prop3_1': 19, 'prop3_2': 22, 'prop3_3': 19, 'prop3_4': 22, 'prop3_5': 19, 'prop3_6': 22, 'prop3_7': 19, 'prop3_8': 22, 'prop3_9': 22,
               'res0_0':   0, 'res0_1':   0, 'res0_2':   0, 'res0_3':   0, 'res0_4':   0, 'res0_5':   0, 'res0_6':   0, 'res0_7':   0, 'res0_8':   0, 'res0_9':   0,
               'res1_0':  11, 'res1_1':  11, 'res1_2':  11, 'res1_3':  11, 'res1_4':  11, 'res1_5':  11, 'res1_6':  11, 'res1_7':  11, 'res1_8':  11, 'res1_9':  11}

lut_bits = {'prop1': 6, 'prop2': 6, 'prop3': 6, 'res0': 6, 'res1': 6}

trash_eta_zones_by_layer = {0: [[208, 240], [352, 512]],
                            1: [[192, 256]],
                            2: [[176, 192]],
                            3: [[32, 48], [176, 192]],
                            4: [[160, 176]],
                            5: [[16, 32], [48, 64], [80, 96], [112, 128], [160, 176]],
                            6: [[16, 32], [144, 160]],
                            7: [[16, 32], [48, 64], [96, 112], [128, 144]],
                            8: [[0, 144]],
                            9: [[0, 144], [272, 320]]}