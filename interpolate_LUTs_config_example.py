# helper config file to be used when running interpolate_LUTs.py
# 
# nbins: number of eta bins to interpolate over
# num_cols: purely cosmetic; number of columns output LUTs will have before going to new line
# factor_bits: dictionary with coeffs as keys telling how many bits to shift the raw parameters when storing them in LUTs
# lut_bits: dictionary with coeffs as keys telling how many bits of precision are to be saved in LUTs
# trash_eta_zones_by_layer: 3D list indexed by layer, then each resulting 2D list contains a sequence of pairs of 
#                            eta values for which the fits should be thrown out and interpolated to instead

nbins = 512
num_cols = 32
factor_bits = {'prop1': 9, 'prop2': 23, 'prop3': 22, 'res0': 0, 'res1': 11}
# factor_bits = {'prop1': 9, 'prop2': 24, 'prop3': 19, 'res0': 0, 'res1': 11}
lut_bits = {'prop1': 6, 'prop2': 6, 'prop3': 6, 'res0': 6, 'res1': 6}

trash_eta_zones_by_layer = [[[208, 256], [352, 512]],
                            [[208, 240]],
                            [[176, 192]],
                            [[176, 192]],
                            [[160, 176]],
                            [[32, 48], [160, 176]],
                            [[144, 160]],
                            [[16, 32], [80, 112]],
                            [[0, 144]],
                            [[0, 144], [272, 320]]]