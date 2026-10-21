# helper config file to be used when running interpolate_LUTs.py
# 
# nbins: number of eta bins to interpolate over
# num_cols: purely cosmetic; number of columns output LUTs will have before going to new line
# factor_bits: dictionary with coeffs as keys telling how many bits to shift the raw parameters when storing them in LUTs
# lut_bits: dictionary with coeffs as keys telling how many bits of precision are to be saved in LUTs
# trash_eta_zones_by_layer_coord1: 3D list indexed by layer, then each resulting 2D list contains a sequence of pairs of 
#                                  eta values for which the fits should be thrown out and interpolated to instead
# trash_eta_zones_by_layer_coord2: 3D list indexed by layer, then each resulting 2D list contains a sequence of pairs of 
#                                  eta values for which the fits should be thrown out and interpolated to instead

nbins = 512
num_cols = 32
factor_bits = {'prop1': 9, 'prop2': 23, 'prop3': 22, 'res0': 0, 'res1': 11}
# factor_bits = {'prop1': 9, 'prop2': 24, 'prop3': 19, 'res0': 0, 'res1': 11}
lut_bits = {'prop1': 6, 'prop2': 6, 'prop3': 6, 'res0': 6, 'res1': 6}

trash_eta_zones_by_layer_coord1 = [[[208, 256]],
                                   [],
                                   [],
                                   [],
                                   [[112, 144], [272, 512]]]

trash_eta_zones_by_layer_coord2 = [[[192, 256], [272, 288]],
                                   [[48, 64], [368, 400]],
                                   [[32, 48], [144, 160]],
                                   [[160, 176], [368, 400]],
                                   [[112, 144], [288, 512]]]

# trash_eta_zones_by_layer_coord1 = [[[208, 256]],
#                                    [],
#                                    [],
#                                    [],
#                                    [[272, 512]]]

# trash_eta_zones_by_layer_coord2 = [[[272, 512]],
#                                    [[368, 400]],
#                                    [],
#                                    [[144, 192], [384, 400]],
#                                    [[112, 144]]]