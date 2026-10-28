import ROOT,itertools,math
import numpy as np
from DataFormats.FWLite import Events, Handle
ROOT.FWLiteEnabler.enable()
import os
import shutil
import argparse
from common import barrel_limits, fitfuncs_info, get_command_and_commit, tf_layers, phib_layers



# let root see the fitfunc we choose
def initialize_fitfunc(fitfunc_name, mode, *args):

    global fitfunc_py
    def fitfunc_py(x, p):
        return fitfuncs_info[fitfunc_name]['fitfunc'](x, p, *args, quantize=False)
    parlimits = fitfuncs_info[fitfunc_name]['parlimits']
    fitfunc = ROOT.TF1(fitfunc_name, fitfunc_py, -16384, 16384, len(parlimits))
    for i, pair in enumerate(parlimits):
        fitfunc.SetParLimits(i, pair[0], pair[1])
    fitfunc.SetParameters(*fitfuncs_info[fitfunc_name]['initparams'])

    # nice for plotting resolution graphs
    if mode == 'res':
        global refl_fitfunc_py
        def refl_fitfunc_py(x, p):
            return -1 * fitfunc_py(x, p)
        refl_fitfunc = ROOT.TF1('refl_fitfunc', refl_fitfunc_py, -16384, 16384, len(parlimits))

    if mode == 'prop':
        return fitfunc
    else:
        return fitfunc, refl_fitfunc



# cut the tails of a 1d histogram with a symmetrization option to cut on both ends according to the faster decaying end
def dphi_hist_tail_bounds(dphi_hist, cutoff_factor=None, counter_limit=1, force_center_peak=False, dphi_override_radius=None, symmetrization_thresh=np.inf):

    peak_index = dphi_hist.GetMaximumBin() if not force_center_peak else (dphi_hist.GetNbinsX() // 2)
    cutoff = 0 if cutoff_factor is None else (dphi_hist.GetMaximum() // cutoff_factor)
    upper_index = peak_index
    lower_index = peak_index

    empty_counter = 0
    while empty_counter < counter_limit and upper_index < dphi_hist.GetNbinsX() + 2:
        if dphi_hist.GetBinContent(upper_index+1) <= cutoff:
            empty_counter = empty_counter + 1
        else:
            empty_counter = -1
        upper_index = upper_index + 1
    
    empty_counter = 0
    while empty_counter < counter_limit and lower_index > 0:
        if dphi_hist.GetBinContent(lower_index-1) <= cutoff:
            empty_counter = empty_counter + 1
        else:
            empty_counter = -1
        lower_index = lower_index - 1

    if dphi_override_radius is not None:
        while (dphi_hist.GetBinCenter(upper_index) - dphi_hist.GetBinCenter(peak_index) < dphi_override_radius) and (upper_index < dphi_hist.GetNbinsX() + 2):
            upper_index = upper_index + 1
        while (dphi_hist.GetBinCenter(peak_index) - dphi_hist.GetBinCenter(lower_index) < dphi_override_radius) and (lower_index > 0):
            lower_index = lower_index - 1

    #do better at getting rid of noise attached to central hist by symmetrizing cuts
    n = dphi_hist.Integral(lower_index, upper_index)
    if n > symmetrization_thresh:
        if peak_index - lower_index > upper_index - peak_index:
            lower_index = peak_index - (upper_index - peak_index)
            if upper_index - peak_index == 1:
                lower_index = lower_index - 1
        elif upper_index - peak_index > peak_index - lower_index:
            upper_index = peak_index + (peak_index - lower_index)
            if peak_index - lower_index == 1:
                upper_index = upper_index + 1

    return upper_index, lower_index


# look for double peaks in a 1d histogram and cut them away (can happen with noise in bending angle input data)
def dphi_hist_doublepeak_bounds(dphi_hist, upper_index, lower_index):

    search_index = dphi_hist.GetMaximumBin()
    while dphi_hist.GetBinContent(search_index) > dphi_hist.GetBinContent(search_index + 1) and search_index < upper_index - 1:
        search_index = search_index + 1
    upper_index = search_index + 1

    search_index = dphi_hist.GetMaximumBin()
    while dphi_hist.GetBinContent(search_index) > dphi_hist.GetBinContent(search_index - 1) and search_index > lower_index + 1:
        search_index = search_index - 1
    lower_index = search_index - 1

    return upper_index, lower_index


# directly compute params from a 1d histogram
def direct_params(dphi_hist, n, mode):

    mean = dphi_hist.GetMean()
    sigma = dphi_hist.GetStdDev()

    if mode == 'prop':
        param = mean
        param_error = sigma / math.sqrt(n)
    elif mode == 'res':
        param = math.sqrt(sigma**2 + mean**2) #use std relative to 0 here since that's what matters
        if n < 100:
            sigma_error = sigma * math.sqrt(1 - (2/(n-1))*((math.gamma(n/2)/math.gamma((n-1)/2))**2)) #standard error of sample standard deviation 
        else:
            sigma_error = sigma / math.sqrt(2*(n-1)) #asymptotic expansion of gamma functions, need it bc math has finite range and the approx is VERY close
        param_error = math.sqrt((sigma*sigma_error)**2 + (mean*sigma/math.sqrt(n))**2) / param

    return param, param_error


# move the stat box to the upper left corner (e.g. for k > 0, where the dphi peak sits on the right)
def place_stats_box(hist, c, left):
    c.Update() # stats box is only created on paint
    stats = hist.GetListOfFunctions().FindObject('stats')
    if not stats or not left:
        return
    width = stats.GetX2NDC() - stats.GetX1NDC()
    x1 = c.GetLeftMargin() + 0.02
    stats.SetX1NDC(x1)
    stats.SetX2NDC(x1 + width)
    c.Modified()
    c.Update()


# remove errors that can appear for the average dphi at high |k| by comparing points to lower |k| points
# THIS SHOULD NOT BE USED FOR BENDING ANGLE SINCE WE EXPECT SOME BACKWARDS BENDING
def high_k_error_removal(fitslices_graph, mode, curvature_thresh, removal_thresh, max_iterations):

    correction_needed = True
    num_iterations = 0
    while(correction_needed and num_iterations < max_iterations): #don't let it go forever
        correction_needed = False
        num_iterations = num_iterations + 1
        removelist = []
        for l in np.concatenate((np.flip(np.arange(0, fitslices_graph.GetN()//2)), np.arange(fitslices_graph.GetN()//2, fitslices_graph.GetN()))):
            l = int(l)
            curvature = fitslices_graph.GetPointX(l)
            if abs(curvature) > curvature_thresh : 
                ytest = fitslices_graph.Eval(curvature)
                ycompare1 = fitslices_graph.Eval(curvature/2) * 2 * removal_thresh
                ycompare2 = fitslices_graph.Eval(curvature/4) * 4 * removal_thresh
                prop_condition = mode == 'prop' and ((abs(ytest) < abs(ycompare1)) or (abs(ytest) < abs(ycompare2)) 
                                                        or (np.sign(ytest) != np.sign(ycompare1)) or (np.sign(ytest) != np.sign(ycompare2)))
                res_condition = mode == 'res' and ((abs(ytest) > abs(ycompare1)) and (abs(ytest) > abs(ycompare2))
                                                        or (curvature > 0 and fitslices_graph.GetPointX(l-1) < curvature/2)
                                                        or (curvature < 0 and fitslices_graph.GetPointX(l+1) > curvature/2)) 
                if prop_condition or res_condition:
                    removelist.append(l)
                    correction_needed = True
        for l in np.flip(np.sort(removelist)):
            fitslices_graph.RemovePoint(int(l))

    return fitslices_graph   


# use a first round fit as a baseline to redo fits better
def window_fitslices_update(fitslices_graph, mode, fitfunc, hist_2d, window_curv_thresh, radius_upper, radius_lower, removal_thresh, window_folder, c, required_offset=100):

    hist_k = hist_2d.ProjectionX()
    for i in range(hist_k.GetNbinsX()+2):

        center = hist_k.GetBinCenter(i)
        if abs(center) < window_curv_thresh:
            continue
                                        
        #build histogram, indiscriminately including a radius around the old fit
        dphi_hist = hist_2d.ProjectionY(f'center k {center}', i, i)
        for ii in range(0, dphi_hist.GetNbinsX()+2):
            guess_dphi = fitfunc(center)
            bin_dphi = dphi_hist.GetBinCenter(ii)
            if (abs(bin_dphi) < (abs(guess_dphi) - radius_lower)) or (abs(bin_dphi) > (abs(guess_dphi) + radius_upper)) or (abs(bin_dphi) < required_offset) or (np.sign(bin_dphi) != np.sign(guess_dphi)):
                dphi_hist.SetBinContent(ii, 0)
    
        update_point = False
        for point_update_index in range(fitslices_graph.GetN()):
            if fitslices_graph.GetPointX(point_update_index) == center:
                update_point = True
                break
    
        #update num entries, forbid just one or zero entries
        n = dphi_hist.Integral()
        if n < 2:
            if update_point:
                fitslices_graph.RemovePoint(point_update_index)
            continue
        
        dphi_hist.Draw()
        place_stats_box(dphi_hist, c, left=(center<0))
        c.Draw()
        c.SaveAs(f'{window_folder}/{i}.png')
        
        param, param_error = direct_params(dphi_hist, n, mode)
        bin_error = 0.5 * hist_k.GetBinWidth(i)
        
        if update_point:
            fitslices_graph.SetPointY(point_update_index, param)
            fitslices_graph.SetPointError(point_update_index, bin_error, param_error)
        else:
            if abs(param) > abs(fitslices_graph.Eval(center/4) * 4 * removal_thresh): #don't add a new point if it sucks
                fitslices_graph.AddPoint(center, param)
                fitslices_graph.SetPointError(fitslices_graph.GetN()-1, bin_error, param_error)

        dphi_hist.Delete()

    hist_k.Delete()

    return fitslices_graph



def parseArgs(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['prop', 'res'], help='Whether we are fitting propagation or resolution')
    parser.add_argument('--inputFiles', nargs='+', help='Files containing input histograms')
    parser.add_argument('--outdir', default='dy_jpsi_calibration', help='Directory to write outputs in (default: dy_jpsi_calibration)')
    parser.add_argument('--layers', nargs='+', type=int, default=tf_layers, help='Which layers to calibrate (default: all)')
    parser.add_argument('--etaRegion', choices=['barrel', 'endcap', 'inclusive'], default='inclusive', help='whether to fit in barrel, endcap, or both (inclusive) (default: inclusive)')
    parser.add_argument('--overwrite', action='store_true', help='if output directories are allowed to be overwritten (otherwise it is impossible to delete things by mistake)')
    parser.add_argument('--etaWidth', type=int, default=16, help='Width of eta slices used to conduct calibration fits in (default 16)')
    parser.add_argument('--fitfunc', type=str, default=None, help='which function to fit dphi vs k to')
    parser.add_argument('--fitfuncArgs', nargs='+', default=[], help='additional args to pass to fitfunc (look in common.py to see what is needed for each fitfunc)')
    parser.add_argument('--resScaling', type=float, default=3, help='multiply sigma by this to get res. Does nothing if in prop mode (default: 3)')
    parser.add_argument('--tryRangeFit', action='store_true', help='First attempt to fit a gaussian to each 1D dphi histogram before falling back on the direct mean and std')
    parser.add_argument('--fitRadius', type=float, default=150, help='radius to attempt the initial gaussian fit in if tryRangeFit is on (default: 150)')
    parser.add_argument('--fitParamErrorCutoff', type=float, nargs='+', default=[np.inf, 0.5], help='cutoffs for fit parameters from gaussian fits if tryRangeFit is on (default: inf, 0.5)')
    parser.add_argument('--numEntriesCutoff', type=int, default=2, help='minimum number of entries a 1D dphi histogram must have to be considered (default: 2)')
    parser.add_argument('--peakCutoff', type=int, default=0, help='minimum height of the peak a 1D dphi histogram must have to be considered (default: 0)')
    parser.add_argument('--initialKcutoff', type=int, default=13000, help='highest |k| bin midpoint to attempt to calibrate over in initial pass (before inclusive) (default: 13000)')
    parser.add_argument('--minFitslicesPoints', type=int, default=4, help='minimum number of points (slices of k) we must successfully fit in order to perform the overall 2D eta slice histogram fit (default: 4)')
    parser.add_argument('--errorRemovalMaxIterations', type=int, default=20, help='how many times are we allowed to loop through the fitslices error removal algo (0 to skip it, for example for bending angle) (default: 20)')
    parser.add_argument('--errorRemovalCurvatureThresh', type=int, default=None, help='|k| above this will be subject to error removal (default: 3000 for prop, infinity for res)')
    parser.add_argument('--errorRemovalThresh', type=float, default=None, help='scale the error removal comparison points (adjusted for k ratio) (default: 0.9 for prop, 3 for res)')
    parser.add_argument('--windowCurvThresh', type=int, default=None, help='|k| above this will find prop datapoints more inclusively than usual to accomodate sparse data (default: 5000 for prop, infinity for res)')
    parser.add_argument('--windowRadiusLower', type=int, default=200, help='dphi radius about initial fit to consider more inclusive dphi points, in the direction pointing toward 0 (default: 200)')
    parser.add_argument('--windowRadiusUpper', type=int, default=200, help='dphi radius about initial fit to consider more inclusive dphi points, in the direction pointing away from 0 (default: 200)')
    parser.add_argument('--symmetrizationThresh', type=int, default=10, help='minimum number of data points (after initial noise removal) that a 1D dphi hist can have for it to undergo symmetrization (cutting a skewed tail) (default: 10)')
    parser.add_argument('--resReflectSym', action='store_true', help='get better stats in res mode by combining k and -k bins')
    parser.add_argument('--dropCentralKbin', action='store_true', help='drop central k bin for res')
    parser.add_argument('--doUnweightedFit', action='store_true', help='do unweighted fit (useful for res)')
    args = parser.parse_args(argv)
    return args



def main():

    #stop the root printout from spamming me
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
    #dont show me canvases
    ROOT.gROOT.SetBatch(1)
    ROOT.gStyle.SetStatX(0.83)
    ROOT.gStyle.SetStatY(0.89)
    ROOT.gStyle.SetTitleOffset(1.4, 'Y')

    args = parseArgs()
    mode = args.mode
    input_files = args.inputFiles
    outdir = args.outdir
    layers = args.layers
    eta_region = args.etaRegion
    overwrite = args.overwrite
    etaWidth = args.etaWidth
    res_scaling = args.resScaling
    coarsening_ratio = [1, 9, 2] if mode == 'prop' else [1, 9, 1] #eta, k, dphi

    if mode == 'prop':
        fitfunc = initialize_fitfunc(args.fitfunc, mode, *args.fitfuncArgs)
    else:
        fitfunc, refl_fitfunc = initialize_fitfunc(args.fitfunc, mode, *args.fitfuncArgs)

    # Specify the approach to be taken in order to compute the "average" dphi over each 1D slice of eta and k (dphi histograms)
    # True: first attempt to fit a gaussian around the histogram's mean. If this fails (outright faiure or high error in fit params), fall back on the following
    # False: considering only the bins touching the central peak (mode), directly compute the mean and std. If there is no data, no value is reported.
    #        Also discard skewed tails if there are many data points and discard secondary peaks for phib.

    try_range_fit = args.tryRangeFit
    fit_radius = args.fitRadius
    fit_param_error_cutoff = args.fitParamErrorCutoff

    num_entries_cutoff = args.numEntriesCutoff
    peak_cutoff = args.peakCutoff
    initial_k_cutoff = args.initialKcutoff
    min_fitslices_points = args.minFitslicesPoints
    error_removal_max_iterations = args.errorRemovalMaxIterations
    error_removal_curvature_thresh = args.errorRemovalCurvatureThresh if args.errorRemovalCurvatureThresh is not None else 3000 if (mode == 'prop') else np.inf
    error_removal_thresh = args.errorRemovalThresh if args.errorRemovalThresh is not None else 0.9 if (mode == 'prop') else 3

    window_curv_thresh = args.windowCurvThresh if args.windowCurvThresh is not None else 5000 if (mode == 'prop') else np.inf
    window_radius_lower = args.windowRadiusLower
    window_radius_upper = args.windowRadiusUpper

    symmetrization_thresh = args.symmetrizationThresh
    res_reflect_sym = args.resReflectSym
    drop_central_k_bin = args.dropCentralKbin
    do_unweighted_fit = args.doUnweightedFit


    c = ROOT.TCanvas('canvas', 'canvastitle', 1200, 800)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.15)
    c.SetLogz()

    rfiles = []
    for input_file in input_files:
        rfiles.append(ROOT.TFile.Open(input_file, 'READ'))
    hist_names = [key.GetName() for key in rfiles[0].GetListOfKeys()]

    # make output directories and throw an error if we would accidentally overwrite one of the layer directories
    out_paths = []
    for layer in layers:
        out_paths.append(f'{outdir}/{mode}/{eta_region}/layer{layer}')
    for out_path in out_paths:
        if os.path.isdir(out_path):
            if overwrite:
                shutil.rmtree(out_path)
            else:
                raise Exception(f'output folder {out_path} already exists and --overwrite was not specified')
        os.makedirs(out_path)
    
    os.chdir(outdir)
    os.chdir(mode)
    os.chdir(eta_region)

    # first set up the output graphs, organized by fit param then by layer
    num_params = fitfunc.GetNpar()
    coeff_graphs = [[] for _ in range(num_params)]
    for layer in tf_layers:
        for i in range(num_params):
            index = i+1 if mode == 'prop' else i
            g = ROOT.TGraphErrors()
            g.SetName(f'{mode}{index}_{layer}')
            g.SetTitle(f'{mode}{index} layer {layer} vs eta')
            g.GetXaxis().SetTitle('prop_eta')
            g.GetYaxis().SetTitle(f'{mode}{index}') 
            g.SetMarkerStyle(20)
            g.SetMarkerSize(1)
            g.SetFillStyle(10)
            coeff_graphs[i].append(g)

    for name in hist_names:

        #check if we are considering this layer (thank god the new layer is still one digit)
        layer = int(name[-1])
        if layer not in layers:
            continue

        print(f'Processing {name} (layer {layer}):') 

        #get data
        h = rfiles[0].Get(name)
        for rfile in rfiles[1:]:
            h_temp = rfile.Get(name)
            h.Add(h_temp)
            h_temp.Delete()

        h.GetXaxis().SetTitle('prop_eta')
        h.GetYaxis().SetTitle('track_curvature')
        #label for 'dphi' axis depends on what exactly we plotted
        if layer not in phib_layers :
            h.GetZaxis().SetTitle('deltaphi')
            #h.GetZaxis().SetTitle(f'{coord} (deta)')
        else:
            h.GetZaxis().SetTitle('deltaphi for endcap, phib for barrel')
            #h.GetZaxis().SetTitle(f'{coord} (deta)')
        if mode == 'res':
            h.GetZaxis().SetTitle(h.GetZaxis().GetTitle() + ' res')
            
        h.Rebin3D(coarsening_ratio[0], coarsening_ratio[1], coarsening_ratio[2]) #not sure about the errors computed with this method
        k_axis = h.GetYaxis()
        os.chdir(f'layer{layer}')
        os.mkdir('overlays')       

        #make some diagnostic plots 
        hproj_eta = h.ProjectionX()
        hproj_k = h.ProjectionY()
        hproj_dphi = h.ProjectionZ()
        hproj_eta.SetTitle(hproj_eta.GetTitle() + ' ' + h.GetXaxis().GetTitle())
        hproj_k.SetTitle(hproj_k.GetTitle() + ' ' + h.GetYaxis().GetTitle())
        hproj_dphi.SetTitle(hproj_dphi.GetTitle() + ' ' + h.GetZaxis().GetTitle())
            
        hproj_eta.Draw()
        c.Draw()
        c.SaveAs(f'{hproj_eta.GetTitle()}.png')
        c.Clear()
        hproj_eta.Delete()
        hproj_k.Draw()
        c.Draw()
        c.SaveAs(f'{hproj_k.GetTitle()}.png')
        c.Clear()
        hproj_k.Delete()
        hproj_dphi.Draw()
        c.Draw()
        c.SaveAs(f'{hproj_dphi.GetTitle()}.png')
        c.Clear()
        hproj_dphi.Delete()

        #now do fits in slices of eta
        #normal layers - width 16 or 32, full range 0-512 (511 is highest)
        etaList = np.arange(0, 512+1, etaWidth)
            
        for j in range(len(etaList) - 1):
        
            etaLower = int(etaList[j])
            etaUpper = int(etaList[j+1] - 1) #inclusive on the right bound for some reason

            #see if we're in barrel, endcap, or including sections of both (skip these sections)
            #if (etaLower < barrel_limits[layer]) and (etaUpper > barrel_limits[layer]):
            #    continue
            #determine with midpoint - plays nicely with the specific barrel limits we have
            isBarrel = (etaLower + etaWidth/2) < barrel_limits[layer]
            if (isBarrel and eta_region == 'endcap') or ((not isBarrel) and eta_region == 'barrel'):
                continue
            
            print(f'    Slice {etaLower} <= eta < {etaUpper+1}')
            
            #extract the slice and skip to the next if it's empty
            h.GetXaxis().SetRange(etaLower+1, etaUpper+1)    
            hist_2d = h.Project3D('zy')
            if hist_2d.GetEntries() == 0:
                continue

            #save some diagnostic plots                
            dirname = f'eta{etaLower}-{etaUpper}'
            os.mkdir(dirname)
            os.chdir(dirname)
        
            hist_2d.Draw()
            c.Draw()
            c.SaveAs(f'{hist_2d.GetName()}.png')
        
            hist_2d_dphi = hist_2d.ProjectionY()
            hist_2d_dphi.Draw()
            c.Draw()
            c.SaveAs(f'{hist_2d_dphi.GetName()}.png')
            hist_2d_dphi.Delete()
                
            hist_2d_k = hist_2d.ProjectionX()     
            hist_2d_k.Draw()
            c.Draw()
            c.SaveAs(f'{hist_2d_k.GetName()}.png')
            hist_2d_k.Delete()      
                
            #set up output folders and graphs for the main computations now
            hist_folder = 'dphi_hists'
            raw_folder = hist_folder + '/raw'
            range_fit_folder = hist_folder + '/range_fit'
            direct_params_folder = hist_folder + '/direct_params'
            window_folder = hist_folder + '/window'
            os.mkdir(hist_folder)
            os.mkdir(raw_folder)
            os.mkdir(range_fit_folder)
            os.mkdir(direct_params_folder)
            os.mkdir(window_folder)
        
            ROOT.gStyle.SetOptFit(111) #displays box with fit parameters

            fitslices_graph = ROOT.TGraphErrors()
            fitslices_graph.GetXaxis().SetTitle('track_curvature')
            fitslices_graph.SetTitle(f'Layer {layer}, eta {etaLower}-{etaUpper} fitted {mode}')
            fitslices_graph.GetYaxis().SetTitle(f'fitted {mode}')

            #more customizables version of fitslices now begins

            #go through each k slice and compute the relevant parameter on the dphi histogram
            for i in range(hist_2d.GetNbinsX()+2):
                    
                if drop_central_k_bin and mode == 'res' and i == hist_2d.GetNbinsX()//2 + 1:
                    continue
                dphi_hist = hist_2d.ProjectionY(f'center k {k_axis.GetBinCenter(i)}', i, i)
                if mode == 'res' and res_reflect_sym:
                    dphi_hist.Add(hist_2d.ProjectionY('', hist_2d.GetNbinsX()+1-i, hist_2d.GetNbinsX()+1-i))
                # n = dphi_hist.GetEntries()
                n = dphi_hist.Integral()
                peakheight = dphi_hist.GetMaximum()
                if n == 0:
                    continue
                dphi_hist.Draw()
                place_stats_box(dphi_hist, c, left=(k_axis.GetBinCenter(i)<0))
                c.Draw()
                c.SaveAs(f'{raw_folder}/{i}.png')
                if n < num_entries_cutoff or peakheight < peak_cutoff:
                    continue

                #for phib, don't look at curvatur too high to be seen (for deltaphi, error removal would take care of this)
                if abs(k_axis.GetBinCenter(i)) > initial_k_cutoff:
                    continue 
                
                #The logic here is to try to fit a gaussian in a range around the peak to hopefully capture the bell curve (need the range to cast aside noise).
                #If the fit fails (either directly or by having a large standard error in its parameters) cut the histogram beyond its central bulge,
                #    then exact mean and stddev and their standard errors.
                #If the fit succeeds, report the fitted mean and standard deviation and their errors.
                #In order to even be considered, a slice needs to satisfy the minimum numEntries and peak size requirements.
                range_fit_status = False
                if try_range_fit:                
                    #peak = dphi_hist.GetMean()
                    peak = dphi_hist.GetBinCenter(dphi_hist.GetMaximumBin())
                    sub_gaus = ROOT.TF1('sub_gaus', 'gaus', peak - fit_radius, peak + fit_radius)
                    #sub_gaus.SetParLimits(1, -100, 100) #needed for deta fitting I think
                    result = dphi_hist.Fit('sub_gaus', option='Q R S')
                    mean = result.Parameter(1)
                    mean_error = result.ParError(1)
                    sigma = result.Parameter(2)
                    sigma_error = result.ParError(2)
                    if mode == 'prop':
                        param = mean
                        param_error = mean_error
                    elif mode =='res':
                        param = sigma
                        param_error = sigma_error
                    dphi_hist.Draw()
                    place_stats_box(dphi_hist, c, left=(k_axis.GetBinCenter(i)<0))
                    c.Draw()
                    c.SaveAs(f'{range_fit_folder}/{i}.png')
                    range_fit_status = not ((str(result) == '<nullptr TFitResult>') or (result.ParError(1) >= fit_param_error_cutoff[0] * abs(result.Parameter(1))) or (result.ParError(2) >= fit_param_error_cutoff[1] * abs(result.Parameter(2))))
                    
                if (not range_fit_status):
                    #fit failed/wasn't tried. directly use sample mean and stddv and the standard errors thereof after trimming past the central bulge
                        
                    #look at everything if curvature is high enough (for phib)
                    trim = not (abs(k_axis.GetBinCenter(i)) > 3000)# and isBarrel)# and coord == 'coord2')
                    lower_index = -1
                    upper_index = dphi_hist.GetNbinsX() + 2

                    if trim:

                        if mode == 'prop':
                            cutoff_factor = None
                            counter_limit = 1
                            force_center_peak = False
                            dphi_override_radius = None
                            symmetrization_thresh = symmetrization_thresh
                        elif mode == 'res':
                            cutoff_factor = 20
                            counter_limit = 5 if (abs(k_axis.GetBinCenter(i)) < 3000) else 20
                            force_center_peak = True
                            dphi_override_radius = 200 if (abs(k_axis.GetBinCenter(i)) > 3000 and isBarrel) else None
                            symmetrization_thresh = np.inf

                        upper_index, lower_index = dphi_hist_tail_bounds(dphi_hist, cutoff_factor=cutoff_factor, counter_limit=counter_limit, force_center_peak=force_center_peak, dphi_override_radius=dphi_override_radius, symmetrization_thresh=symmetrization_thresh)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

                        #try to remedy double peaks from noise and signal being close in phib
                        if mode == 'prop' and (layer in phib_layers):
                            upper_index, lower_index = dphi_hist_doublepeak_bounds(dphi_hist, upper_index, lower_index)

                        #now actually trim the data
                        for ii in range(0, lower_index+1): #changed to lower_index+1
                            dphi_hist.SetBinContent(ii, 0)
                        for ii in range(upper_index, dphi_hist.GetNbinsX()+2): #changed to upper_index from upper_index+1
                            dphi_hist.SetBinContent(ii, 0)

                        n = dphi_hist.Integral(lower_index, upper_index)

                    #in res, histograms with just one nonempty bin are problematic because they give 0 error (root treats all entries in a bin as having the same value).
                    #should really fix the error calculation for these cases in prop as well. throw them out for now.
                    #also require a minimum number of entries since noise getting mischaracterized always throws trend slope UP - doesn't average out like in prop
                    if mode == 'res' and ((trim and (lower_index + 2 == upper_index)) or n < 5):
                        continue

                    # forbid just one or zero entries
                    if n < 2:
                        continue
    
                    #compute params
                    param, param_error = direct_params(dphi_hist, n, mode)

                    dphi_hist.Draw()
                    place_stats_box(dphi_hist, c, left=(k_axis.GetBinCenter(i)<0))
                    c.Draw()
                    c.SaveAs(f'{direct_params_folder}/{i}.png')

                dphi_hist.Delete()
                                
                #fill the graph of (mean vs k) or (res vs k) within the eta slice with this newest point
                center = k_axis.GetBinCenter(i)
                bin_error = 0.5 * k_axis.GetBinWidth(i)
                if mode == 'res':
                    param = res_scaling * param
                    param_error = res_scaling * param_error
                fitslices_graph.AddPoint(center, param)
                fitslices_graph.SetPointError(fitslices_graph.GetN()-1, bin_error, param_error)

            #sparsely populated high curvature points can sometimes have bad values, characterized by being way too low in magnitude (peak happened to be noise)
            #here we remove them
            fitslices_graph = high_k_error_removal(fitslices_graph, mode, error_removal_curvature_thresh, error_removal_thresh, error_removal_max_iterations)

            fitslices_graph.Draw('AP')
            options_str = 'Q S EX0' #need to use EX0 for fitting procedure to actually work - ignores x errors
            if mode == 'res' and do_unweighted_fit:
                options_str = options_str + ' W' #large |k| are very important to get right
            fitslices_fit = fitslices_graph.Fit(fitfunc, options_str, '', -16384, 16384)
                
            #update high curvature points using this fit (using more inclusive window criteria), then fit again
            if mode == 'prop':# and isBarrel:
                fitslices_graph = window_fitslices_update(fitslices_graph, mode, fitfunc, hist_2d, window_curv_thresh, window_radius_upper, window_radius_lower, error_removal_thresh, window_folder, c, required_offset=100)
                try:
                    fitslices_graph.GetFunction(fitfunc.GetName()).Delete()
                except:
                    pass
                fitslices_fit = fitslices_graph.Fit(fitfunc, options_str, '', -16384, 16384)
                
            #fitfunc.Draw()
            c.Draw()
            c.SaveAs(f'fitslices_{mode}.png')

            #add these fitted params to the output graphs if they have enough points 
            if (fitslices_graph.GetN() >= min_fitslices_points):
                for k, layer_coeff_graphs in enumerate(coeff_graphs):
                    layer_coeff_graphs[layer].AddPoint(etaLower + etaWidth/2, fitslices_fit.Parameter(k))
                    layer_coeff_graphs[layer].SetPointError(layer_coeff_graphs[layer].GetN()-1, etaWidth/2, fitslices_fit.ParError(k))
        
            #get the 2d hist and fit graph overlaid to check that it makes sense
            hist_2d.SetTitle(f'Layer {layer}, {etaLower} <= eta < {etaUpper+1}')
            hist_2d.SetStats(0)
            hist_2d.Draw('COLZ')
            # fitslices_graph.SetLineColor(ROOT.TColor.GetColorTransparent(9, 1))
            fitslices_graph.Draw('s')
            if mode == 'res':
                try:
                    refl_fitfunc.SetParameters(*fitslices_fit.Parameters())
                    refl_fitfunc.Draw('same')
                except:
                    pass
            c.Draw()    
            c.SaveAs(f'{mode}_hist_overlay.png')
            c.SaveAs(f'../overlays/{j}.png')            

            hist_2d.Delete()
            os.chdir('..') # now finished for this eta slice


        
        #make param vs eta graphs
        for k, layer_coeff_graphs in enumerate(coeff_graphs):
            index = k+1 if (mode == 'prop') else k
            layer_coeff_graphs[layer].Draw('AP')
            c.Draw()
            c.SaveAs(f'{mode}{index}_{layer}.png')
        
        h.Delete()
        os.chdir('..') # now finished for this hist
        

        
    f=ROOT.TFile(f"{mode}_graphs.root","UPDATE")
    f.cd()
    command, commit = get_command_and_commit()
    for layer_coeff_graphs in coeff_graphs:
        for layer in layers:
            g = layer_coeff_graphs[layer]
            f.Delete(f'{g.GetName()};*')
            g.Write()
            f.Delete(f'command_{g.GetName()};*')
            f.Delete(f'commit_{g.GetName()};*')
            ROOT.TNamed(f'command_{g.GetName()}', command).Write()
            ROOT.TNamed(f'commit_{g.GetName()}', commit).Write()
    f.Close()



if __name__ == "__main__":
    main()
