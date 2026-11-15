import ROOT,itertools,math      
from array import array          
from DataFormats.FWLite import Events, Handle
ROOT.FWLiteEnabler.enable()
import math
verbose=False
import numpy as np
import os
import sys
import argparse
import shlex
from common import (etaLUT, etaShifts, ptShifts, ptLUT, etaLookup, ptLookup, pt_to_gev_factor, fitfuncs_info, 
                    tf_layers, phi_layers, phib_layers, phib_to_phi_layer_map, barrel_limits, load_module_from_file, 
                    prop_coord_axis_titles, get_command_and_commit)



def strAppend(s):
     return "root://cmsxrootd.fnal.gov/"+s

#def strAppend(s):
#     return "root://cmseos.fnal.gov/"+s

def twos_comp(val, bits):
    """compute the 2's complement of int value val"""
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val          

def deltaPhiINT(p1,p2,bits=13):
  dp=p1-p2
  if dp>(1<<(bits-1)):
    return dp - (1<<bits)
  if dp<-(1<<(bits-1)):
    return dp+(1<<bits)
  return dp

def deltaPhi( p1, p2):
     '''Computes delta phi, handling periodic limit conditions.'''
     res = p1 - p2
     while res > math.pi:
         res -= 2*math.pi
     while res < -math.pi:
         res += 2*math.pi
     return res

def deltaR( *args ):
     return math.sqrt( deltaR2(*args) )

def deltaR2( e1, p1, e2, p2):
     de = e1 - e2
     dp = deltaPhi(p1, p2)
     return de*de + dp*dp


def fetchTracks(event):
     trackHandle = Handle('vector<TTTrack<edm::Ref<edm::DetSetVector<Phase2TrackerDigi>,Phase2TrackerDigi,edm::refhelper::FindForDetSetVector<Phase2TrackerDigi> > > >')
     event.getByLabel("l1tTTTracksFromTrackletEmulation:Level1TTTracks",trackHandle)
     return trackHandle.product()
     

def fetchSA(event,tag):
     phiSeg2    = Handle  ('vector<l1t::SAMuon>')
     event.getByLabel(tag,phiSeg2)
     return phiSeg2.product()

def fetchStubs(event,tag):
     phiSeg2    = Handle  ('vector<l1t::MuonStub>')
     event.getByLabel(tag,phiSeg2)
     return phiSeg2.product()


def fetchGEN(event,etaMax=3.0):
     genH  = Handle  ('vector<reco::GenParticle>')
     event.getByLabel('genParticles',genH)
     genMuons=list(filter(lambda x: abs(x.pdgId())==13 and x.status()==1 and abs(x.eta())<etaMax,genH.product()))
     return genMuons

def fetchMatchedTracks(event):
     tracks = fetchTracks(event)
     gen=fetchGEN(event,3.5)
     out=[]
     for t in tracks:
          track_curvature = twos_comp(t.getRinvBits(),15)
          trackpt = pt_to_gev_factor*ptLookup(abs(track_curvature))
          matched=False
          for g in gen:
               if (deltaR(t.eta(),t.phi(),g.eta(),g.phi())<0.05) and (abs(trackpt - g.pt()) < 0.5):
                    matched=True
                    break
          if matched:
               out.append(t)
     return out

def fetchTPS(event,tag,etamax=3.0):
     phiSeg2    = Handle  ('vector<l1t::TrackerMuon>')
     event.getByLabel(tag,phiSeg2)
     return list(filter(lambda x: abs(x.eta())<etamax,phiSeg2.product()))
 


class FuncLoader:

    def __init__(self, func_file, prop_coord, mode):

        self.prop_coord = prop_coord
        assert mode in ['prop', 'res']
        self.mode = mode
        
        if func_file is None:
            self.trivial_func = True   
        else:
            self.trivial_func = False
            func_module = load_module_from_file(func_file)
            self.factor_bits = func_module.factor_bits
            barrel_commands = func_module.barrel_commands
            endcap_commands = func_module.endcap_commands
            
            self.grand_func_LUT = {}
            self.barrel_fitfuncs = {}
            self.barrel_fitfunc_args = {}
            self.barrel_fitfunc_params = {}
            self.endcap_fitfuncs = {}
            self.endcap_fitfunc_args = {}
            self.endcap_fitfunc_params = {}

            from calibration_fit import parseArgs as command_reader

            for layer, command in barrel_commands.items():
                cal_args = command_reader(shlex.split(command)[2:])
                self.barrel_fitfuncs[layer] = fitfuncs_info[cal_args.fitfunc]['fitfunc']
                self.barrel_fitfunc_args[layer] = cal_args.fitfuncArgs
                nparams = len(fitfuncs_info[cal_args.fitfunc]['initparams'])
                if self.mode == 'prop':
                    self.barrel_fitfunc_params[layer] = [f'prop{i+1}_{layer}' for i in range(nparams)]
                elif self.mode == 'res':
                    self.barrel_fitfunc_params[layer] = [f'res{i}_{layer}' for i in range(nparams)]
                for param in self.barrel_fitfunc_params[layer]:
                    self.grand_func_LUT[param] = getattr(func_module, f'lt_{param}')

            for layer, command in endcap_commands.items():
                cal_args = command_reader(shlex.split(command)[2:])
                self.endcap_fitfuncs[layer] = fitfuncs_info[cal_args.fitfunc]['fitfunc']
                self.endcap_fitfunc_args[layer] = cal_args.fitfuncArgs
                nparams = len(fitfuncs_info[cal_args.fitfunc]['initparams'])
                if self.mode == 'prop':
                    self.endcap_fitfunc_params[layer] = [f'prop{i+1}_{layer}' for i in range(nparams)]
                elif self.mode == 'res':
                    self.endcap_fitfunc_params[layer] = [f'res{i}_{layer}' for i in range(nparams)]
                for param in self.endcap_fitfunc_params[layer]:
                    self.grand_func_LUT[param] = getattr(func_module, f'lt_{param}') # repeating some here is fine since LUTs span all eta

            assert 512 % len(getattr(func_module, f'lt_{next(iter(self.grand_func_LUT))}')) == 0
            self.eta_divider = 512 // len(getattr(func_module, f'lt_{next(iter(self.grand_func_LUT))}'))


    def func_eval(self, prop_eta, k, layer, isBarrel):

        if self.trivial_func:
            if self.mode == 'prop':
                return 0
            elif self.mode == 'res':
                return 10000
            else:
                raise Exception 

        params = self.barrel_fitfunc_params[layer] if isBarrel else self.endcap_fitfunc_params[layer]
        p = []
        for param in params:
            p.append(self.grand_func_LUT[param][prop_eta // self.eta_divider] / (2**self.factor_bits[param]))

        if isBarrel:
            out = self.barrel_fitfuncs[layer]([k], p, *self.barrel_fitfunc_args[layer], quantize=True)
        else:
            out = self.endcap_fitfuncs[layer]([k], p, *self.endcap_fitfunc_args[layer], quantize=True)

        return out


    def propagate(self, track_coord_int, prop_eta, k, layer, isBarrel):

        assert self.mode == 'prop', 'propagate only makes sense if the FuncLoader is in prop mode'

        deltacoord = self.func_eval(prop_eta, k, layer, isBarrel)
        
        if self.prop_coord == 'phi':
            if isBarrel and (layer in phib_layers):
                return deltaPhiINT(deltacoord, 0)
            else:
                return deltaPhiINT(track_coord_int + deltacoord, 0)
        else:
            if track_coord_int < 0:
                deltacoord = -deltacoord # etalike coords propagate in their magnitudes
            return track_coord_int + deltacoord



def parseArgs(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['prop', 'res'], help='Whether the input generated is for fitting propagation or resolution')
    parser.add_argument('--propCoord', choices=list(prop_coord_axis_titles.keys()), help='Which variables we are deriving the propagation/resolution for')
    parser.add_argument('--inFile', help='A root file or a python file with an array of root file paths, each containing tracks and stubs to build the inputs from')
    parser.add_argument('--noEosPrefix', action='store_true', help='whether to use raw input file directories (default is to give them a prefix that is needed if they are in fnal eos)')
    parser.add_argument('--inSample', help='which element of inFile to use (if inFile is a python dictionary)')
    parser.add_argument('--phiPropFile', default=None, help='python file (output from interpolate_LUTs.py) containing prop coefficients and functions to update the matching condition for phi')
    parser.add_argument('--phiResFile', default=None, help='python file (output from interpolate_LUTs.py) containing res coefficients and functions to update the matching condition for phi')
    parser.add_argument('--phiResMultiplier', type=float, default=1, help='scale the max allowed phi residual for a match by this value (default: 1)')
    parser.add_argument('--etalikePropFile', default=None, help='python file (output from interpolate_LUTs.py) containing prop coefficients and functions to use for eta-like in res mode (matching condiiton is still wide eta for now)')
    parser.add_argument('--useEtaDiffMatching', action='store_true', help='Use an eta difference as the etalike matching condition (this is default if etalikePropFile is not specified, otherwise etalike matching is done with that propagation)')
    parser.add_argument('--matchEtalikeWindow', type=int, default=261, help='Window of delta etalike that a track-stub pair must fall in to be considered a match (necessary but not sufficeint) (default: 261)')
    parser.add_argument('--usePreviouslyMatchedMuons', action='store_true', help='Only consider stubs from and require tracks to be matched to preexisting tps muons (ignores matching conditions)')
    parser.add_argument('--outdir', default = './', help='output directory')
    parser.add_argument('--postfix', default= '', help='string to append to output root file name')
    parser.add_argument('--maxCounter', type=int, default=240000, help='maximum number of events to consider')
    parser.add_argument('--propEtaBins', type=int, default=512, help='number of bins for prop_eta (default: 512)')
    parser.add_argument('--kBins', type=int, default=81, help='number of bins for k (default: 81) (note: use a power of 3 for rebinning to still be centered at 0)')
    parser.add_argument('--propCoordBins', type=int, default=128, help='number of bins for prop coord (default: 128)')
    parser.add_argument('--propZ0Adjustment', action='store_true', help='whether to adjust for track z0 for eta-like propagations')
    parser.add_argument('--requiredEtaQualities', type=int, nargs='+', help='only allow stubs with these eta qualities')
    args = parser.parse_args(argv)
    return args



def main():

    args = parseArgs()
    mode = args.mode
    prop_coord = args.propCoord
    in_file = args.inFile
    eos_prefix = not args.noEosPrefix
    in_sample = args.inSample
    phi_prop_file = args.phiPropFile
    phi_res_file = args.phiResFile
    phi_res_multiplier = args.phiResMultiplier
    etalike_prop_file = args.etalikePropFile
    use_eta_diff_matching = args.useEtaDiffMatching
    match_etalike_window = args.matchEtalikeWindow
    use_previously_matched_muons = args.usePreviouslyMatchedMuons
    outdir = args.outdir
    postfix = args.postfix
    maxCounter = args.maxCounter
    prop_z0_adjustment = args.propZ0Adjustment
    required_eta_qualities = args.requiredEtaQualities
    
    # read input events
    try:
        in_module = load_module_from_file(in_file)
        to_process = getattr(in_module, in_sample)
        files = []
        for p in to_process:
            if eos_prefix:
                files.append(strAppend(p))
            else:
                files.append(p)
        events = Events(files)
    except:
        events = Events([strAppend(in_file)]) if eos_prefix else Events([f'file:{in_file}'])

    # read prop information if given
    phi_propagator = FuncLoader(phi_prop_file, 'phi', 'prop')
    phi_res_calculator = FuncLoader(phi_res_file, 'phi', 'res')
    etalike_propagator = FuncLoader(etalike_prop_file, 'etalike', 'prop')

    # build etalike propagation and residual helper function
    def etalike_helper(gmt_eta, track_tanl, stub, track_z0, track_z0_cm):
        if prop_coord in ['phi', 'eta']:
            track_etalike = gmt_eta
            stub_etalike = int((1<<12)*stub.offline_eta1()/math.pi)
            if prop_z0_adjustment:
                z0_adj = int((1<<12) * track_z0_cm / np.sqrt(stub.offline_rho()**2 + stub.offline_z()**2) / math.pi)
            else:
                z0_adj = 0
        elif prop_coord == 'sinheta':
            track_etalike = int(track_tanl / 8) #13 bit signed and rounded towards 0 (so this is 2**9 times the phys value, highest allowable phys value is 8)
            stub_etalike = int((1<<9) * stub.offline_z() / stub.offline_rho())
            if prop_z0_adjustment:
                z0_adj = int((1<<9) * track_z0_cm / stub.offline_rho())
            else:
                z0_adj = 0
        elif prop_coord == 'z':
            # match the lsb of track z0 for now
            track_etalike = int((1<11) * stub.offline_rho() * track_tanl_phys / 20.46912512)
            stub_etalike = int((1<11) * stub.offline_z() / 20.46912512)
            if prop_z0_adjustment:
                z0_adj = track_z0
            else:
                z0_adj = 0
        return track_etalike, stub_etalike, z0_adj

    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)

    if prop_coord == 'phi':
        bound = 2048 if mode =='prop' else 512
    elif prop_coord == 'eta':
        bound = 2048 if mode =='prop' else 512
    elif prop_coord == 'sinheta':
        bound = 2048 if mode =='prop' else 512
    elif prop_coord == 'z':
        bound = 128

    calibrationH = {}
    for layer in tf_layers:
        calibrationH[layer] = ROOT.TH3I(f"calibrationH_{layer}", "", args.propEtaBins, 0, 512, args.kBins, -16384, 16384, args.propCoordBins, -bound, bound)
        calibrationH[layer].GetXaxis().SetTitle('prop_eta')
        calibrationH[layer].GetYaxis().SetTitle('track curvature #it{k}')
        calibrationH[layer].GetZaxis().SetTitle(prop_coord_axis_titles[prop_coord]) # this will be changed to bending angle on the fly during fits depending on coord, layer, and eta region
        if mode == 'res':
            calibrationH[layer].GetZaxis().SetTitle(f'res of {h.GetZaxis().GetTitle()}')

    print(f'Running in {mode} mode, outputting to {outdir}')

    totalNumEvents = events.size()
    print(f'Total events in sample: {totalNumEvents}')
    if totalNumEvents > maxCounter:
        print(f'Truncating at {maxCounter} events')

    counter=0
    for event in events:
        if verbose:
            print ("EVENT:",counter)
        counter=counter+1
        tracks=fetchMatchedTracks(event)
        if (counter-1)%10000==0:
            print(f'    Processing events {counter} - {min(counter-1+10000, maxCounter, totalNumEvents)}...')
        if counter==maxCounter:
            break

        for track in tracks:

            #get phi as a signed integer in global coords
            trackerPhi = twos_comp(track.getPhiBits(), 12)
            gmt_phi = deltaPhiINT(track.phiSector()*910 + int(trackerPhi/2), 0)
            
            #calculate signed track curvature integer
            track_curvature = twos_comp(track.getRinvBits(),15)
            
            # now the eta that will be read by LUTs
            track_tanl = (twos_comp(track.getTanlBits(), 16))
            track_tanl_phys = track_tanl / (1<<12)
            abs_gmt_eta = etaLookup(abs(track_tanl))
            gmt_eta = abs_gmt_eta if track_tanl > 0 else -abs_gmt_eta
            prop_eta = abs_gmt_eta // 8
            # propagation isBarrel depends on the track only
            track_is_barrel = {layer: (prop_eta < barrel_limits[layer]) for layer in tf_layers}
            # also may need z0 (maybe need to shift by 0.5 for bin centers when converting to cm, but for now leave it out to emulate hardware)
            track_z0 = twos_comp(track.getZ0Bits(), 12)
            track_z0_cm = track_z0 * (20.46912512 / (1<<11))
            
            # now handle track-stub matching
            # technique is to first define a window of delta etalike that the pair must satisfy, then within all candidates take the one with the closest phi match
            # "delta etalike" here is most intuitively a difference in eta, but if we have a prop etalike file then it is by default the residual in that
            # if a phi res module is given, the candidates are also required to satisfy a delta phi requirement
            # for now, we do not think about using an eta res module here and always use only phi to differentiate candidates
            best_stubs_by_layer = {layer: None for layer in tf_layers}
            phi_diffs_by_layer = {layer: 10000 for layer in tf_layers}
            etalike_diffs_by_layer = {layer: 10000 for layer in tf_layers}
            phi_resids_by_layer = {layer: phi_res_calculator.func_eval(prop_eta, track_curvature, layer, track_is_barrel[layer])*phi_res_multiplier for layer in tf_layers}
            etalike_resids_by_layer = {layer: 10000 for layer in tf_layers}

            if use_previously_matched_muons:

                matchedmu = None
                min_dr = 0.2
                trackpt = pt_to_gev_factor*ptLookup(abs(track_curvature))

                tpsmu=fetchTPS(event,'l1tTkMuonsGmt', 3.5)
                for mu in tpsmu:
                    if abs(trackpt - mu.pt()) > 0.5:
                        continue
                    dr = deltaR(mu.eta(), mu.phi(), track.eta(), track.phi())
                    if dr < min_dr:
                        min_dr = dr
                        matchedmu = mu

                if matchedmu is None:
                    continue

                for stub in matchedmu.stubs():
                    if (stub.quality() & 0x1) != 0:

                        best_stubs_by_layer[stub.tfLayer()] = stub
                        stub_phi = int((1<<12)*stub.offline_coord1()/math.pi)
                        phi_diffs_by_layer[stub.tfLayer()] = deltaPhiINT(stub_phi, gmt_phi)
                        propagated_phi = phi_propagator.propagate(gmt_phi, prop_eta, track_curvature, stub.tfLayer(), track_is_barrel[stub.tfLayer()])
                        phi_resids_by_layer[stub.tfLayer()] = deltaPhiINT(stub_phi, propagated_phi)

                        track_etalike, stub_etalike, z0_adj = etalike_helper(gmt_eta, track_tanl, stub, track_z0, track_z0_cm)
                        etalike_diffs_by_layer[stub.tfLayer()] = (stub_etalike - track_etalike - z0_adj) if gmt_eta > 0 else (-stub_etalike + track_etalike + z0_adj)
                        propagated_etalike = etalike_propagator.propagate(track_etalike, prop_eta, track_curvature, stub.tfLayer(), track_is_barrel[stub.tfLayer()])
                        etalike_resids_by_layer[stub.tfLayer()] = (stub_etalike - propagated_etalike - z0_adj) if gmt_eta > 0 else (-stub_etalike + propagated_etalike + z0_adj)            

            else:

                stubs = fetchStubs(event,'l1tStubsGmt:tps')

                # need to sort them by layer to handle bending angle layers being on their own now
                stubs_by_layer = {i: [] for i in tf_layers}
                for stub in stubs:
                    stubs_by_layer[stub.tfLayer()].append(stub)
                # with current scheme standard ordering would work but this is more general
                ordered_layers = phi_layers + phib_layers

                for layer in ordered_layers:
                    for stub in stubs_by_layer[layer]:
                        if (stub.quality() & 0x1) != 0:

                            stub_phi = int((1<<12)*stub.offline_coord1()/math.pi)
                            propagated_phi = phi_propagator.propagate(gmt_phi, prop_eta, track_curvature, layer, track_is_barrel[layer])
                            phi_resid = deltaPhiINT(stub_phi, propagated_phi)

                            track_etalike, stub_etalike, z0_adj = etalike_helper(gmt_eta, track_tanl, stub, track_z0, track_z0_cm)
                            propagated_etalike = etalike_propagator.propagate(track_etalike, prop_eta, track_curvature, stub.tfLayer(), track_is_barrel[stub.tfLayer()])
                            etalike_resid = (stub_etalike - propagated_etalike - z0_adj) if gmt_eta > 0 else (-stub_etalike + propagated_etalike + z0_adj)
                            if use_eta_diff_matching:
                                match_etalike_resid = int((1<<12)*(stub.offline_eta1() - track.eta())/math.pi)
                            else:
                                match_etalike_resid = etalike_resid
                            
                            if not ((layer in phib_layers) and stub.isBarrel()):
                                if (abs(phi_resid) < abs(phi_resids_by_layer[layer])) and (abs(match_etalike_resid) < match_etalike_window):
                                    best_stubs_by_layer[layer] = stub
                                    phi_diffs_by_layer[layer] = deltaPhiINT(stub_phi, gmt_phi)
                                    phi_resids_by_layer[layer] = phi_resid
                                    etalike_diffs_by_layer[stub.tfLayer()] = (stub_etalike - track_etalike - z0_adj) if gmt_eta > 0 else (-stub_etalike + track_etalike + z0_adj)
                                    etalike_resids_by_layer[stub.tfLayer()] = etalike_resid
                            else:
                                # if this stub has bending angle, check if the corresponding stub carrying phi was matched
                                matched_phi_stub = best_stubs_by_layer[phib_to_phi_layer_map[layer]]
                                if matched_phi_stub is not None:
                                    if stub.address() == matched_phi_stub.address():
                                        best_stubs_by_layer[layer] = stub
                                        phi_diffs_by_layer[layer] = stub_phi
                                        phi_resids_by_layer[layer] = phi_resid
                                        etalike_diffs_by_layer[stub.tfLayer()] = (stub_etalike - track_etalike - z0_adj) if gmt_eta > 0 else (-stub_etalike + track_etalike + z0_adj)
                                        etalike_resids_by_layer[stub.tfLayer()] = etalike_resid
                                        break

            # fill output histograms with results of matching
            for layer, stub in best_stubs_by_layer.items():
                if stub is None:
                    continue
                if stub.etaQuality() not in required_eta_qualities:
                    continue
                if mode == 'prop':
                    if prop_coord == 'phi':
                        deltacoord = phi_diffs_by_layer[layer]
                    else:
                        deltacoord = etalike_diffs_by_layer[layer]
                else:
                    if prop_coord == 'phi':
                        deltacoord = phi_resids_by_layer[layer]
                    else:
                        deltacoord = etalike_resids_by_layer[layer]
                calibrationH[layer].Fill(prop_eta, track_curvature, deltacoord)


    fname = f"calibrationResults_{mode}_{prop_coord}.root" if postfix == '' else f"calibrationResults_{mode}_{prop_coord}_{postfix}.root"
    f=ROOT.TFile(fname,"RECREATE")
    f.cd()
    command, commit = get_command_and_commit()
    ROOT.TNamed('command', command).Write()
    ROOT.TNamed('commit', commit).Write()
    for layer in tf_layers:
        calibrationH[layer].Write()
    f.Close()

if __name__ == "__main__":
    main()