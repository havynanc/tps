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
from common import etaLUT, etaShifts, ptShifts, ptLUT, etaLookup, ptLookup, pt_to_gev_factor, fitfuncs_info, tf_layers, phi_layers, phib_layers, phib_to_phi_layer_map, barrel_limits, load_module_from_file



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
 


class Propagator:

    def __init__(self, prop_file, prop_coord):

        self.prop_coord = prop_coord
        
        if prop_file is None:
            self.trivial_propagate = True   
        else:
            self.trivial_propagate = False
            prop_module = load_module_from_file(prop_file)
            self.factor_bits = prop_module.factor_bits
            barrel_commands = prop_module.barrel_commands
            endcap_commands = prop_module.endcap_commands

            self.coeffs = sorted(list(set([key.split('_')[0] for key in barrel_commands.keys() if not key.startswith(('command', 'commit'))])))
            
            self.grand_prop_LUT = {}
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
                self.barrel_fitfunc_params[layer] = [f'prop{i+1}_{layer}' for i in range(nparams)]
                for param in self.barrel_fitfunc_params[layer]:
                    self.grand_prop_LUT[param] = getattr(prop_module, f'lt_{param}')

            for layer, command in endcap_commands.items():
                cal_args = command_reader(shlex.split(command)[2:])
                self.endcap_fitfuncs[layer] = fitfuncs_info[cal_args.fitfunc]['fitfunc']
                self.endcap_fitfunc_args[layer] = cal_args.fitfuncArgs
                nparams = len(fitfuncs_info[cal_args.fitfunc]['initparams'])
                self.endcap_fitfunc_params[layer] = [f'prop{i+1}_{layer}' for i in range(nparams)]
                for param in self.endcap_fitfunc_params[layer]:
                    self.grand_prop_LUT[param] = getattr(prop_module, f'lt_{param}') # repeating some here is fine since LUTs span all eta

            assert 512 % len(getattr(prop_module, f'lt_{next(iter(self.grand_prop_LUT))}')) == 0
            self.eta_divider = 512 // len(getattr(prop_module, f'lt_{next(iter(self.grand_prop_LUT))}'))


    def propagate(self, track_coord_int, prop_eta, k, layer, isBarrel):

        if self.trivial_propagate:
            return track_coord_int

        params = self.barrel_fitfunc_params if isBarrel else self.endcap_fitfunc_params
        p = []
        for param in params:
            p.append(self.grand_prop_LUT[param][prop_eta // self.eta_divider] / (2**self.factor_bits[param]))

        if isBarrel:
            deltacoord = self.barrel_fitfuncs[layer]([k], p, *self.barrel_fitfunc_args[layer], quantize=True)
        else:
            deltacoord = self.endcap_fitfuncs[layer]([k], p, *self.endcap_fitfunc_args[layer], quantize=True)
        
        if self.prop_coord == 'phi':
            return deltaPhiINT(track_coord_int + deltacoord, 0)
        else:
            if track_coord_int < 0:
                deltacoord = -deltacoord
            return track_coord_int + deltacoord






def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['prop', 'res'], help='Whether the input generated is for fitting propagation or resolution')
    parser.add_argument('--propCoord', choices=['phi', 'eta'], help='Which variables we are deriving the propagation/resolution for')
    parser.add_argument('--inFile', help='A root file or a python file with an array of root file paths, each containing tracks and stubs to build the inputs from')
    parser.add_argument('--noEosPrefix', action='store_true', help='whether to use raw input file directories (default is to give them a prefix that is needed if they are in fnal eos)')
    parser.add_argument('--inSample', help='which element of inFile to use (if inFile is a python dictionary)')
    parser.add_argument('--phiPropFile', default=None, help='python file (output from interpolate_LUTs.py) containing prop coefficients and functions to update the matching condition for phi')
    parser.add_argument('--etaPropFile', default=None, help='python file (output from interpolate_LUTs.py) containing prop coefficients and functions to update the matching condition for eta')
    parser.add_argument('--usePreviouslyMatchedMuons', action='store_true', help='Only consider stubs from and require tracks to be matched to preexisting tps muons')
    parser.add_argument('--outdir', default = './', help='output directory')
    parser.add_argument('--postfix', default= '', help='string to append to output root file name')
    parser.add_argument('--maxCounter', default=240000, help='maximum number of events to consider')
    parser.add_argument('--propEtaBins', type=int, default=512, help='number of bins for prop_eta (default: 512)')
    parser.add_argument('--kBins', type=int, default=81, help='number of bins for k (default: 81) (note: use a power of 3 for rebinning to still be centered at 0)')
    parser.add_argument('--propCoordBins', type=int, default=128, help='number of bins for prop coord (default: 128)')
    args = parser.parse_args()
    return args



def main():

    args = parseArgs()
    mode = args.mode
    prop_coord = args.propCoord
    in_file = args.inFile
    eos_prefix = not args.noEosPrefix
    in_sample = args.inSample
    phi_prop_file = args.phiPropFile
    eta_prop_file = args.etaPropFile
    use_previously_matched_muons = args.usePreviouslyMatchedMuons
    outdir = args.outdir
    postfix = args.postfix
    maxCounter = int(args.maxCounter)
    
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
    phi_propagator = Propagator(phi_prop_file, 'phi')
    eta_propagator = Propagator(eta_prop_file, 'eta')

    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)

    if prop_coord == 'phi':
        bound = 2048 if mode =='prop' else 512
    elif prop_coord == 'eta':
        bound = 64
    calibrationH = {layer: ROOT.TH3I(f"calibrationH_{layer}", "", args.propEtaBins, 0, 512, args.kBins, -16384, 16384, args.propCoordBins, -bound, bound) for layer in tf_layers}

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
            abs_gmt_eta = etaLookup(abs(track_tanl))
            gmt_eta = abs_gmt_eta if track_tanl > 0 else -abs_gmt_eta
            prop_eta = abs_gmt_eta // 8
            # propagation isBarrel depends on the track only
            track_is_barrel = {layer: (prop_eta < barrel_limits[layer]) for layer in tf_layers}
            
            # now handle track-stub matching
            best_stubs_by_layer = {layer: None for layer in tf_layers}
            phi_diffs_by_layer = {layer: 10000 for layer in tf_layers}
            eta_diffs_by_layer = {layer: 10000 for layer in tf_layers}
            phi_resids_by_layer = {layer: 10000 for layer in tf_layers}
            eta_resids_by_layer = {layer: 10000 for layer in tf_layers}

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
                        stub_eta = int((1<<12)*stub.offline_eta1()/math.pi)
                        eta_diffs_by_layer[stub.tfLayer()] = (stub_eta - gmt_eta) if gmt_eta > 0 else (-stub_eta + gmt_eta)
                        propagted_eta = eta_propagator.propagate(gmt_eta, prop_eta, track_curvature, stub.tfLayer(), track_is_barrel[stub.tfLayer()])
                        eta_resids_by_layer[stub.tfLayer()] = (stub_eta - propagated_eta) if gmt_eta > 0 else (-stub_eta + propagated_eta)

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
                            stub_eta = int((1<<12)*stub.offline_eta1()/math.pi)
                            propagated_eta = eta_propagator.propagate(gmt_eta, prop_eta, track_curvature, stub.tfLayer(), track_is_barrel[stub.tfLayer()])
                            eta_resid = (stub_eta - propagated_eta) if gmt_eta > 0 else (-stub_eta + propagated_eta)
                            if not ((layer in phib_layers) and stub.isBarrel() and (prop_coord == 'phi')):
                                if (abs(phi_resid) < abs(phi_resids_by_layer[layer])) and (abs(eta_resid) < (0.2*(1<<12)/math.pi)):
                                    best_stubs_by_layer[layer] = stub
                                    phi_diffs_by_layer[layer] = deltaPhiINT(stub_phi, gmt_phi)
                                    eta_diffs_by_layer[layer] = (stub_eta - gmt_eta) if gmt_eta > 0 else (-stub_eta + gmt_eta)
                                    phi_resids_by_layer[layer] = phi_resid
                                    eta_resids_by_layer[layer] = eta_resid
                            else:
                                # if this stub has bending angle, check if the corresponding stub carrying phi was matched
                                matched_phi_stub = best_stubs_by_layer[phib_to_phi_layer_map[layer]]
                                if matched_phi_stub is not None:
                                    if stub.address() == matched_phi_stub.address():
                                        best_stubs_by_layer[layer] = stub
                                        phi_diffs_by_layer[layer] = stub_phi
                                        eta_diffs_by_layer[layer] = (stub_eta - gmt_eta) if gmt_eta > 0 else (-stub_eta + gmt_eta)
                                        phi_resids_by_layer[layer] = phi_resid
                                        eta_resids_by_layer[layer] = eta_resid
                                        break

            # fill output histograms with results of matching
            for layer, stub in best_stubs_by_layer.items():
                if stub is None:
                    continue
                if mode == 'prop':
                    if prop_coord == 'phi':
                        deltacoord = phi_diffs_by_layer[layer]
                    elif prop_coord == 'eta':
                        deltacoord = eta_diffs_by_layer[layer]
                else:
                    if prop_coord == 'phi':
                        deltacoord = phi_resids_by_layer[layer]
                    elif prop_coord == 'eta':
                        deltacoord = eta_resids_by_layer[layer]
                calibrationH[layer].Fill(prop_eta, track_curvature, deltacoord)


    fname = f"calibrationResults_{mode}_{prop_coord}.root" if postfix == '' else f"calibrationResults_{mode}_{prop_coord}_{postfix}.root"
    f=ROOT.TFile(fname,"RECREATE")
    f.cd()
    for layer in tf_layers:
        calibrationH[layer].Write()
    f.Close()

if __name__ == "__main__":
    main()