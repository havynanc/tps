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
from common import etaLUT, etaShifts, ptShifts, ptLUT, etaLookup, ptLookup, pt_to_gev_factor, fitfuncs_info, tf_layers, barrel_limits, load_module_from_file



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
     gen=fetchGEN(event,2.4)
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
 

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['prop', 'res'], help='Whether the input generated is for fitting propagation or resolution')
    parser.add_argument('--propcoord', choices=['phi', 'eta'], help='Which variables we are deriving the propagation/resolution for')
    parser.add_argument('--inFile', help='A root file or a python file with an array of root file paths, each containing tracks and stubs to build the inputs from')
    parser.add_argument('--noEosPrefix', action='store_true', help='whether to use raw input file directories (default is to give them a prefix that is needed if they are in fnal eos)')
    parser.add_argument('--inSample', help='which element of inFile to use (if inFile is a python dictionary)')
    parser.add_argument('--propFileCoord1', default=None, help='python file (output from interpolate_LUTs.py) containing prop coefficients and functions to update the matching condition')
    parser.add_argument('--propFileCoord2', default=None, help='python file (output from interpolate_LUTs.py) containing prop coefficients and functions to update the matching condition')
    parser.add_argument('--usePreviouslyMatchedMuons', action='store_true', help='Only consider stubs from and require tracks to be matched to preexisting tps muons')
    parser.add_argument('--outdir', default = './', help='output directory')
    parser.add_argument('--postfix', default= '', help='string to append to output root file name')
    parser.add_argument('--maxCounter', default=240000, help='maximum number of events to consider')
    args = parser.parse_args()
    return args



def main():

    args = parseArgs()
    mode = args.mode
    propcoord = args.propcoord
    in_file = args.inFile
    eos_prefix = not args.noEosPrefix
    in_sample = args.inSample
    prop_file_coord1 = args.propFileCoord1
    prop_file_coord2 = args.propFileCoord2
    assert (((prop_file_coord1 is None) and (prop_file_coord2 is None)) or ((prop_file_coord1 is not None) and (prop_file_coord2 is not None)))
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
    if prop_file_coord1 is not None:

        prop_module_coord1 = load_module_from_file(prop_file_coord1)
        factor_bits_coord1 = prop_module_coord1.factor_bits
        barrel_commands_coord1 = prop_module_coord1.barrel_commands
        endcap_commands_coord1 = prop_module_coord1.endcap_commands

        prop_module_coord2 = load_module_from_file(prop_file_coord2)
        factor_bits_coord2 = prop_module_coord2.factor_bits
        barrel_commands_coord2 = prop_module_coord2.barrel_commands
        endcap_commands_coord2 = prop_module_coord2.endcap_commands

        coeffs = sorted(list(set([key.split('_')[0] for key in barrel_commands_coord1.keys() if not key.startswith(('command', 'commit'))])))
        
        grand_prop_LUT = {}
        barrel_fitfuncs_coord1 = {}
        barrel_fitfunc_args_coord1 = {}
        endcap_fitfuncs_coord1 = {}
        endcap_fitfunc_args_coord1 = {}
        barrel_fitfuncs_coord2 = {}
        barrel_fitfunc_args_coord2 = {}
        endcap_fitfuncs_coord2 = {}
        endcap_fitfunc_args_coord2 = {}
        from calibration_fit import parseArgs as command_reader
        for name, command in barrel_commands_coord1.items():
            # structure is such that barrel and endcap commands lists both have all the same names
            grand_prop_LUT[name] = getattr(prop_module_coord1, f'lt_{name}')
            cal_args = command_reader(shlex.split(command)[2:])
            barrel_fitfuncs_coord1[name] = fitfuncs_info[cal_args.fitfunc]['fitfunc']
            barrel_fitfunc_args_coord1[name] = cal_args.fitfuncArgs
        for name, command in endcap_commands_coord1.items():
            cal_args = command_reader(shlex.split(command)[2:])
            endcap_fitfuncs_coord1[name] = fitfuncs_info[cal_args.fitfunc]['fitfunc']
            endcap_fitfunc_args_coord1[name] = cal_args.fitfuncArgs
        for name, command in barrel_commands_coord2.items():
            # structure is such that barrel and endcap commands lists both have all the same names
            grand_prop_LUT[name] = getattr(prop_module_coord2, f'lt_{name}')
            cal_args = command_reader(shlex.split(command)[2:])
            barrel_fitfuncs_coord2[name] = fitfuncs_info[cal_args.fitfunc]['fitfunc']
            barrel_fitfunc_args_coord2[name] = cal_args.fitfuncArgs
        for name, command in endcap_commands_coord2.items():
            cal_args = command_reader(shlex.split(command)[2:])
            endcap_fitfuncs_coord2[name] = fitfuncs_info[cal_args.fitfunc]['fitfunc']
            endcap_fitfunc_args_coord2[name] = cal_args.fitfuncArgs

        assert 512 % len(getattr(prop_module_coord1, f'lt_{next(iter(grand_prop_LUT))}')) == 0
        eta_divider = 512 // len(getattr(prop_module_coord1, f'lt_{next(iter(grand_prop_LUT))}'))

        def propagate(track_phi_int, prop_eta, k, layer, coord, isBarrel):

            p = []
            for i, coeff in enumerate(coeffs):
                name = f'{coeff}_{coord}_{layer}'
                this_factor_bits = factor_bits_coord1[coeff] if coord == 'coord1' else factor_bits_coord2[coeff]
                p.append(grand_prop_LUT[name][prop_eta // eta_divider] / (2**this_factor_bits))

            # print('NEW PROPAGATION')
            # latest value of name is fine here
            if isBarrel:
                if coord == 'coord1':
                    deltaphi = barrel_fitfuncs_coord1[name]([k], p, *barrel_fitfunc_args_coord1[name], quantize=True)
                    # print('fitfunc args: ', barrel_fitfunc_args_coord1[name])
                else:
                    deltaphi = barrel_fitfuncs_coord2[name]([k], p, *barrel_fitfunc_args_coord2[name], quantize=True)
                    # print('fitfunc args: ', barrel_fitfunc_args_coord2[name])
            else:
                if coord == 'coord1':
                    deltaphi = endcap_fitfuncs_coord1[name]([k], p, *endcap_fitfunc_args_coord1[name], quantize=True)
                    # print('fitfunc args: ', endcap_fitfunc_args_coord1[name])
                else:
                    deltaphi = endcap_fitfuncs_coord2[name]([k], p, *endcap_fitfunc_args_coord2[name], quantize=True)
                    # print('fitfunc args: ', endcap_fitfunc_args_coord2[name])
                
            # print('k: ', k)
            # print('eta: ', prop_eta)
            # print('layer: ', layer)
            # print('coord: ', coord)
            # print('isBarrel: ', isBarrel)
            # print('p: ', p)
            # print('deltaphi: ', deltaphi)

            return deltaPhiINT(track_phi_int + deltaphi, 0)

    else:

        def propagate(track_phi_int, prop_eta, k, layer, coord, isBarrel):
            return track_phi_int

    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)

    #changed lower limit from -16383 to -16384 and -1024 to 1024 bc it makes the bins soooo much nicer
    #now do the kbinning here explicitly - 64+1 bins to center at limit 16384+512/2 to have bin widths of 512 (original k value were 512, -16383, 1684)
    #sike I realized the best way to keep it centered at 0 is to use powers of 3
    #expanded dphi hist in hopes of catching higher values for higher curvatures (1024->2048, keep 256 bc memory issues (could fix by doing 1 layer at a time)

    if propcoord == 'phi':
        bound = 2048 if mode =='prop' else 512
        nbins = 256
    elif propcoord == 'eta':
        bound = 64
        nbins = 128
    calibrationH1 = {layer: ROOT.TH3I(f"calibrationH1_{layer}","",512,0,512,729,-16384,16384,nbins,-bound-0.5,bound-0.5) for layer in tf_layers}
    calibrationH2 = {layer: ROOT.TH3I(f"calibrationH2_{layer}","",512,0,512,729,-16384,16384,nbins,-bound-0.5,bound-0.5) for layer in tf_layers}

    #eta identity prop
    '''
    calibrationH1={
        0:ROOT.TH2I("calibrationH1_0","",1024,0,4096,1024,0,4096),
        1:ROOT.TH2I("calibrationH1_1","",1024,0,4096,1024,0,4096),
        2:ROOT.TH2I("calibrationH1_2","",1024,0,4096,1024,0,4096),
        3:ROOT.TH2I("calibrationH1_3","",1024,0,4096,1024,0,4096),
        4:ROOT.TH2I("calibrationH1_4","",1024,0,4096,1024,0,4096),
    }
    calibrationH2={
        0:ROOT.TH2I("calibrationH2_0","",1024,0,4096,1024,0,4096),
        1:ROOT.TH2I("calibrationH2_1","",1024,0,4096,1024,0,4096),
        2:ROOT.TH2I("calibrationH2_2","",1024,0,4096,1024,0,4096),
        3:ROOT.TH2I("calibrationH2_3","",1024,0,4096,1024,0,4096),
        4:ROOT.TH2I("calibrationH2_4","",1024,0,4096,1024,0,4096),
    }


    #res in units of stub eta - lower range, higher resolution
    calibrationH1={
        0:ROOT.TH3I("calibrationH1_0","",512,0,512,729,-16384,16384,64,-32,32),
        1:ROOT.TH3I("calibrationH1_1","",512,0,512,729,-16384,16384,64,-32,32),
        2:ROOT.TH3I("calibrationH1_2","",512,0,512,729,-16384,16384,64,-32,32),
        3:ROOT.TH3I("calibrationH1_3","",512,0,512,729,-16384,16384,64,-32,32),
        4:ROOT.TH3I("calibrationH1_4","",512,0,512,729,-16384,16384,64,-32,32)
    }
    calibrationH2={
        0:ROOT.TH3I("calibrationH2_0","",512,0,512,729,-16384,16384,64,-32,32),
        1:ROOT.TH3I("calibrationH2_1","",512,0,512,729,-16384,16384,64,-32,32),
        2:ROOT.TH3I("calibrationH2_2","",512,0,512,729,-16384,16384,64,-32,32),
        3:ROOT.TH3I("calibrationH2_3","",512,0,512,729,-16384,16384,64,-32,32),
        4:ROOT.TH3I("calibrationH2_4","",512,0,512,729,-16384,16384,64,-32,32)
    }
    '''

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
            gmtphi = deltaPhiINT(track.phiSector()*910 + int(trackerPhi/2), 0)
            
            #calculate signed track curvature integer
            track_curvature = twos_comp(track.getRinvBits(),15)
            
            # now the eta that will be read by LUTs
            track_tanl = (twos_comp(track.getTanlBits(), 16))
            gmt_eta = etaLookup(abs(track_tanl))
            prop_eta = gmt_eta // 8
            # propagation isBarrel depends on the track only
            track_is_barrel = {layer: (prop_eta < barrel_limits[layer]) for layer in tf_layers}
            
            # now handle track-stub matching
            best_stubs_by_layer_coord1 = {layer: None for layer in tf_layers}
            best_stubs_by_layer_coord2 = {layer: None for layer in tf_layers}
            phi_resids_by_layer_coord1 = {layer: 10000 for layer in tf_layers}
            phi_resids_by_layer_coord2 = {layer: 10000 for layer in tf_layers}

            if use_previously_matched_muons:

                matchedmu = None
                min_dr = 0.2
                trackpt = pt_to_gev_factor*ptLookup(abs(track_curvature))

                tpsmu=fetchTPS(event,'l1tTkMuonsGmt',2.5)
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
                        best_stubs_by_layer_coord1[stub.tfLayer()] = stub
                        stub_phi = int((1<<12)*stub.offline_coord1()/math.pi)
                        prop_phi = propagate(gmtphi, prop_eta, track_curvature, stub.tfLayer(), 'coord1', track_is_barrel[stub.tfLayer()])
                        phi_resids_by_layer_coord1[stub.tfLayer()] = deltaPhiINT(stub_phi, prop_phi)
                    if (stub.quality() & 0x2) != 0:
                        best_stubs_by_layer_coord2[stub.tfLayer()] = stub
                        stub_phi = int((1<<12)*stub.offline_coord2()/math.pi)
                        prop_phi = propagate(gmtphi, prop_eta, track_curvature, stub.tfLayer(), 'coord2', track_is_barrel[stub.tfLayer()])
                        phi_resids_by_layer_coord2[stub.tfLayer()] = deltaPhiINT(stub_phi, prop_phi)

            else:

                stubs = fetchStubs(event,'l1tStubsGmt:tps')

                for stub in stubs:

                    if (stub.quality() & 0x1) != 0:
                        stub_phi = int((1<<12)*stub.offline_coord1()/math.pi)
                        prop_phi = propagate(gmtphi, prop_eta, track_curvature, stub.tfLayer(), 'coord1', track_is_barrel[stub.tfLayer()])
                        phi_resid = deltaPhiINT(stub_phi, prop_phi)
                        float_stub_eta = stub.offline_eta1()
                        float_track_eta = track.eta()
                        float_eta_resid = float_stub_eta - float_track_eta
                        if (abs(phi_resid) < abs(phi_resids_by_layer_coord1[stub.tfLayer()])) and (abs(float_eta_resid) < 0.2):
                            phi_resids_by_layer_coord1[stub.tfLayer()] = phi_resid
                            best_stubs_by_layer_coord1[stub.tfLayer()] = stub

                    if (stub.quality() & 0x2) != 0:
                        if stub.isBarrel():
                            stub_phi = int((1<<12)*stub.offline_coord1()/math.pi)
                            prop_phi = propagate(gmtphi, prop_eta, track_curvature, stub.tfLayer(), 'coord1', track_is_barrel[stub.tfLayer()])
                        else:
                            stub_phi = int((1<<12)*stub.offline_coord2()/math.pi)
                            prop_phi = propagate(gmtphi, prop_eta, track_curvature, stub.tfLayer(), 'coord2', track_is_barrel[stub.tfLayer()])
                        phi_resid = deltaPhiINT(stub_phi, prop_phi)
                        float_stub_eta = stub.offline_eta1()
                        float_track_eta = track.eta()
                        float_eta_resid = float_stub_eta - float_track_eta
                        if (abs(phi_resid) < abs(phi_resids_by_layer_coord2[stub.tfLayer()])) and (abs(float_eta_resid) < 0.2):
                            phi_resids_by_layer_coord2[stub.tfLayer()] = phi_resid
                            best_stubs_by_layer_coord2[stub.tfLayer()] = stub

            # fill output histograms with results of matching
            for layer, stub in best_stubs_by_layer_coord1.items():
                if stub is None:
                    continue
                deltaphi = deltaPhiINT(int((1<<12)*stub.offline_coord1()/math.pi), gmtphi) if mode == 'prop' else phi_resids_by_layer_coord1[layer]
                calibrationH1[layer].Fill(prop_eta, track_curvature, deltaphi)
            for layer, stub in best_stubs_by_layer_coord2.items():
                if stub is None:
                    continue
                comp_phi = 0 if stub.isBarrel() else gmtphi
                deltaphi = deltaPhiINT(int((1<<12)*stub.offline_coord2()/math.pi), comp_phi) if mode == 'prop' else phi_resids_by_layer_coord2[layer]
                calibrationH2[layer].Fill(prop_eta, track_curvature, deltaphi)
                
                '''
                #eta res or prop, on the absolute value of eta
                track_eta = int(gmt_eta>>5)# if track.eta()>0 else int((-gmt_eta)>>5)
                if abs(track.eta())<0.83:
                    continue #just do endcap for now
                #track_eta = gmt_eta if track.eta()>0 else -gmt_eta
                if (s.etaQuality()&0x1)!=0:
                    #stub_eta = int((1<<12)*s.offline_eta1()/math.pi)
                    stub_eta = abs(s.eta1())#*32
                    #calibrationH1[s.tfLayer()].Fill(gmt_eta, stub_eta)
                    calibrationH1[s.tfLayer()].Fill(prop_eta,track_curvature,stub_eta-track_eta)
                    #if abs(stub_eta - track_eta) > 150 and abs(track_curvature) < 3000:
                    #    print(f'large delta eta found in layer {s.tfLayer()}')
                    #    print(f'muon: pt = {matchedmu.pt()}, phi = {matchedmu.phi()}, eta = {matchedmu.eta()}')
                    #    print(f'stubs:')
                    #    for sss in matchedmu.stubs():
                    #        print(f'    layer = {sss.tfLayer()}, c1 = {sss.offline_coord1()}, c2 = {sss.offline_coord2()}, eta1 = {sss.offline_eta1()}, eta2 = {sss.offline_eta2()}')
                    #    print('\n\n\n')
                        
                if (s.etaQuality()&0x2)!=0:
                    #stub_eta = int((1<<12)*s.offline_eta2()/math.pi)
                    stub_eta = abs(s.eta2())#*32
                    #calibrationH1[s.tfLayer()].Fill(gmt_eta, stub_eta)
                    calibrationH2[s.tfLayer()].Fill(prop_eta,track_curvature,stub_eta-track_eta)
                '''

    fname = f"calibrationResults_{mode}.root" if postfix == '' else f"calibrationResults_{mode}_{postfix}.root"
    f=ROOT.TFile(fname,"RECREATE")
    f.cd()
    for layer in tf_layers:
        calibrationH1[layer].Write()
        calibrationH2[layer].Write()
    f.Close()

if __name__ == "__main__":
    main()