import numpy as np
import ROOT
import os
import argparse
from common import barrel_limits, load_module_from_file, get_command_and_commit



def format_graph_command_info(input_graph_regions, barrel_file, endcap_file):

    barrel_commands_str = 'barrel_commands = {'
    barrel_commits_str = 'barrel_commits = {'
    endcap_commands_str = 'endcap_commands = {'
    endcap_commits_str = 'endcap_commits = {'

    for name, region in input_graph_regions.items():

        if input_graph_regions[name] != 'endcap':
            barrel_command = barrel_file.Get(f'command_{name}').GetTitle()
            barrel_commit = barrel_file.Get(f'commit_{name}').GetTitle()
        else:
            barrel_command = endcap_file.Get(f'command_{name}').GetTitle()
            barrel_commit = endcap_file.Get(f'commit_{name}').GetTitle()
        
        if input_graph_regions[name] != 'barrel':
            endcap_command = endcap_file.Get(f'command_{name}').GetTitle()
            endcap_commit = endcap_file.Get(f'commit_{name}').GetTitle()
        else:
            endcap_command = barrel_file.Get(f'command_{name}').GetTitle()
            endcap_commit = barrel_file.Get(f'commit_{name}').GetTitle()
        
        barrel_commands_str = barrel_commands_str + f"'{name}': '{barrel_command}',\n                   "
        endcap_commands_str = endcap_commands_str + f"'{name}': '{endcap_command}',\n                   "
        barrel_commits_str = barrel_commits_str + f"'{name}': '{barrel_commit}',\n                  "
        endcap_commits_str = endcap_commits_str + f"'{name}': '{endcap_commit}',\n                  "
    
    barrel_commands_str = barrel_commands_str + '}\n\n'
    endcap_commands_str = endcap_commands_str + '}\n\n'
    barrel_commits_str = barrel_commits_str + '}\n\n'
    endcap_commits_str = endcap_commits_str + '}\n\n'

    return barrel_commands_str, endcap_commands_str, barrel_commits_str, endcap_commits_str



def clean_and_interpolate_graphs(barrel_graph, endcap_graph, trash_eta_zones, sample_eta, barrel_limit):

    if barrel_graph is None:
        barrel_graph = endcap_graph
    if endcap_graph is None:
        endcap_graph = barrel_graph

    eta_width = barrel_graph.GetPointX(1) - barrel_graph.GetPointX(0) #assumes uniform eta_width

    barrel_graph = barrel_graph.Clone(f'{barrel_graph.GetName()}_clone')
    endcap_graph = endcap_graph.Clone(f'{endcap_graph.GetName()}_clone')

    #remove erroneous points, determined manually
    for eta_region, graph in zip(['barrel', 'endcap'], [barrel_graph, endcap_graph]):
        for i in np.flip(np.arange(graph.GetN())):
            eta = graph.GetPointX(int(i))
            if (eta_region == 'barrel' and eta >= barrel_limit) or (eta_region == 'endcap' and eta < barrel_limit):
                graph.RemovePoint(int(i))
                continue
            if trash_eta_zones != []:
                for limit_pair in trash_eta_zones:
                    if (limit_pair[0] <= eta and eta < limit_pair[1]):
                        graph.RemovePoint(int(i))
                        break
    
    lower_eta_graph = barrel_graph if barrel_graph.GetN() != 0 else endcap_graph
    upper_eta_graph = endcap_graph if endcap_graph.GetN() != 0 else barrel_graph

    min_eta_eval = min(lower_eta_graph.GetX()) - eta_width/2
    max_eta_eval = max(upper_eta_graph.GetX()) + eta_width/2

    y_min_eta_eval = lower_eta_graph.Eval(min_eta_eval)
    y_max_eta_eval = upper_eta_graph.Eval(max_eta_eval)

    ymin = min(min(lower_eta_graph.GetY()), min(upper_eta_graph.GetY()), y_min_eta_eval, y_max_eta_eval)
    ymax = max(max(lower_eta_graph.GetY()), max(upper_eta_graph.GetY()), y_min_eta_eval, y_max_eta_eval)

    table = np.zeros(len(sample_eta))
    for i, eta in enumerate(sample_eta):
        if eta < min_eta_eval:
            table[i] = y_min_eta_eval
        elif eta < barrel_limit:
            table[i] = barrel_graph.Eval(eta)
        elif eta < max_eta_eval:
            table[i] = endcap_graph.Eval(eta)
        else:
            table[i] = y_max_eta_eval

    interpo_graph = ROOT.TGraph()
    interpo_graph.SetName(lower_eta_graph.GetName())
    interpo_graph.SetTitle(lower_eta_graph.GetTitle())
    interpo_graph.GetXaxis().SetTitle(lower_eta_graph.GetXaxis().GetTitle())
    interpo_graph.GetYaxis().SetTitle(lower_eta_graph.GetYaxis().GetTitle()) 
    for (eta, interpo_value) in zip(sample_eta, table):
        interpo_graph.AddPoint(eta, interpo_value)

    barrel_graph.GetXaxis().SetRangeUser(min_eta_eval, max_eta_eval)
    barrel_graph.GetYaxis().SetRangeUser(ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin))
    endcap_graph.GetXaxis().SetRangeUser(min_eta_eval, max_eta_eval)
    endcap_graph.GetYaxis().SetRangeUser(ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin))
    interpo_graph.GetXaxis().SetRangeUser(min_eta_eval, max_eta_eval)
    interpo_graph.GetYaxis().SetRangeUser(ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin))

    return table, interpo_graph, barrel_graph, endcap_graph # return graphs for diagnostic plots and table for making LUT
    

    
def parseArgs(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['prop', 'res'], help='Whether we are handling propagation or resolution coefficients')
    parser.add_argument('--inputBarrelFile', help='Input root file containing fitted calibration parameter graphs for barrel')
    parser.add_argument('--inputEndcapFile', help='Input root file containing fitted calibration parameter graphs for endcap')
    parser.add_argument('--outdir', default='./output_LUTs/', help='Directory to write outputs in (default: ./output_LUts/)')
    parser.add_argument('--overwrite', action='store_true', help='overwrite existing LUTs that this would collide with')
    parser.add_argument('--config', help='python file containing config for interpolation and LUT format')
    args = parser.parse_args(argv)
    return args



def main():

    #make root not spam me
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
    ROOT.gStyle.SetOptFit(0)

    args = parseArgs()
    mode = args.mode
    barrel_file = ROOT.TFile.Open(args.inputBarrelFile)
    endcap_file = ROOT.TFile.Open(args.inputEndcapFile)
    outdir = args.outdir
    overwrite = args.overwrite
    os.makedirs(outdir, exist_ok=True)

    config = load_module_from_file(args.config)
    nbins = config.nbins
    num_cols = config.num_cols
    factor_bits = config.factor_bits
    lut_bits = config.lut_bits
    trash_eta_zones_by_layer = config.trash_eta_zones_by_layer

    outfile_txt_name = f'{outdir}/{mode}_LUTs.txt'
    outfile_py_name = f'{outdir}/{mode}_LUTs.py'
    if not args.overwrite and (os.path.isfile(outfile_txt_name) or os.path.isfile(outfile_py_name)):
        raise Exception('output files already exist and --overwrite was not specified')
    outfile_txt = open(outfile_txt_name, 'w')
    outfile_py = open(outfile_py_name, 'w')
    
    # save the command info
    command, commit = get_command_and_commit()
    outfile_py.write(f"command = '{command}'\n")
    outfile_py.write(f"commit = '{commit}'\n\n")
    outfile_txt.write(f'// command = {command}\n')
    outfile_txt.write(f'// commit = {commit}\n\n')

    assert 512 % nbins == 0
    step = 512 / nbins
    sample_eta = np.arange(step/2, 512, step) # midpoint sampling

    c = ROOT.TCanvas('canvas', 'canvastitle', 1200, 800)
    c.SetLeftMargin(0.15)

    print(f'Working in {mode} parameter mode')

    # load inputs and organize each parameter by if it has a barrel graph only, endcap graph only, or both
    barrel_graph_names = [key.GetName() for key in barrel_file.GetListOfKeys() if not key.GetName().startswith(('command', 'commit'))]
    endcap_graph_names = [key.GetName() for key in endcap_file.GetListOfKeys() if not key.GetName().startswith(('command', 'commit'))]

    input_graph_regions = {}
    for name in barrel_graph_names:
        if name in endcap_graph_names:
            input_graph_regions[name] = 'both'
        else:
            input_graph_regions[name] = 'barrel'
    for name in endcap_graph_names:
        if name in barrel_graph_names:
            continue
        else:
            input_graph_regions[name] = 'endcap'
    # sort it so that output files print LUTs in the right order
    input_graph_regions = dict(sorted(input_graph_regions.items()))

    # save the command info for each input graph at the top of the output files
    barrel_commands_str, endcap_commands_str, barrel_commits_str, endcap_commits_str = format_graph_command_info(input_graph_regions, barrel_file, endcap_file)
    outfile_py.write(barrel_commands_str)
    outfile_py.write(endcap_commands_str)
    outfile_py.write(barrel_commits_str)
    outfile_py.write(endcap_commits_str)
    outfile_py.write('\n\n')
    outfile_txt.write(barrel_commands_str)
    outfile_txt.write(endcap_commands_str)
    outfile_txt.write(barrel_commits_str)
    outfile_txt.write(endcap_commits_str)
    outfile_txt.write('\n\n')

    outfile_py.write(f'factor_bits = {str(factor_bits)}\n\n')
    outfile_txt.write(f'// factor_bits = {str(factor_bits)}\n\n')

    # now read from the graphs and build the LUTs
    for name, region in input_graph_regions.items():

        print(f'Processing {name}')

        param, layer = name.split('_')
        layer = int(layer)
        
        barrel_graph = barrel_file.Get(name) if input_graph_regions[name] != 'endcap' else None
        endcap_graph = endcap_file.Get(name) if input_graph_regions[name] != 'barrel' else None
        trash_eta_zones = trash_eta_zones_by_layer[layer]
        table, interpo_graph, barrel_graph, endcap_graph = clean_and_interpolate_graphs(barrel_graph, endcap_graph, trash_eta_zones, sample_eta, barrel_limits[layer])
        
        # save diagnostic plot
        interpo_graph.Draw('AL')
        barrel_graph.Draw('P')
        endcap_graph.Draw('P')
        c.Draw()
        c.SaveAs(f'{outdir}/interpo_{name}.png')

        # now make the LUTs
        factor = 2**(factor_bits[param])
        max_entry = 2**lut_bits[param] - 1
        table = 2**(factor_bits[param]) * table + 0.5
        table = table.astype(int)
        for index, entry in enumerate(table):
            if entry > max_entry:
                print(f'    saturating {param} from {entry} to {max_entry} at eta = {sample_eta[index]}')
                table[index] = max_entry
            if entry < 0:
                print(f'    overwriting {param} < 0 at eta = {sample_eta[index]}')
                table[index] = 0
        
        #headers for output file tables
        if mode == 'prop':
            header_txt = f'  const ap_uint<BITSPROPCOORD> lt_{name}[{nbins}] = ' + '{\n'
        elif mode == 'res':
            if param == 'res0':
                header_txt = f'  const ap_uint<BITSPROPSIGMACOORD_A> lt_{name}[{nbins}] = ' + '{\n'
            elif param == 'res1':
                header_txt = f'  const ap_uint<BITSPROPSIGMACOORD_B> lt_{name}[{nbins}] = ' + '{\n'

        #python output file as well for use in res propagation and for sanity checks
        header_py = f'lt_{name} = [' + '\n'
        
        #format the tables nicely and store them in a string
        num_digits = int(np.log10(max_entry)) + 1

        sub_outstring = '    '
        for j, val in enumerate(table):
            val = str(val) + ','
            sub_outstring = sub_outstring + f'{val:<{num_digits+2}}'
            if (j+1)%num_cols == 0:
                sub_outstring = sub_outstring + '\n    '

        comma_index = len(sub_outstring) - 1
        while sub_outstring[comma_index] != ',':
            comma_index = comma_index - 1
        sub_outstring = sub_outstring[:comma_index] 

        outstring_txt = header_txt + sub_outstring + '};\n\n'
        outstring_py = header_py + sub_outstring + ']\n\n'

        outfile_txt.write(outstring_txt)
        outfile_py.write(outstring_py)
        
        print('    done')

    outfile_txt.close()
    outfile_py.close()



if __name__ == "__main__":
    main()
