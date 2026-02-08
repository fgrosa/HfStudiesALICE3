#!/usr/bin/env python3

import argparse
import ROOT

def set_obj_style(obj, color, marker=ROOT.kFullCircle, fillstyle=-1, alpha=1.1):
    """
    Set style
    """
    if isinstance(obj, ROOT.TH1):
        obj.SetDirectory(0)
    obj.SetLineWidth(2)
    obj.SetLineColor(color)
    obj.SetMarkerColor(color)
    obj.SetMarkerStyle(marker)
    if fillstyle > -1:
        obj.SetFillStyle(fillstyle)
        obj.SetFillColorAlpha(color, alpha)


def compare(input_files, leg_labels, outfile_name):
    """
    Main function for efficiency comparison

    --------------------------------
    PARAMETERS
    input_files: list of input file names
    """

    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadTopMargin(0.035)
    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetTitleSize(0.045, "xy")
    ROOT.gStyle.SetLabelSize(0.045, "xy")
    ROOT.gStyle.SetTitleOffset(1.3, "x")
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    if len(leg_labels) != len(input_files):
        print("ERROR: number of legend lables does not match number of input files")
        for _ in range(len(leg_labels), len(input_files)):
            leg_labels.append("")

    colors = [ROOT.kRed-7, ROOT.kRed+2, ROOT.kAzure+2, ROOT.kAzure+4]
    markers = [ROOT.kOpenSquare, ROOT.kFullSquare, ROOT.kOpenCircle, ROOT.kFullCircle]

    graph_eff_vspt, graph_eff_vseta = [], []
    for ifile, infile_name in enumerate(input_files):
        infile = ROOT.TFile.Open(infile_name)
        graph_eff_vspt.append(infile.Get("trackeff_vs_pT"))
        graph_eff_vseta.append(infile.Get("trackeff_vs_eta"))
        set_obj_style(graph_eff_vspt[ifile], colors[ifile], markers[ifile])
        set_obj_style(graph_eff_vseta[ifile], colors[ifile], markers[ifile])

    canv = ROOT.TCanvas("canv", "canv", 1000, 500)
    leg = ROOT.TLegend(0.4, 0.2, 0.9, 0.5)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetMargin(0.1)
    for graph, label in zip(graph_eff_vspt, leg_labels):
        leg.AddEntry(graph, label, "lp")
    canv.Divide(2, 1)
    frame = canv.cd(1).DrawFrame(0.05, 1.e-3, 4., 1.2, ";#it{p}_{T} (GeV/#it{c}); tracking efficiency")
    canv.cd(1).SetLogx()
    frame.GetYaxis().SetDecimals()
    frame.GetXaxis().SetMoreLogLabels()
    for graph in graph_eff_vspt:
        graph.Draw("pzsame")
    leg.Draw()
    frame = canv.cd(2).DrawFrame(-4, 1.e-3, 4, 1.2, ";#it{#eta}; tracking efficiency")
    frame.GetYaxis().SetDecimals()
    for graph in graph_eff_vseta:
        graph.Draw("pzsame")
    canv.SaveAs(outfile_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--infiles", "-i", nargs="+", type=str,
                        default="", help="List of input file names")
    parser.add_argument("--labels", "-l", nargs="+", type=str,
                        default="", help="Legend labels")
    parser.add_argument("--outfile", "-o", metavar="text",
                        default="efficiency.pdf", help="Output file name")
    args = parser.parse_args()

    compare(args.infiles, args.labels, args.outfile)
